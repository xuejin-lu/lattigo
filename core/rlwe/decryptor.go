package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// Decryptor is a structure used to decrypt [Ciphertext]. It stores the secret-key.
type Decryptor struct {
	params                       Parameters
	ringQ                        *ring.Ring
	poolQ                        *ring.BufferPool
	sk                           *SecretKey
	fastCKKSZeroSecretSimulation bool
}

// NewDecryptor instantiates a new generic RLWE [Decryptor].
func NewDecryptor(params ParameterProvider, sk *SecretKey) *Decryptor {

	p := params.GetRLWEParameters()
	_, fastCKKSZeroSecretSimulation := params.(fastCKKSZeroSecretSimulationProvider)

	if sk.Value.Q.N() != p.N() {
		panic(fmt.Errorf("cannot NewDecryptor: secret_key ring degree does not match parameters ring degree"))
	}

	return &Decryptor{
		params:                       *p,
		ringQ:                        p.RingQ(),
		sk:                           sk,
		poolQ:                        ring.NewPool(p.ringQ),
		fastCKKSZeroSecretSimulation: fastCKKSZeroSecretSimulation,
	}
}

// GetRLWEParameters returns the underlying [Parameters].
func (d Decryptor) GetRLWEParameters() *Parameters {
	return &d.params
}

// DecryptNew decrypts the [Ciphertext] and returns the result in a new [Plaintext].
// Output pt [MetaData] will match the input ct [MetaData]. For a validated
// compact Fast CKKS zero-secret ciphertext, it returns an experimental
// plaintext projection at min(ct.Level(), 3); the input ciphertext's logical
// Level is unchanged and the result is not a full high-Level decryption.
func (d Decryptor) DecryptNew(ct *Ciphertext) (pt *Plaintext) {
	if d.fastCKKSZeroSecretSimulation {
		projectionLevel, useProjection, err := d.fastCKKSZeroSecretProjectionLevel(ct)
		if err != nil {
			panic(err)
		}
		if useProjection {
			pt = NewPlaintext(d.params, projectionLevel)
			d.decryptFastCKKSZeroSecretProjection(ct, pt, projectionLevel)
			return
		}
	}
	pt = NewPlaintext(d.params, ct.Level())
	d.Decrypt(ct, pt)
	return
}

// Decrypt decrypts the [Ciphertext] and writes the result in pt.
// The level of the output [Plaintext] is min(ct.Level(), pt.Level()), except
// that a validated compact Fast CKKS zero-secret ciphertext is projected to
// min(ct.Level(), 3) before applying the output-plaintext level limit. Output
// pt [MetaData] will match the input ct [MetaData].
func (d Decryptor) Decrypt(ct *Ciphertext, pt *Plaintext) {
	if d.fastCKKSZeroSecretSimulation {
		projectionLevel, useProjection, err := d.fastCKKSZeroSecretProjectionLevel(ct)
		if err != nil {
			panic(err)
		}
		if useProjection {
			d.decryptFastCKKSZeroSecretProjection(ct, pt, projectionLevel)
			return
		}
	}

	level := utils.Min(ct.Level(), pt.Level())

	ringQ := d.ringQ.AtLevel(level)

	pt.Resize(0, level)

	*pt.MetaData = *ct.MetaData

	if ct.IsNTT {
		pt.Value.CopyLvl(level, ct.Value[ct.Degree()])
	} else {
		ringQ.NTTLazy(ct.Value[ct.Degree()], pt.Value)
	}

	for i := ct.Degree(); i > 0; i-- {

		ringQ.MulCoeffsMontgomery(pt.Value, d.sk.Value.Q, pt.Value)

		if !ct.IsNTT {
			poolQ := d.poolQ.AtLevel(level)
			buff := poolQ.GetBuffPoly()
			ringQ.NTTLazy(ct.Value[i-1], *buff)
			ringQ.Add(pt.Value, *buff, pt.Value)
			poolQ.RecycleBuffPoly(buff)
		} else {
			ringQ.Add(pt.Value, ct.Value[i-1], pt.Value)
		}

		if i&7 == 7 {
			ringQ.Reduce(pt.Value, pt.Value)
		}
	}

	if (ct.Degree())&7 != 7 {
		ringQ.Reduce(pt.Value, pt.Value)
	}

	if !ct.IsNTT {
		ringQ.INTT(pt.Value, pt.Value)
	}
}

// fastCKKSZeroSecretProjectionLevel validates a compact Fast CKKS ciphertext
// and reports the low-level plaintext projection that can be decoded from its
// authoritative Q prefix. Complete-Q ciphertexts deliberately remain on the
// ordinary generic RLWE decryption path. This experimental projection does not
// claim to decrypt the ciphertext at its higher logical Level.
func (d Decryptor) fastCKKSZeroSecretProjectionLevel(ct *Ciphertext) (projectionLevel int, useProjection bool, err error) {
	if !d.fastCKKSZeroSecretSimulation {
		return 0, false, nil
	}
	if ct == nil {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: ciphertext cannot be nil")
	}
	if ct.MetaData == nil {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: ciphertext metadata cannot be nil")
	}
	if len(ct.Value) == 0 {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: ciphertext has no polynomial components")
	}
	if len(ct.Value[0].Coeffs) == 0 || len(ct.Value[0].Coeffs[0]) == 0 {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: c0 has no logical rows")
	}
	level := len(ct.Value[0].Coeffs) - 1
	if level < 0 || level > d.params.MaxLevel() {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: logical Level %d is outside configured Q", level)
	}
	if ct.Value[0].N() != d.params.N() {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: c0 ring degree %d does not match parameters N=%d", ct.Value[0].N(), d.params.N())
	}

	// Preserve the generic path whenever every logical Q row is fully backed.
	// This inspects only row-slice lengths above the Q-prefix; it never reads
	// dormant residue values.
	completeQBacking := true
	for component := range ct.Value {
		poly := ct.Value[component]
		if len(poly.Coeffs) != level+1 {
			return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: component %d has inconsistent logical Q rows", component)
		}
		if poly.N() != d.params.N() {
			return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: component %d ring degree does not match parameters N=%d", component, d.params.N())
		}
		for row := 0; row <= level; row++ {
			if len(poly.Coeffs[row]) != d.params.N() {
				completeQBacking = false
			}
		}
	}
	if completeQBacking {
		return 0, false, nil
	}

	if d.params.RingType() != ring.Standard {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: only the Standard CKKS ring is supported")
	}
	if ct.IsMontgomery {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: Montgomery-domain plaintext is unsupported by the CKKS Decoder")
	}
	if ct.Scale.Value.Sign() <= 0 || ct.Scale.Value.IsInf() || (ct.Scale.Mod != nil && ct.Scale.Mod.Sign() <= 0) {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: ciphertext Scale metadata is invalid")
	}
	if ct.LogDimensions.Rows != 0 || ct.LogDimensions.Cols < 0 || ct.LogDimensions.Cols > d.params.LogN()-1 {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: CKKS slot-dimension metadata is invalid")
	}
	rows := min(level+1, 4)
	if rows < 1 || rows > len(d.params.RingQ().SubRings) {
		return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: required Q-prefix width %d is unavailable", rows)
	}

	for component := range ct.Value {
		poly := ct.Value[component]
		for row := 0; row < rows; row++ {
			if len(poly.Coeffs[row]) != d.params.N() {
				return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: component %d is missing authoritative q%d backing", component, row)
			}
			modulus := d.params.RingQ().SubRings[row].Modulus
			for coefficient, residue := range poly.Coeffs[row] {
				if residue >= modulus {
					return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: component %d q%d coefficient %d is not a canonical residue", component, row, coefficient)
				}
				if component > 0 && residue != 0 {
					return 0, false, fmt.Errorf("cannot Fast CKKS decrypt projection: non-c0 component %d is nonzero at q%d coefficient %d", component, row, coefficient)
				}
			}
		}
	}

	return min(level, 3), true, nil
}

// decryptFastCKKSZeroSecretProjection copies only the already-validated c0
// prefix into a lower-Level plaintext. Validation is completed before the
// caller's output is resized or otherwise modified.
func (d Decryptor) decryptFastCKKSZeroSecretProjection(ct *Ciphertext, pt *Plaintext, projectionLevel int) {
	if pt == nil || pt.MetaData == nil || len(pt.Element.Value) != 1 || len(pt.Element.Value[0].Coeffs) == 0 || pt.Element.Value[0].N() != d.params.N() {
		panic(fmt.Errorf("cannot Fast CKKS decrypt projection: output plaintext structure is invalid"))
	}
	if pt.MetaData == ct.MetaData {
		panic(fmt.Errorf("cannot Fast CKKS decrypt projection: output plaintext metadata aliases input ciphertext"))
	}
	level := min(projectionLevel, pt.Level())
	rows := level + 1
	for row := 0; row < rows; row++ {
		if row < len(pt.Element.Value[0].Coeffs) {
			for component := range ct.Value {
				if row < len(ct.Value[component].Coeffs) && utils.Alias1D(pt.Element.Value[0].Coeffs[row], ct.Value[component].Coeffs[row]) {
					panic(fmt.Errorf("cannot Fast CKKS decrypt projection: output plaintext aliases ciphertext component %d q%d", component, row))
				}
			}
		}
	}
	pt.Resize(0, level)
	pt.Value = pt.Element.Value[0]
	*pt.MetaData = *ct.MetaData
	for row := 0; row < rows; row++ {
		copy(pt.Value.Coeffs[row], ct.Value[0].Coeffs[row])
	}
}

// WithKey creates a shallow copy of [Decryptor] with a new decryption key, in which all the
// data-structures are shared with the receiver.
// The receiver and the returned [Decryptor] can be used concurrently.
func (d Decryptor) WithKey(sk *SecretKey) *Decryptor {
	return &Decryptor{
		params:                       d.params,
		ringQ:                        d.ringQ,
		sk:                           sk,
		fastCKKSZeroSecretSimulation: d.fastCKKSZeroSecretSimulation,
	}
}
