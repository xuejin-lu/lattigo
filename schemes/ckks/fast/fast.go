package fast

import (
	"errors"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// RNSBackend provides a standalone Fast CKKS q0/q1 reconstruction
// and redistribution utility in the coefficient domain.
type RNSBackend struct {
	RingQ *ring.Ring
	Level int
}

// NewRNSBackend initializes a new RNSBackend.
func NewRNSBackend(ringQ *ring.Ring, level int) *RNSBackend {
	return &RNSBackend{
		RingQ: ringQ,
		Level: level,
	}
}

// ReconstructQ0Q1 reconstructs the signed coefficients from q0 and q1 residues.
// It explicitly requires coefficient-domain RNS polynomial input as an API precondition.
// It ONLY reads src.Coeffs[0] and src.Coeffs[1].
// It returns the unique centered representative modulo Q01.
func (backend *RNSBackend) ReconstructQ0Q1(src *ring.Poly) ([]*big.Int, error) {
	if backend.RingQ == nil {
		return nil, errors.New("RingQ cannot be nil")
	}
	if backend.Level < 1 || backend.RingQ.Level() < 1 {
		return nil, errors.New("requires at least q0 and q1 (Level >= 1)")
	}
	if backend.Level > backend.RingQ.Level() {
		return nil, errors.New("invalid level: exceeds RingQ.Level()")
	}
	if src == nil || src.Level() < 1 {
		return nil, errors.New("source polynomial requires at least q0 and q1")
	}
	N := backend.RingQ.N()
	if len(src.Coeffs[0]) < N || len(src.Coeffs[1]) < N {
		return nil, errors.New("source polynomial coefficient length is insufficient")
	}

	q0 := backend.RingQ.SubRings[0].Modulus
	q1 := backend.RingQ.SubRings[1].Modulus

	q0Big := new(big.Int).SetUint64(q0)
	q1Big := new(big.Int).SetUint64(q1)

	// Q01 = q0 * q1
	Q01 := new(big.Int).Mul(q0Big, q1Big)
	Q01Half := new(big.Int).Rsh(Q01, 1)

	// q0InvModq1 = q0^-1 mod q1
	// Using explicitly secure modulo inverse instead of Fermat's Little Theorem
	// because q1 primality is not formally guaranteed at this strict layer.
	q0InvModq1 := new(big.Int).ModInverse(q0Big, q1Big)
	if q0InvModq1 == nil {
		return nil, errors.New("modular inverse of q0 mod q1 does not exist")
	}

	res := make([]*big.Int, N)

	for k := 0; k < N; k++ {
		c0 := src.Coeffs[0][k]
		c1 := src.Coeffs[1][k]

		if c0 >= q0 || c1 >= q1 {
			return nil, errors.New("invalid residue: exceeds modulus bounds")
		}

		c0Big := new(big.Int).SetUint64(c0)
		c1Big := new(big.Int).SetUint64(c1)

		// u = (c1 - c0) * q0Inv mod q1
		u := new(big.Int).Sub(c1Big, c0Big)
		u.Mod(u, q1Big)
		u.Mul(u, q0InvModq1)
		u.Mod(u, q1Big)

		// x = c0 + u * q0
		x := new(big.Int).Mul(u, q0Big)
		x.Add(x, c0Big)

		// Center reconstruction around Q01 (-Q01/2 <= x < Q01/2)
		if x.Cmp(Q01Half) >= 0 {
			x.Sub(x, Q01)
		}

		res[k] = x
	}

	return res, nil
}

// Redistribute takes recovered signed coefficients and distributes them
// to all active RNS limbs (q0 ... qL) in the destination polynomial.
func (backend *RNSBackend) Redistribute(coeffs []*big.Int, dst *ring.Poly) error {
	if backend.RingQ == nil {
		return errors.New("RingQ cannot be nil")
	}
	if backend.Level < 0 || backend.Level > backend.RingQ.Level() {
		return errors.New("invalid level for redistribution")
	}

	N := backend.RingQ.N()
	if len(coeffs) < N {
		return errors.New("insufficient number of coefficients provided")
	}

	if dst == nil || dst.Level() < backend.Level {
		return errors.New("destination polynomial has insufficient RNS limbs")
	}

	for i := 0; i <= backend.Level; i++ {
		if len(dst.Coeffs[i]) < N {
			return errors.New("destination polynomial coefficient length is insufficient")
		}
	}

	for i := 0; i <= backend.Level; i++ {
		qi := backend.RingQ.SubRings[i].Modulus
		qiBig := new(big.Int).SetUint64(qi)

		for k := 0; k < N; k++ {
			tmp := new(big.Int).Mod(coeffs[k], qiBig)
			if tmp.Sign() < 0 {
				tmp.Add(tmp, qiBig)
			}
			dst.Coeffs[i][k] = tmp.Uint64()
		}
	}

	return nil
}

// ReconstructAndRedistribute performs both Phase 1A reconstruction from q0/q1
// and full RNS basis redistribution sequentially.
// It requires coefficient-domain RNS polynomial input.
func (backend *RNSBackend) ReconstructAndRedistribute(src *ring.Poly, dst *ring.Poly) error {
	coeffs, err := backend.ReconstructQ0Q1(src)
	if err != nil {
		return err
	}

	return backend.Redistribute(coeffs, dst)
}

// CheckCoefficientBounds explicitly checks if a reconstructed coefficient
// resides within a caller-provided absolute bound (-bound <= x <= bound).
func (backend *RNSBackend) CheckCoefficientBounds(x *big.Int, bound *big.Int) error {
	negBound := new(big.Int).Neg(bound)
	if x.Cmp(bound) > 0 || x.Cmp(negBound) < 0 {
		return errors.New("reconstructed coefficient is outside the provided bound")
	}
	return nil
}
