package fast

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// The Fast evaluator deliberately does not implement schemes.Evaluator. Its
// embedded rlwe.EvaluatorProvider requires DecomposeNTT,
// CheckAndGetGaloisKey, GadgetProductLazy, GadgetProductHoistedLazy,
// AutomorphismHoistedLazy, ModDownQPtoQNTT and AutomorphismIndex. Exposing
// Standard implementations of those methods would make full-Q/QP and
// key-switching execution available as a hidden fallback. This file is the
// smallest explicit surface needed by the current polynomial caller.

// Add adds an RLWE element or scalar to op0 in the NTT domain using q0/q1 only.
func (eval *Evaluator) Add(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) error {
	return eval.addSub(op0, op1, opOut, false)
}

func (eval *Evaluator) AddNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (*rlwe.Ciphertext, error) {
	degree, level := op0.Degree(), op0.Level()
	if el, ok := op1.(rlwe.ElementInterface[ring.Poly]); ok {
		degree, level = utils.Max(degree, el.Degree()), utils.Min(level, el.Level())
	}
	out := ckks.NewCiphertext(eval.Parameters, degree, level)
	return out, eval.Add(op0, op1, out)
}

// Sub subtracts an RLWE element or scalar from op0 in the NTT domain using q0/q1 only.
func (eval *Evaluator) Sub(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) error {
	return eval.addSub(op0, op1, opOut, true)
}

func (eval *Evaluator) SubNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (*rlwe.Ciphertext, error) {
	degree, level := op0.Degree(), op0.Level()
	if el, ok := op1.(rlwe.ElementInterface[ring.Poly]); ok {
		degree, level = utils.Max(degree, el.Degree()), utils.Min(level, el.Level())
	}
	out := ckks.NewCiphertext(eval.Parameters, degree, level)
	return out, eval.Sub(op0, op1, out)
}

func (eval *Evaluator) addSub(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext, sub bool) error {
	if err := eval.validateUnary(op0, opOut); err != nil {
		return err
	}
	if el, ok := op1.(rlwe.ElementInterface[ring.Poly]); ok {
		return eval.addSubElement(op0, el.El(), opOut, sub)
	}
	c, err := eval.scalar(op1)
	if err != nil {
		return err
	}
	level := utils.Min(op0.Level(), opOut.Level())
	opOut.Resize(op0.Degree(), level)
	*opOut.MetaData = *op0.MetaData
	if op0 != opOut {
		for d := range op0.Value {
			copyQ01(op0.Value[d], opOut.Value[d])
		}
	}
	values := eval.scalarNTT(c, &op0.Scale.Value, op0.IsMontgomery)
	for limb := 0; limb < 2; limb++ {
		s := eval.Parameters.RingQ().SubRings[limb]
		half := eval.Parameters.N() >> 1
		if sub {
			s.SubScalar(opOut.Value[0].Coeffs[limb][:half], values[limb][0], opOut.Value[0].Coeffs[limb][:half])
			s.SubScalar(opOut.Value[0].Coeffs[limb][half:], values[limb][1], opOut.Value[0].Coeffs[limb][half:])
		} else {
			s.AddScalar(opOut.Value[0].Coeffs[limb][:half], values[limb][0], opOut.Value[0].Coeffs[limb][:half])
			s.AddScalar(opOut.Value[0].Coeffs[limb][half:], values[limb][1], opOut.Value[0].Coeffs[limb][half:])
		}
	}
	return nil
}

func (eval *Evaluator) addSubElement(op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], opOut *rlwe.Ciphertext, sub bool) error {
	if err := eval.validateBinary(op0.El(), op1, opOut); err != nil {
		return err
	}
	if !op0.Scale.Equal(op1.Scale) {
		return errors.New("Fast Add/Sub requires equal operand scales")
	}
	level := utils.Min(utils.Min(op0.Level(), op1.Level()), opOut.Level())
	maxDegree, minDegree := utils.Max(op0.Degree(), op1.Degree()), utils.Min(op0.Degree(), op1.Degree())
	opOut.Resize(maxDegree, level)
	*opOut.MetaData = *op0.MetaData
	opOut.LogDimensions.Rows = utils.Max(op0.LogDimensions.Rows, op1.LogDimensions.Rows)
	opOut.LogDimensions.Cols = utils.Max(op0.LogDimensions.Cols, op1.LogDimensions.Cols)
	for d := 0; d <= minDegree; d++ {
		for limb := 0; limb < 2; limb++ {
			s := eval.Parameters.RingQ().SubRings[limb]
			if sub {
				s.Sub(op0.Value[d].Coeffs[limb], op1.Value[d].Coeffs[limb], opOut.Value[d].Coeffs[limb])
			} else {
				s.Add(op0.Value[d].Coeffs[limb], op1.Value[d].Coeffs[limb], opOut.Value[d].Coeffs[limb])
			}
		}
	}
	if op0.Degree() > minDegree && opOut != op0 {
		for d := minDegree + 1; d <= op0.Degree(); d++ {
			copyQ01(op0.Value[d], opOut.Value[d])
		}
	} else if op1.Degree() > minDegree {
		for d := minDegree + 1; d <= op1.Degree(); d++ {
			if sub {
				negQ01(eval.Parameters.RingQ(), op1.Value[d], opOut.Value[d])
			} else if opOut.El() != op1 {
				copyQ01(op1.Value[d], opOut.Value[d])
			}
		}
	}
	return nil
}

// Mul multiplies in NTT representation. Degree-one ciphertext products retain c2.
func (eval *Evaluator) Mul(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) error {
	if err := eval.validateUnary(op0, opOut); err != nil {
		return err
	}
	if el, ok := op1.(rlwe.ElementInterface[ring.Poly]); ok {
		return eval.mulElement(op0, el.El(), opOut, false)
	}
	c, err := eval.scalar(op1)
	if err != nil {
		return err
	}
	return eval.mulScalar(op0, c, opOut)
}

func (eval *Evaluator) MulNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (*rlwe.Ciphertext, error) {
	degree, level := op0.Degree(), op0.Level()
	if el, ok := op1.(rlwe.ElementInterface[ring.Poly]); ok {
		degree, level = op0.Degree()+el.Degree(), utils.Min(level, el.Level())
	}
	out := ckks.NewCiphertext(eval.Parameters, degree, level)
	return out, eval.Mul(op0, op1, out)
}

// MulRelin omits z2 because current zero-secret relinearization discards it.
func (eval *Evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) error {
	if err := eval.validateUnary(op0, opOut); err != nil {
		return err
	}
	if el, ok := op1.(rlwe.ElementInterface[ring.Poly]); ok {
		return eval.mulElement(op0, el.El(), opOut, true)
	}
	return eval.Mul(op0, op1, opOut)
}

func (eval *Evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (*rlwe.Ciphertext, error) {
	level := op0.Level()
	if el, ok := op1.(rlwe.ElementInterface[ring.Poly]); ok {
		level = utils.Min(level, el.Level())
	}
	out := ckks.NewCiphertext(eval.Parameters, 1, level)
	return out, eval.MulRelin(op0, op1, out)
}

func (eval *Evaluator) mulElement(op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], opOut *rlwe.Ciphertext, relin bool) error {
	if err := eval.validateBinary(op0.El(), op1, opOut); err != nil {
		return err
	}
	if op0.Degree() > 1 || op1.Degree() > 1 {
		return errors.New("Fast Mul requires operands of degree at most one")
	}
	degree := op0.Degree() + op1.Degree()
	if relin && degree == 2 {
		degree = 1
	}
	level := utils.Min(utils.Min(op0.Level(), op1.Level()), opOut.Level())
	if op0.Degree() == 1 && op1.Degree() == 1 {
		t0, t1, t2 := eval.nttScratch[0], eval.nttScratch[1], eval.nttScratch[2]
		eval.pointMul(op0.Value[0], op1.Value[0], t0, op0.IsMontgomery)
		if op0.El() == op1 {
			eval.pointMul(op0.Value[0], op0.Value[1], t1, op0.IsMontgomery)
			for limb := 0; limb < 2; limb++ {
				eval.Parameters.RingQ().SubRings[limb].Add(t1.Coeffs[limb], t1.Coeffs[limb], t1.Coeffs[limb])
			}
		} else {
			eval.pointMul(op0.Value[0], op1.Value[1], t1, op0.IsMontgomery)
			eval.pointMulThenAdd(op0.Value[1], op1.Value[0], t1, op0.IsMontgomery)
		}
		if !relin {
			eval.pointMul(op0.Value[1], op1.Value[1], t2, op0.IsMontgomery)
		}
		opOut.Resize(degree, level)
		copyQ01(t0, opOut.Value[0])
		copyQ01(t1, opOut.Value[1])
		if !relin {
			copyQ01(t2, opOut.Value[2])
		}
	} else {
		var ct []ring.Poly
		var pt ring.Poly
		if op0.Degree() == 0 {
			pt, ct = op0.Value[0], op1.Value
		} else {
			pt, ct = op1.Value[0], op0.Value
		}
		for d := range ct {
			eval.pointMul(pt, ct[d], eval.nttScratch[d], op0.IsMontgomery)
		}
		opOut.Resize(degree, level)
		for d := range ct {
			copyQ01(eval.nttScratch[d], opOut.Value[d])
		}
	}
	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale.Mul(op1.Scale)
	opOut.LogDimensions.Rows = utils.Max(op0.LogDimensions.Rows, op1.LogDimensions.Rows)
	opOut.LogDimensions.Cols = utils.Max(op0.LogDimensions.Cols, op1.LogDimensions.Cols)
	return nil
}

// Relinearize applies the current Stage-A zero-secret c2 truncation.
func (eval *Evaluator) Relinearize(op0, opOut *rlwe.Ciphertext) error {
	if err := eval.validateUnary(op0, opOut); err != nil {
		return err
	}
	if op0.Level() != opOut.Level() {
		return errors.New("Fast Relinearize requires equal input/output levels")
	}
	return FastTruncateDegree2To1(op0, opOut)
}

// MulThenAdd adds the product directly into q0/q1 output storage.
func (eval *Evaluator) MulThenAdd(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) error {
	if err := eval.validateUnary(op0, opOut); err != nil {
		return err
	}
	if el, ok := op1.(rlwe.ElementInterface[ring.Poly]); ok {
		if err := eval.validateBinary(op0.El(), el.El(), opOut); err != nil {
			return err
		}
		if op0 == opOut || el.El() == opOut.El() {
			return errors.New("Fast MulThenAdd RLWE output must not alias an input")
		}
		if op0.Degree() > 1 || el.Degree() > 1 || opOut.Degree() < op0.Degree()+el.Degree() {
			return errors.New("invalid Fast MulThenAdd RLWE degrees")
		}
		if !opOut.Scale.Equal(op0.Scale.Mul(el.El().Scale)) {
			return errors.New("Fast MulThenAdd RLWE requires output scale equal to product scale")
		}
		opOut.Resize(opOut.Degree(), utils.Min(utils.Min(op0.Level(), el.Level()), opOut.Level()))
		return eval.mulElementThenAdd(op0, el.El(), opOut)
	}
	c, err := eval.scalar(op1)
	if err != nil {
		return err
	}
	if opOut.Degree() != op0.Degree() {
		return errors.New("Fast scalar MulThenAdd requires matching input/output degrees")
	}
	level := utils.Min(op0.Level(), opOut.Level())
	source := op0.Value
	var scale rlwe.Scale
	if cmp := op0.Scale.Cmp(opOut.Scale); cmp == 0 {
		if c.IsInt() {
			scale = rlwe.NewScale(1)
		} else {
			// Preserve op0 when the accumulator aliases it: scaling opOut below
			// must not also change the multiplicand used by the fused add.
			if op0 == opOut {
				for d := range op0.Value {
					copyQ01(op0.Value[d], eval.nttScratch[d])
				}
				source = eval.nttScratch[:len(op0.Value)]
			}
			scale, err = eval.coefficientScale(level)
			if err != nil {
				return err
			}
			if err = eval.mulScalarAtScale(opOut, bignum.ToComplex(scale.BigInt(), eval.Parameters.EncodingPrecision()), rlwe.NewScale(1), opOut); err != nil {
				return err
			}
			opOut.Scale = opOut.Scale.Mul(scale)
		}
	} else if cmp < 0 {
		scale = opOut.Scale.Div(op0.Scale)
	} else {
		// Promote the accumulator to the term scale before the fused add.
		// Mod1's Standard target-scale schedule can produce a higher-scale
		// power term; multiplying the maintained residues and updating the
		// metadata preserves the represented value without a level transition.
		ratio := op0.Scale.Div(opOut.Scale).BigInt()
		if ratio.Sign() <= 0 {
			return errors.New("invalid Fast MulThenAdd scale promotion ratio")
		}
		if err := eval.MulIntegerMaintained(opOut, ratio, opOut); err != nil {
			return err
		}
		opOut.Scale = op0.Scale
		if c.IsInt() {
			scale = rlwe.NewScale(1)
		} else {
			if op0 == opOut {
				for d := range op0.Value {
					copyQ01(op0.Value[d], eval.nttScratch[d])
				}
				source = eval.nttScratch[:len(op0.Value)]
			}
			scale, err = eval.coefficientScale(level)
			if err != nil {
				return err
			}
			if err = eval.mulScalarAtScale(opOut, bignum.ToComplex(scale.BigInt(), eval.Parameters.EncodingPrecision()), rlwe.NewScale(1), opOut); err != nil {
				return err
			}
			opOut.Scale = opOut.Scale.Mul(scale)
		}
	}
	values := eval.scalarNTT(c, &scale.Value, false)
	for d := range source {
		for limb := 0; limb < 2; limb++ {
			s := eval.Parameters.RingQ().SubRings[limb]
			half := eval.Parameters.N() >> 1
			s.MulScalarMontgomeryThenAdd(source[d].Coeffs[limb][:half], ring.MForm(values[limb][0], s.Modulus, s.BRedConstant), opOut.Value[d].Coeffs[limb][:half])
			s.MulScalarMontgomeryThenAdd(source[d].Coeffs[limb][half:], ring.MForm(values[limb][1], s.Modulus, s.BRedConstant), opOut.Value[d].Coeffs[limb][half:])
		}
	}
	return nil
}

func (eval *Evaluator) mulElementThenAdd(op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], out *rlwe.Ciphertext) error {
	if op0.Degree() == 1 && op1.Degree() == 1 {
		eval.pointMulThenAdd(op0.Value[0], op1.Value[0], out.Value[0], op0.IsMontgomery)
		eval.pointMul(op0.Value[0], op1.Value[1], eval.nttScratch[0], op0.IsMontgomery)
		eval.pointMulThenAdd(op0.Value[1], op1.Value[0], eval.nttScratch[0], op0.IsMontgomery)
		for limb := 0; limb < 2; limb++ {
			eval.Parameters.RingQ().SubRings[limb].Add(out.Value[1].Coeffs[limb], eval.nttScratch[0].Coeffs[limb], out.Value[1].Coeffs[limb])
		}
		eval.pointMulThenAdd(op0.Value[1], op1.Value[1], out.Value[2], op0.IsMontgomery)
	} else {
		var ct []ring.Poly
		var pt ring.Poly
		if op0.Degree() == 0 {
			pt, ct = op0.Value[0], op1.Value
		} else {
			pt, ct = op1.Value[0], op0.Value
		}
		for d := range ct {
			eval.pointMulThenAdd(pt, ct[d], out.Value[d], op0.IsMontgomery)
		}
	}
	return nil
}

func (eval *Evaluator) mulScalar(op0 *rlwe.Ciphertext, c *bignum.Complex, out *rlwe.Ciphertext) error {
	scale := rlwe.NewScale(1)
	if !c.IsInt() {
		var err error
		scale, err = eval.coefficientScale(utils.Min(op0.Level(), out.Level()))
		if err != nil {
			return err
		}
	}
	return eval.mulScalarAtScale(op0, c, scale, out)
}

func (eval *Evaluator) mulScalarAtScale(op0 *rlwe.Ciphertext, c *bignum.Complex, scale rlwe.Scale, out *rlwe.Ciphertext) error {
	level := utils.Min(op0.Level(), out.Level())
	values := eval.scalarNTT(c, &scale.Value, false)
	out.Resize(op0.Degree(), level)
	for d := range op0.Value {
		for limb := 0; limb < 2; limb++ {
			s := eval.Parameters.RingQ().SubRings[limb]
			half := eval.Parameters.N() >> 1
			s.MulScalarMontgomery(op0.Value[d].Coeffs[limb][:half], ring.MForm(values[limb][0], s.Modulus, s.BRedConstant), out.Value[d].Coeffs[limb][:half])
			s.MulScalarMontgomery(op0.Value[d].Coeffs[limb][half:], ring.MForm(values[limb][1], s.Modulus, s.BRedConstant), out.Value[d].Coeffs[limb][half:])
		}
	}
	*out.MetaData = *op0.MetaData
	out.Scale = op0.Scale.Mul(scale)
	return nil
}

func (eval *Evaluator) pointMul(a, b, out ring.Poly, montgomery bool) {
	for limb := 0; limb < 2; limb++ {
		s := eval.Parameters.RingQ().SubRings[limb]
		if montgomery {
			s.MulCoeffsMontgomery(a.Coeffs[limb], b.Coeffs[limb], out.Coeffs[limb])
		} else {
			s.MulCoeffsBarrett(a.Coeffs[limb], b.Coeffs[limb], out.Coeffs[limb])
		}
	}
}

func (eval *Evaluator) pointMulThenAdd(a, b, out ring.Poly, montgomery bool) {
	for limb := 0; limb < 2; limb++ {
		s := eval.Parameters.RingQ().SubRings[limb]
		if montgomery {
			s.MulCoeffsMontgomeryThenAdd(a.Coeffs[limb], b.Coeffs[limb], out.Coeffs[limb])
		} else {
			s.MulCoeffsBarrettThenAdd(a.Coeffs[limb], b.Coeffs[limb], out.Coeffs[limb])
		}
	}
}

func (eval *Evaluator) validateUnary(op0, out *rlwe.Ciphertext) error {
	if eval == nil || op0 == nil || out == nil {
		return errors.New("Fast evaluator and operands cannot be nil")
	}
	if op0.MetaData == nil || out.MetaData == nil {
		return errors.New("Fast operand metadata cannot be nil")
	}
	if op0.N() != eval.Parameters.N() || out.N() != eval.Parameters.N() {
		return errors.New("Fast operand dimensions do not match parameters")
	}
	if op0.Level() < 1 || out.Level() < 1 {
		return errors.New("Fast evaluator requires q0 and q1")
	}
	if !op0.IsNTT || !out.IsNTT {
		return errors.New("Fast evaluator arithmetic requires NTT-domain operands")
	}
	if op0.IsMontgomery != out.IsMontgomery {
		return errors.New("Fast evaluator requires matching Montgomery representations")
	}
	return nil
}

func (eval *Evaluator) validateBinary(op0, op1 *rlwe.Element[ring.Poly], out *rlwe.Ciphertext) error {
	if op1 == nil || op1.MetaData == nil {
		return errors.New("Fast second operand and metadata cannot be nil")
	}
	if op1.N() != eval.Parameters.N() || op1.Level() < 1 {
		return errors.New("Fast second operand has invalid dimensions or level")
	}
	if !op1.IsNTT || op0.IsNTT != op1.IsNTT || op0.IsMontgomery != op1.IsMontgomery || op0.IsMontgomery != out.IsMontgomery {
		return errors.New("Fast operands require matching NTT/Montgomery representations")
	}
	if op0.IsBatched != op1.IsBatched {
		return errors.New("Fast operands require matching batching metadata")
	}
	return nil
}

func (eval *Evaluator) scalar(value rlwe.Operand) (*bignum.Complex, error) {
	if v, ok := value.(uint); ok {
		value = uint64(v)
	}
	switch value.(type) {
	case complex128, float64, int, int64, uint64, *big.Int, *big.Float, *bignum.Complex:
		return bignum.ToComplex(value, eval.Parameters.EncodingPrecision()), nil
	default:
		return nil, fmt.Errorf("unsupported Fast scalar type %T", value)
	}
}

func (eval *Evaluator) coefficientScale(level int) (rlwe.Scale, error) {
	consumed := eval.Parameters.LevelsConsumedPerRescaling()
	if level-consumed+1 < 0 {
		return rlwe.Scale{}, errors.New("insufficient level for Fast scalar encoding scale")
	}
	scale := rlwe.NewScale(1)
	for i := 0; i < consumed; i++ {
		scale = scale.Mul(rlwe.NewScale(eval.Parameters.RingQ().SubRings[level-i].Modulus))
	}
	return scale, nil
}

func (eval *Evaluator) scalarNTT(c *bignum.Complex, scale *big.Float, montgomery bool) [2][2]uint64 {
	var out [2][2]uint64
	toInt := func(x *big.Float) *big.Int {
		z := new(big.Int)
		if x == nil {
			return z
		}
		v := new(big.Float).Mul(x, scale)
		if x.Sign() > 0 {
			v.Add(v, new(big.Float).SetFloat64(0.5))
		} else if x.Sign() < 0 {
			v.Sub(v, new(big.Float).SetFloat64(0.5))
		}
		v.Int(z)
		return z
	}
	real, imag := toInt(c[0]), toInt(c[1])
	for limb := 0; limb < 2; limb++ {
		s := eval.Parameters.RingQ().SubRings[limb]
		r := new(big.Int).Mod(new(big.Int).Set(real), new(big.Int).SetUint64(s.Modulus)).Uint64()
		i := new(big.Int).Mod(new(big.Int).Set(imag), new(big.Int).SetUint64(s.Modulus)).Uint64()
		i = ring.MRed(i, s.RootsForward[1], s.Modulus, s.MRedConstant)
		out[limb][0], out[limb][1] = ring.CRed(r+i, s.Modulus), ring.CRed(r+s.Modulus-i, s.Modulus)
		if montgomery {
			out[limb][0] = ring.MForm(out[limb][0], s.Modulus, s.BRedConstant)
			out[limb][1] = ring.MForm(out[limb][1], s.Modulus, s.BRedConstant)
		}
	}
	return out
}
