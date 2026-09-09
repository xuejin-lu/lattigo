package polynomial

import (
	"errors"
	"fmt"
	"math/bits"

	commonpolynomial "github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// FastEvaluator evaluates a single Chebyshev polynomial with the explicit
// q0/q1-authoritative Fast CKKS evaluator. It is intentionally a bounded
// polynomial surface and does not implement schemes.Evaluator.
//
// FastEvaluator is not safe for concurrent use. Its workspace is reused by
// successive calls to Evaluate.
type FastEvaluator struct {
	Parameters ckks.Parameters
	Evaluator  *fastckks.Evaluator
	workspace  fastPolynomialWorkspace
}

// NewFastEvaluator returns a Fast polynomial evaluator. If eval is nil, a
// Fast CKKS evaluator is created for params.
func NewFastEvaluator(params ckks.Parameters, eval *fastckks.Evaluator) *FastEvaluator {
	if eval == nil {
		eval = fastckks.NewEvaluator(params)
	}
	return &FastEvaluator{
		Parameters: params,
		Evaluator:  eval,
		workspace: fastPolynomialWorkspace{
			powerBuffers: make(map[int]*rlwe.Ciphertext),
			powers:       make(map[int]*rlwe.Ciphertext),
		},
	}
}

// Evaluate evaluates the single polynomial p on input. The current Fast
// polynomial surface accepts Chebyshev polynomials in the Standard ring with
// NTT/Montgomery degree-one input at a level containing q0 and q1.
func (eval *FastEvaluator) Evaluate(input *rlwe.Ciphertext, p bignum.Polynomial, targetScale rlwe.Scale) (*rlwe.Ciphertext, error) {
	if eval == nil || eval.Evaluator == nil {
		return nil, errors.New("Fast polynomial evaluator cannot be nil")
	}
	if input == nil {
		return nil, errors.New("Fast polynomial input cannot be nil")
	}
	if input.MetaData == nil {
		return nil, errors.New("Fast polynomial input metadata cannot be nil")
	}
	if len(p.Coeffs) == 0 {
		return nil, errors.New("Fast polynomial cannot be empty")
	}
	if p.Basis != bignum.Chebyshev {
		return nil, fmt.Errorf("Fast polynomial evaluator only supports the Chebyshev basis, got %v", p.Basis)
	}
	if eval.Parameters.RingType() != ring.Standard {
		return nil, fmt.Errorf("Fast polynomial evaluator requires the Standard ring, got %s", eval.Parameters.RingType())
	}
	if input.N() != eval.Parameters.N() {
		return nil, errors.New("Fast polynomial input dimensions do not match parameters")
	}
	if input.Degree() != 1 {
		return nil, fmt.Errorf("Fast polynomial evaluator requires degree-one input, got degree %d", input.Degree())
	}
	if input.Level() < 1 {
		return nil, errors.New("Fast polynomial evaluator requires input level containing q0 and q1")
	}
	if !input.IsNTT {
		return nil, errors.New("Fast polynomial evaluator requires NTT-domain input")
	}
	if !input.IsMontgomery {
		return nil, errors.New("Fast polynomial evaluator requires Montgomery input")
	}
	if targetScale.Value.Sign() <= 0 {
		return nil, errors.New("Fast polynomial target scale must be positive")
	}

	levelsConsumed := eval.Parameters.LevelsConsumedPerRescaling()
	commonPoly := commonpolynomial.NewPolynomial(p)
	if input.Level() < levelsConsumed*commonPoly.Depth() {
		return nil, fmt.Errorf("%d levels < %d log(d) -> cannot evaluate poly", input.Level(), levelsConsumed*commonPoly.Depth())
	}

	ws := &eval.workspace
	ws.reset(eval.Parameters, input)

	// Degree zero has no power-generation or PS planning work. This also keeps
	// the bounded Fast surface well-defined for constant Chebyshev polynomials.
	if commonPoly.Degree() == 0 {
		if p.Coeffs[0] == nil {
			return nil, errors.New("Fast polynomial constant coefficient cannot be nil")
		}
		out := ws.babyBuffer(eval.Parameters, 0, 1, input.Level())
		zeroMaintained(out)
		*out.MetaData = *input.MetaData
		out.Scale = targetScale
		if err := eval.Evaluator.Add(out, p.Coeffs[0], out); err != nil {
			return nil, fmt.Errorf("Fast polynomial constant: %w", err)
		}
		return cloneMaintainedResult(eval.Parameters, out), nil
	}

	if err := ws.generatePowers(eval.Parameters, eval.Evaluator, p, commonPoly); err != nil {
		return nil, err
	}

	// Reuse the common PS planner and simulation so that Fast and Standard
	// polynomial evaluation retain the same decomposition and metadata plan.
	plan := commonPoly.PatersonStockmeyerPolynomial(
		eval.Evaluator.GetRLWEParameters(),
		ws.powers[1].Level(),
		ws.powers[1].Scale,
		targetScale,
		simEvaluator{params: eval.Parameters, levelsConsumedPerRescaling: levelsConsumed},
	)
	result, err := ws.evaluatePlan(eval.Parameters, eval.Evaluator, plan, ws.powers)
	if err != nil {
		return nil, err
	}
	return cloneMaintainedResult(eval.Parameters, result), nil
}

type fastPolynomialWorkspace struct {
	x1           *rlwe.Ciphertext
	powers       map[int]*rlwe.Ciphertext
	powerBuffers map[int]*rlwe.Ciphertext
	babyBuffers  []*rlwe.Ciphertext
	babySteps    []*fastBabyStep
	giantSteps   []int
	scaleScratch *rlwe.Ciphertext
}

type fastBabyStep struct {
	Degree int
	Value  *rlwe.Ciphertext
}

func (ws *fastPolynomialWorkspace) reset(params ckks.Parameters, input *rlwe.Ciphertext) {
	ws.x1 = ws.ensureCiphertext(params, ws.x1, 1, input.Level())
	copyMaintained(input, ws.x1)
	for key := range ws.powers {
		delete(ws.powers, key)
	}
	ws.powers[1] = ws.x1
}

func (ws *fastPolynomialWorkspace) ensureCiphertext(params ckks.Parameters, ct *rlwe.Ciphertext, degree, level int) *rlwe.Ciphertext {
	if ct == nil {
		ct = ckks.NewCiphertext(params, degree, level)
	} else {
		ct.Resize(degree, level)
	}
	return ct
}

func (ws *fastPolynomialWorkspace) powerBuffer(params ckks.Parameters, n, degree, level int, input *rlwe.Ciphertext) *rlwe.Ciphertext {
	ct := ws.powerBuffers[n]
	ct = ws.ensureCiphertext(params, ct, degree, level)
	ct.IsNTT = input.IsNTT
	ct.IsMontgomery = input.IsMontgomery
	*ct.MetaData = *input.MetaData
	ws.powerBuffers[n] = ct
	return ct
}

func (ws *fastPolynomialWorkspace) babyBuffer(params ckks.Parameters, index, degree, level int) *rlwe.Ciphertext {
	for len(ws.babyBuffers) <= index {
		ws.babyBuffers = append(ws.babyBuffers, nil)
	}
	ws.babyBuffers[index] = ws.ensureCiphertext(params, ws.babyBuffers[index], degree, level)
	return ws.babyBuffers[index]
}

func (ws *fastPolynomialWorkspace) generatePowers(params ckks.Parameters, eval *fastckks.Evaluator, p bignum.Polynomial, commonPoly commonpolynomial.Polynomial) error {
	logDegree := bits.Len64(uint64(commonPoly.Degree()))
	pb := fastPowerBasis{basis: p.Basis, values: ws.powers, workspace: ws, params: params, eval: eval}
	if err := pb.genPower(1<<(logDegree-1), false); err != nil {
		return err
	}
	logSplit := bignum.OptimalSplit(logDegree)
	for i := (1 << logSplit) - 1; i > 2; i-- {
		if !(p.IsEven || p.IsOdd) || (i&1 == 0 && p.IsEven) || (i&1 == 1 && p.IsOdd) {
			if err := pb.genPower(i, commonPoly.Lazy); err != nil {
				return err
			}
		}
	}
	return nil
}

type fastPowerBasis struct {
	basis     bignum.Basis
	values    map[int]*rlwe.Ciphertext
	workspace *fastPolynomialWorkspace
	params    ckks.Parameters
	eval      *fastckks.Evaluator
}

func (pb *fastPowerBasis) genPower(n int, lazy bool) error {
	if pb.values[n] != nil {
		return nil
	}
	if _, err := pb.genPowerInternal(n, lazy, true); err != nil {
		return err
	}
	if err := pb.eval.Rescale(pb.values[n], pb.values[n]); err != nil {
		return fmt.Errorf("Fast power %d: rescale: %w", n, err)
	}
	return nil
}

func (pb *fastPowerBasis) genPowerInternal(n int, lazy, rescale bool) (bool, error) {
	if pb.values[n] != nil {
		return false, nil
	}
	a, b := commonpolynomial.SplitDegree(n)
	isPow2 := n&(n-1) == 0

	rescaleA, err := pb.genPowerInternal(a, lazy && !isPow2, rescale)
	if err != nil {
		return false, fmt.Errorf("Fast power %d: power %d: %w", n, a, err)
	}
	rescaleB, err := pb.genPowerInternal(b, lazy && !isPow2, rescale)
	if err != nil {
		return false, fmt.Errorf("Fast power %d: power %d: %w", n, b, err)
	}

	left, right := pb.values[a], pb.values[b]
	degree := 1
	if lazy {
		if left.Degree() == 2 {
			if err := pb.eval.Relinearize(left, left); err != nil {
				return false, fmt.Errorf("Fast power %d: relinearize left: %w", n, err)
			}
		}
		if right.Degree() == 2 {
			if err := pb.eval.Relinearize(right, right); err != nil {
				return false, fmt.Errorf("Fast power %d: relinearize right: %w", n, err)
			}
		}
		if rescaleA {
			if err := pb.eval.Rescale(left, left); err != nil {
				return false, fmt.Errorf("Fast power %d: rescale left: %w", n, err)
			}
		}
		if rescaleB {
			if err := pb.eval.Rescale(right, right); err != nil {
				return false, fmt.Errorf("Fast power %d: rescale right: %w", n, err)
			}
		}
		degree = 2
	} else {
		if rescaleA {
			if err := pb.eval.Rescale(left, left); err != nil {
				return false, fmt.Errorf("Fast power %d: rescale left: %w", n, err)
			}
		}
		if rescaleB {
			if err := pb.eval.Rescale(right, right); err != nil {
				return false, fmt.Errorf("Fast power %d: rescale right: %w", n, err)
			}
		}
		if left.Degree() > 1 || right.Degree() > 1 {
			return false, fmt.Errorf("Fast power %d: non-relinearized operand reached strict multiplication", n)
		}
	}

	if !lazy {
		degree = 1
	}
	level := utils.Min(left.Level(), right.Level())
	out := pb.workspace.powerBuffer(pb.params, n, degree, level, pb.values[1])
	if lazy {
		err = pb.eval.Mul(left, right, out)
	} else {
		err = pb.eval.MulRelin(left, right, out)
	}
	if err != nil {
		return false, fmt.Errorf("Fast power %d: multiply: %w", n, err)
	}

	if pb.basis == bignum.Chebyshev {
		c := a - b
		if c < 0 {
			c = -c
		}
		if err = pb.eval.Add(out, out, out); err != nil {
			return false, fmt.Errorf("Fast power %d: double: %w", n, err)
		}
		if c == 0 {
			if err = pb.eval.Add(out, -1, out); err != nil {
				return false, fmt.Errorf("Fast power %d: subtract one: %w", n, err)
			}
		} else {
			if err = pb.genPower(c, lazy); err != nil {
				return false, fmt.Errorf("Fast power %d: difference power %d: %w", n, c, err)
			}
			if err = pb.workspace.subAligned(pb.params, pb.eval, out, pb.values[c]); err != nil {
				return false, fmt.Errorf("Fast power %d: subtract difference power: %w", n, err)
			}
		}
	}

	pb.values[n] = out
	return true, nil
}

func (ws *fastPolynomialWorkspace) evaluatePlan(params ckks.Parameters, eval *fastckks.Evaluator, plan commonpolynomial.PatersonStockmeyerPolynomial, powers map[int]*rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	split := len(plan.Value)
	if len(ws.babySteps) < split {
		ws.babySteps = append(ws.babySteps, make([]*fastBabyStep, split-len(ws.babySteps))...)
	}
	ws.babySteps = ws.babySteps[:split]

	for i := range ws.babySteps {
		step, err := ws.evaluateBabyStep(params, eval, plan.Value[i], powers, i)
		if err != nil {
			return nil, fmt.Errorf("Fast polynomial baby step %d: %w", i, err)
		}
		ws.babySteps[split-i-1] = step
	}
	for len(ws.babySteps) != 1 {
		if cap(ws.giantSteps) < len(ws.babySteps) {
			ws.giantSteps = make([]int, len(ws.babySteps))
		} else {
			ws.giantSteps = ws.giantSteps[:len(ws.babySteps)]
		}
		for i := range ws.giantSteps {
			ws.giantSteps[i] = 0
		}
		for i := 0; i < len(ws.babySteps); i++ {
			if i == len(ws.babySteps)-1 {
				ws.giantSteps[i] = 2
			} else if ws.babySteps[i].Degree == ws.babySteps[i+1].Degree {
				ws.giantSteps[i] = 1
				i++
			}
		}
		for i := 0; i < len(ws.babySteps); i++ {
			if ws.giantSteps[i] == 2 {
				ws.babySteps[i].Degree = ws.babySteps[i-1].Degree
			} else if ws.giantSteps[i] == 1 {
				even, odd := ws.babySteps[i], ws.babySteps[i+1]
				deg := 1 << bits.Len64(uint64(even.Degree))
				if err := ws.evaluateMonomial(params, eval, even.Value, odd.Value, powers[deg]); err != nil {
					return nil, fmt.Errorf("Fast polynomial giant step %d: %w", i, err)
				}
				odd.Degree = 2*deg - 1
				ws.babySteps[i] = nil
				i++
			}
		}
		idx := 0
		for _, step := range ws.babySteps {
			if step != nil {
				ws.babySteps[idx] = step
				idx++
			}
		}
		ws.babySteps = ws.babySteps[:idx]
	}

	result := ws.babySteps[0].Value
	if result.Degree() == 2 {
		if err := eval.Relinearize(result, result); err != nil {
			return nil, fmt.Errorf("Fast polynomial final relinearization: %w", err)
		}
	}
	if err := eval.Rescale(result, result); err != nil {
		return nil, fmt.Errorf("Fast polynomial final rescale: %w", err)
	}
	return result, nil
}

func (ws *fastPolynomialWorkspace) evaluateBabyStep(params ckks.Parameters, eval *fastckks.Evaluator, poly commonpolynomial.Polynomial, powers map[int]*rlwe.Ciphertext, index int) (*fastBabyStep, error) {
	if poly.Degree() < 0 {
		return nil, errors.New("invalid empty baby-step polynomial")
	}
	level := poly.Level
	scale := poly.Scale
	out := ws.babyBuffer(params, index, 1, level)
	zeroMaintained(out)
	*out.MetaData = *powers[1].MetaData
	out.Scale = scale

	if poly.IsEven {
		if poly.Coeffs[0] == nil {
			return nil, errors.New("nil even polynomial constant coefficient")
		}
		if err := eval.Add(out, poly.Coeffs[0], out); err != nil {
			return nil, fmt.Errorf("constant: %w", err)
		}
	}
	for key := poly.Degree(); key > 0; key-- {
		if !(poly.IsEven || poly.IsOdd) || (key&1 == 0 && poly.IsEven) || (key&1 == 1 && poly.IsOdd) {
			if poly.Coeffs[key] == nil {
				continue
			}
			if powers[key] == nil {
				return nil, fmt.Errorf("missing Fast power %d", key)
			}
			if err := eval.MulThenAdd(powers[key], poly.Coeffs[key], out); err != nil {
				return nil, fmt.Errorf("coefficient %d: %w", key, err)
			}
		}
	}
	return &fastBabyStep{Degree: poly.Degree(), Value: out}, nil
}

func (ws *fastPolynomialWorkspace) evaluateMonomial(params ckks.Parameters, eval *fastckks.Evaluator, a, b, xpow *rlwe.Ciphertext) error {
	if xpow == nil {
		return errors.New("missing Fast giant-step power")
	}
	if b.Degree() == 2 {
		if err := eval.Relinearize(b, b); err != nil {
			return fmt.Errorf("relinearize: %w", err)
		}
	}
	if err := eval.Rescale(b, b); err != nil {
		return fmt.Errorf("rescale: %w", err)
	}
	if err := eval.Mul(b, xpow, b); err != nil {
		return fmt.Errorf("multiply: %w", err)
	}
	if !a.Scale.InDelta(b.Scale, float64(rlwe.ScalePrecision-12)) {
		return fmt.Errorf("scale discrepancy: (rescale(b) * X^n).Scale = %v != a.Scale = %v", &b.Scale.Value, &a.Scale.Value)
	}
	if err := ws.addAligned(params, eval, a, b); err != nil {
		return fmt.Errorf("add: %w", err)
	}
	return nil
}

func copyMaintained(src, dst *rlwe.Ciphertext) {
	dst.Resize(1, src.Level())
	*dst.MetaData = *src.MetaData
	dst.IsNTT = src.IsNTT
	dst.IsMontgomery = src.IsMontgomery
	dst.Scale = src.Scale
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			copy(dst.Value[d].Coeffs[limb], src.Value[d].Coeffs[limb])
		}
	}
}

func (ws *fastPolynomialWorkspace) subAligned(params ckks.Parameters, eval *fastckks.Evaluator, out, sub *rlwe.Ciphertext) error {
	if out.Scale.Equal(sub.Scale) {
		return eval.Sub(out, sub, out)
	}

	if out.Scale.Cmp(sub.Scale) > 0 {
		ws.scaleScratch = ws.ensureCiphertext(params, ws.scaleScratch, sub.Degree(), utils.Min(out.Level(), sub.Level()))
		ws.scaleScratch.IsNTT = sub.IsNTT
		ws.scaleScratch.IsMontgomery = sub.IsMontgomery
		*ws.scaleScratch.MetaData = *sub.MetaData
		ratio := out.Scale.Div(sub.Scale).BigInt()
		if ratio.Sign() <= 0 {
			return errors.New("invalid Fast polynomial scale alignment ratio")
		}
		if err := eval.MulIntegerMaintained(sub, ratio, ws.scaleScratch); err != nil {
			return err
		}
		ws.scaleScratch.Scale = out.Scale
		return eval.Sub(out, ws.scaleScratch, out)
	}

	ws.scaleScratch = ws.ensureCiphertext(params, ws.scaleScratch, out.Degree(), utils.Min(out.Level(), sub.Level()))
	ws.scaleScratch.IsNTT = out.IsNTT
	ws.scaleScratch.IsMontgomery = out.IsMontgomery
	*ws.scaleScratch.MetaData = *out.MetaData
	ratio := sub.Scale.Div(out.Scale).BigInt()
	if ratio.Sign() <= 0 {
		return errors.New("invalid Fast polynomial scale alignment ratio")
	}
	if err := eval.MulIntegerMaintained(out, ratio, ws.scaleScratch); err != nil {
		return err
	}
	ws.scaleScratch.Scale = sub.Scale
	if err := eval.Sub(ws.scaleScratch, sub, ws.scaleScratch); err != nil {
		return err
	}
	copyMaintainedElement(ws.scaleScratch, out)
	return nil
}

func (ws *fastPolynomialWorkspace) addAligned(params ckks.Parameters, eval *fastckks.Evaluator, a, b *rlwe.Ciphertext) error {
	if a.Scale.Equal(b.Scale) {
		return eval.Add(b, a, b)
	}

	if b.Scale.Cmp(a.Scale) > 0 {
		ws.scaleScratch = ws.ensureCiphertext(params, ws.scaleScratch, a.Degree(), utils.Min(a.Level(), b.Level()))
		ws.scaleScratch.IsNTT = a.IsNTT
		ws.scaleScratch.IsMontgomery = a.IsMontgomery
		*ws.scaleScratch.MetaData = *a.MetaData
		ratio := b.Scale.Div(a.Scale).BigInt()
		if ratio.Sign() <= 0 {
			return errors.New("invalid Fast polynomial scale alignment ratio")
		}
		if err := eval.MulIntegerMaintained(a, ratio, ws.scaleScratch); err != nil {
			return err
		}
		ws.scaleScratch.Scale = b.Scale
		return eval.Add(b, ws.scaleScratch, b)
	}

	ws.scaleScratch = ws.ensureCiphertext(params, ws.scaleScratch, b.Degree(), utils.Min(a.Level(), b.Level()))
	ws.scaleScratch.IsNTT = b.IsNTT
	ws.scaleScratch.IsMontgomery = b.IsMontgomery
	*ws.scaleScratch.MetaData = *b.MetaData
	ratio := a.Scale.Div(b.Scale).BigInt()
	if ratio.Sign() <= 0 {
		return errors.New("invalid Fast polynomial scale alignment ratio")
	}
	if err := eval.MulIntegerMaintained(b, ratio, ws.scaleScratch); err != nil {
		return err
	}
	ws.scaleScratch.Scale = a.Scale
	if err := eval.Add(ws.scaleScratch, a, ws.scaleScratch); err != nil {
		return err
	}
	copyMaintainedElement(ws.scaleScratch, b)
	return nil
}

func copyMaintainedElement(src, dst *rlwe.Ciphertext) {
	dst.Resize(src.Degree(), src.Level())
	*dst.MetaData = *src.MetaData
	dst.IsNTT = src.IsNTT
	dst.IsMontgomery = src.IsMontgomery
	dst.Scale = src.Scale
	for d := range src.Value {
		for limb := 0; limb < 2; limb++ {
			copy(dst.Value[d].Coeffs[limb], src.Value[d].Coeffs[limb])
		}
	}
}

// cloneMaintainedResult creates an independently-owned public ciphertext
// without reading dormant q2...qL rows from the Fast workspace.
func cloneMaintainedResult(params ckks.Parameters, src *rlwe.Ciphertext) *rlwe.Ciphertext {
	dst := ckks.NewCiphertext(params, src.Degree(), src.Level())
	copyMaintainedElement(src, dst)
	return dst
}

func zeroMaintained(ct *rlwe.Ciphertext) {
	for d := range ct.Value {
		for limb := 0; limb < 2 && limb < len(ct.Value[d].Coeffs); limb++ {
			ring.ZeroVec(ct.Value[d].Coeffs[limb])
		}
	}
}
