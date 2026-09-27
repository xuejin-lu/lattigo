package polynomial

import (
	"errors"
	"fmt"
	"math/big"
	"math/bits"
	"sort"

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

// GuardSelectionEvidence reports the last explicit guarded-plan selection.
// It is intentionally observational; it does not alter polynomial evaluation.
type GuardSelectionEvidence struct {
	PlanScaleOverride bool
	PlanIndex         int
	PlanCount         int
	ScalarDegree      int
	Operations        int
	Contraction       bool
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

// LastGuardSelectionEvidence returns compact evidence for the last guarded
// evaluation. A zero value means no guarded evaluation has completed.
func (eval *FastEvaluator) LastGuardSelectionEvidence() GuardSelectionEvidence {
	if eval == nil {
		return GuardSelectionEvidence{}
	}
	ws := &eval.workspace
	return GuardSelectionEvidence{
		PlanScaleOverride: ws.planScaleOverride,
		PlanIndex:         ws.guardedPlanIndex,
		PlanCount:         ws.guardedPlanCount,
		ScalarDegree:      ws.guardedScalarDegree,
		Operations:        ws.guardedOperations,
		Contraction:       ws.guardedContraction,
	}
}

// Evaluate evaluates the single polynomial p on input. The current Fast
// polynomial surface accepts Chebyshev polynomials in the Standard ring with
// NTT/Montgomery degree-one input at a level containing q0 and q1.
func (eval *FastEvaluator) Evaluate(input *rlwe.Ciphertext, p bignum.Polynomial, targetScale rlwe.Scale) (*rlwe.Ciphertext, error) {
	return eval.evaluate(input, p, targetScale, nil, false)
}

// EvaluateWithPlanScale evaluates the same Paterson-Stockmeyer plan as
// Evaluate while overriding only the planned polynomial-value scales. The
// planner still receives targetScale and evaluatePlan retains its existing
// final relinearization and single Rescale boundary.
func (eval *FastEvaluator) EvaluateWithPlanScale(input *rlwe.Ciphertext, p bignum.Polynomial, targetScale, planScale rlwe.Scale) (*rlwe.Ciphertext, error) {
	if planScale.Value.Sign() <= 0 {
		return nil, errors.New("Fast polynomial plan scale must be positive")
	}
	return eval.evaluate(input, p, targetScale, &planScale, false)
}

// EvaluateWithPlanScaleFinalParentOneBitScalarGuard evaluates the normalized
// LogN13 plan with one explicit final-parent scalar guard. The ordinary
// Evaluate and EvaluateWithPlanScale paths remain unguarded.
func (eval *FastEvaluator) EvaluateWithPlanScaleFinalParentOneBitScalarGuard(input *rlwe.Ciphertext, p bignum.Polynomial, targetScale, planScale rlwe.Scale) (*rlwe.Ciphertext, error) {
	if planScale.Value.Sign() <= 0 {
		return nil, errors.New("Fast polynomial guarded plan scale must be positive")
	}
	return eval.evaluate(input, p, targetScale, &planScale, true)
}

func (eval *FastEvaluator) evaluate(input *rlwe.Ciphertext, p bignum.Polynomial, targetScale rlwe.Scale, planScale *rlwe.Scale, guardFinalParent bool) (*rlwe.Ciphertext, error) {
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

	rows, err := fastckks.QPrefixWidth(input.Level())
	if err != nil {
		return nil, fmt.Errorf("Fast polynomial input Q-prefix: %w", err)
	}
	ws := &eval.workspace
	ws.reset(eval.Parameters, input, rows)
	if err := eval.Evaluator.ObserveQPrefixCapacity("polynomial-entry", input, rows); err != nil {
		return nil, err
	}

	// Degree zero has no power-generation or PS planning work. This also keeps
	// the bounded Fast surface well-defined for constant Chebyshev polynomials.
	if commonPoly.Degree() == 0 {
		if p.Coeffs[0] == nil {
			return nil, errors.New("Fast polynomial constant coefficient cannot be nil")
		}
		out := ws.babyBuffer(eval.Parameters, 0, 1, input.Level())
		zeroQPrefix(out, rows)
		*out.MetaData = *input.MetaData
		out.Scale = targetScale
		if err := eval.Evaluator.AddScalarQPrefixRows(out, p.Coeffs[0], rows, out); err != nil {
			return nil, fmt.Errorf("Fast polynomial constant: %w", err)
		}
		return cloneQPrefixResult(eval.Parameters, out, rows), nil
	}

	if err := ws.generatePowers(eval.Parameters, eval.Evaluator, p, commonPoly); err != nil {
		return nil, err
	}
	powerKeys := make([]int, 0, len(ws.powers))
	for power := range ws.powers {
		powerKeys = append(powerKeys, power)
	}
	sort.Ints(powerKeys)
	for _, power := range powerKeys {
		powerRows, err := ws.rowsAt(ws.powers[power].Level())
		if err != nil {
			return nil, err
		}
		if err := eval.Evaluator.ObserveQPrefixCapacity(fmt.Sprintf("generated-power-%d", power), ws.powers[power], powerRows); err != nil {
			return nil, err
		}
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
	if planScale != nil {
		ws.planScaleOverride = true
		for i := range plan.Value {
			plan.Value[i].Scale = *planScale
		}
	}
	result, err := ws.evaluatePlan(eval.Parameters, eval.Evaluator, plan, ws.powers, guardFinalParent)
	if err != nil {
		return nil, err
	}
	resultRows, err := ws.rowsAt(result.Level())
	if err != nil {
		return nil, err
	}
	return cloneQPrefixResult(eval.Parameters, result, resultRows), nil
}

type fastPolynomialWorkspace struct {
	x1                  *rlwe.Ciphertext
	rows                int
	powers              map[int]*rlwe.Ciphertext
	powerBuffers        map[int]*rlwe.Ciphertext
	balancedLeft        *rlwe.Ciphertext
	balancedRight       *rlwe.Ciphertext
	babyBuffers         []*rlwe.Ciphertext
	babySteps           []*fastBabyStep
	giantSteps          []int
	scaleScratch        *rlwe.Ciphertext
	planScaleOverride   bool
	guardedPlanIndex    int
	guardedPlanCount    int
	guardedScalarDegree int
	guardedOperations   int
	guardedContraction  bool
}

type fastBabyStep struct {
	Degree              int
	Value               *rlwe.Ciphertext
	GuardApplied        bool
	GuardedScalarDegree int
}

func (ws *fastPolynomialWorkspace) reset(params ckks.Parameters, input *rlwe.Ciphertext, rows int) {
	ws.planScaleOverride = false
	ws.guardedPlanIndex = -1
	ws.guardedPlanCount = 0
	ws.guardedScalarDegree = 0
	ws.guardedOperations = 0
	ws.guardedContraction = false
	ws.rows = rows
	ws.x1 = ws.ensureCiphertext(params, ws.x1, 1, input.Level())
	copyQPrefix(params, input, ws.x1, rows)
	for key := range ws.powers {
		delete(ws.powers, key)
	}
	ws.powers[1] = ws.x1
}

func (ws *fastPolynomialWorkspace) rowsAt(level int) (int, error) {
	width, err := fastckks.QPrefixWidth(level)
	if err != nil {
		return 0, err
	}
	if ws.rows == 0 || ws.rows > width {
		return width, nil
	}
	return ws.rows, nil
}

func (ws *fastPolynomialWorkspace) ensureCiphertext(params ckks.Parameters, ct *rlwe.Ciphertext, degree, level int) *rlwe.Ciphertext {
	if ct == nil {
		ct = fastckks.NewCiphertext(params, degree, level)
	} else {
		fastckks.Resize(ct, degree, level, params.N())
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

const (
	balancedMinScaleBits     = 20
	balancedScaleMarginBits  = 112
	balancedFactorCandidates = 2
)

type balancedFactors struct {
	left  uint64
	right uint64
}

// balancedScaleToleranceBits intentionally leaves a fixed precision margin
// for the integer near-square factorization. It is checked with Scale.InDelta,
// not a float64 comparison, and retains at least 16 relative bits under the
// 128-bit CKKS scale representation.
const balancedScaleToleranceBits = float64(rlwe.ScalePrecision - balancedScaleMarginBits)

func integerSqrt(value uint64) uint64 {
	if value < 2 {
		return value
	}
	root := uint64(1) << ((bits.Len64(value) + 1) / 2)
	for {
		next := (root + value/root) >> 1
		if next >= root {
			break
		}
		root = next
	}
	for (root + 1) <= value/(root+1) {
		root++
	}
	for root > value/root {
		root--
	}
	return root
}

func balancedFactorPair(q uint64) (balancedFactors, error) {
	if q < 2 {
		return balancedFactors{}, errors.New("balanced factor divisor must be at least two")
	}
	root := integerSqrt(q)
	var best balancedFactors
	var bestAbs *big.Int
	var bestMax uint64
	for delta := -balancedFactorCandidates; delta <= balancedFactorCandidates; delta++ {
		candidate := int64(root) + int64(delta)
		if candidate <= 0 {
			continue
		}
		m1 := uint64(candidate)
		quotient, remainder := q/m1, q%m1
		candidates := []uint64{quotient}
		if remainder != 0 {
			candidates = append(candidates, quotient+1)
		}
		for _, m2 := range candidates {
			for _, pair := range [][2]uint64{{m1, m2}, {m2, m1}} {
				product := new(big.Int).Mul(new(big.Int).SetUint64(pair[0]), new(big.Int).SetUint64(pair[1]))
				difference := new(big.Int).Sub(new(big.Int).Set(product), new(big.Int).SetUint64(q))
				absDifference := new(big.Int).Abs(difference)
				pairMax := pair[0]
				if pair[1] > pairMax {
					pairMax = pair[1]
				}
				better := bestAbs == nil || absDifference.Cmp(bestAbs) < 0
				if !better && bestAbs != nil && absDifference.Cmp(bestAbs) == 0 {
					better = pairMax < bestMax || (pairMax == bestMax && (pair[0] < best.left || (pair[0] == best.left && pair[1] < best.right)))
				}
				if better {
					best = balancedFactors{left: pair[0], right: pair[1]}
					bestAbs = absDifference
					bestMax = pairMax
				}
			}
		}
	}
	if bestAbs == nil || best.left == 0 || best.right == 0 {
		return balancedFactors{}, fmt.Errorf("cannot select balanced factors for divisor %d", q)
	}
	return best, nil
}

type balancedSchedule struct {
	factors     balancedFactors
	leftScale   rlwe.Scale
	rightScale  rlwe.Scale
	targetScale rlwe.Scale
	balanced    bool
}

func (pb *fastPowerBasis) balancedScheduleFor(left, right *rlwe.Ciphertext, commonLevel int) (balancedSchedule, error) {
	levelsConsumed := pb.params.LevelsConsumedPerRescaling()
	if commonLevel < levelsConsumed {
		return balancedSchedule{}, fmt.Errorf("common level %d is too low for %d rescale levels", commonLevel, levelsConsumed)
	}
	if levelsConsumed != 1 {
		return balancedSchedule{}, nil
	}
	factors, err := balancedFactorPair(pb.params.Q()[commonLevel])
	if err != nil {
		return balancedSchedule{}, err
	}
	qScale := rlwe.NewScale(pb.params.Q()[commonLevel])
	leftScale := left.Scale.Mul(rlwe.NewScale(factors.left)).Div(qScale)
	rightScale := right.Scale.Mul(rlwe.NewScale(factors.right)).Div(qScale)
	targetScale := left.Scale.Mul(right.Scale).Div(qScale)
	minimumScale := rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), balancedMinScaleBits))
	return balancedSchedule{factors: factors, leftScale: leftScale, rightScale: rightScale, targetScale: targetScale, balanced: leftScale.Cmp(minimumScale) >= 0 && rightScale.Cmp(minimumScale) >= 0}, nil
}

// postProductLogN13Schedule selects the proven LogN13 Chebyshev power
// schedule. In this bounded profile the generated power must keep the full
// product scale through the Chebyshev recurrence and consume exactly one
// level only after the recurrence is complete. Other profiles retain the
// established balanced/low-scale schedule.
func (pb *fastPowerBasis) postProductLogN13Schedule(commonLevel int) bool {
	q := pb.params.Q()
	isP93Schedule := len(q) >= 3 && bits.Len64(q[0]) == 56 && bits.Len64(q[1]) <= 40 && bits.Len64(q[2]) <= 40
	return pb.basis == bignum.Chebyshev &&
		pb.params.LogN() == 13 &&
		pb.params.RingType() == ring.Standard &&
		pb.params.LevelsConsumedPerRescaling() == 1 &&
		commonLevel >= 2 &&
		isP93Schedule
}

func (pb *fastPowerBasis) genPower(n int, lazy bool) error {
	if pb.values[n] != nil {
		return nil
	}
	return pb.genPowerInternal(n, lazy)
}

func (pb *fastPowerBasis) genPowerInternal(n int, lazy bool) error {
	if pb.values[n] != nil {
		return nil
	}
	a, b := commonpolynomial.SplitDegree(n)
	isPow2 := n&(n-1) == 0

	if err := pb.genPowerInternal(a, lazy && !isPow2); err != nil {
		return fmt.Errorf("Fast power %d: power %d: %w", n, a, err)
	}
	if err := pb.genPowerInternal(b, lazy && !isPow2); err != nil {
		return fmt.Errorf("Fast power %d: power %d: %w", n, b, err)
	}

	left, right := pb.values[a], pb.values[b]
	degree := 1
	var err error
	if lazy {
		degree = 2
	} else {
		if left.Degree() > 1 || right.Degree() > 1 {
			return fmt.Errorf("Fast power %d: non-relinearized operand reached strict multiplication", n)
		}
	}

	if !lazy {
		degree = 1
	}
	commonLevel := utils.Min(left.Level(), right.Level())
	rows, err := pb.workspace.rowsAt(commonLevel)
	if err != nil {
		return fmt.Errorf("Fast power %d Q-prefix: %w", n, err)
	}
	schedule, err := pb.balancedScheduleFor(left, right, commonLevel)
	if err != nil {
		return fmt.Errorf("Fast power %d: balanced schedule: %w", n, err)
	}
	var out *rlwe.Ciphertext
	postProduct := pb.postProductLogN13Schedule(commonLevel)
	balanced := schedule.balanced && !postProduct
	if balanced {
		leftCopy, err := pb.workspace.balancedCopy(pb.params, left, commonLevel, rows, true)
		if err != nil {
			return fmt.Errorf("Fast power %d: copy balanced left: %w", n, err)
		}
		rightCopy, err := pb.workspace.balancedCopy(pb.params, right, commonLevel, rows, false)
		if err != nil {
			return fmt.Errorf("Fast power %d: copy balanced right: %w", n, err)
		}
		if lazy {
			if leftCopy.Degree() == 2 {
				if err := pb.eval.RelinearizeQPrefixRows(leftCopy, leftCopy, rows); err != nil {
					return fmt.Errorf("Fast power %d: relinearize balanced left: %w", n, err)
				}
			}
			if rightCopy.Degree() == 2 {
				if err := pb.eval.RelinearizeQPrefixRows(rightCopy, rightCopy, rows); err != nil {
					return fmt.Errorf("Fast power %d: relinearize balanced right: %w", n, err)
				}
			}
		}
		if err := pb.eval.MulIntegerQPrefixRows(leftCopy, new(big.Int).SetUint64(schedule.factors.left), rows, leftCopy); err != nil {
			return fmt.Errorf("Fast power %d: scale balanced left: %w", n, err)
		}
		if err := pb.eval.MulIntegerQPrefixRows(rightCopy, new(big.Int).SetUint64(schedule.factors.right), rows, rightCopy); err != nil {
			return fmt.Errorf("Fast power %d: scale balanced right: %w", n, err)
		}
		leftCopy.Scale = left.Scale.Mul(rlwe.NewScale(schedule.factors.left))
		rightCopy.Scale = right.Scale.Mul(rlwe.NewScale(schedule.factors.right))
		if err := pb.eval.RescaleQPrefixRows(leftCopy, rows, leftCopy); err != nil {
			return fmt.Errorf("Fast power %d: balanced left rescale: %w", n, err)
		}
		if err := pb.eval.RescaleQPrefixRows(rightCopy, rows, rightCopy); err != nil {
			return fmt.Errorf("Fast power %d: balanced right rescale: %w", n, err)
		}
		if leftCopy.Level() != commonLevel-1 || rightCopy.Level() != commonLevel-1 {
			return fmt.Errorf("Fast power %d: balanced operands consumed unexpected levels", n)
		}
		out = pb.workspace.powerBuffer(pb.params, n, degree, commonLevel-1, pb.values[1])
		productRows, err := pb.workspace.rowsAt(commonLevel - 1)
		if err != nil {
			return fmt.Errorf("Fast power %d post-rescale Q-prefix: %w", n, err)
		}
		if lazy {
			err = pb.eval.MulElementQPrefixRows(leftCopy, rightCopy.El(), productRows, out)
		} else {
			err = pb.eval.MulRelinElementQPrefixRows(leftCopy, rightCopy.El(), productRows, out)
		}
	} else {
		// Low-scale inputs retain the established post-product schedule because
		// balanced pre-Rescale would otherwise take an operand below the
		// documented precision floor.
		if lazy {
			if left.Degree() == 2 {
				if err := pb.eval.RelinearizeQPrefixRows(left, left, rows); err != nil {
					return fmt.Errorf("Fast power %d: relinearize left: %w", n, err)
				}
			}
			if right.Degree() == 2 {
				if err := pb.eval.RelinearizeQPrefixRows(right, right, rows); err != nil {
					return fmt.Errorf("Fast power %d: relinearize right: %w", n, err)
				}
			}
		}
		out = pb.workspace.powerBuffer(pb.params, n, degree, commonLevel, pb.values[1])
		if lazy {
			err = pb.eval.MulElementQPrefixRows(left, right.El(), rows, out)
		} else {
			err = pb.eval.MulRelinElementQPrefixRows(left, right.El(), rows, out)
		}
	}
	if err != nil {
		return fmt.Errorf("Fast power %d: multiply: %w", n, err)
	}

	if pb.basis == bignum.Chebyshev {
		outRows, rowErr := pb.workspace.rowsAt(out.Level())
		if rowErr != nil {
			return fmt.Errorf("Fast power %d doubling Q-prefix: %w", n, rowErr)
		}
		if err = pb.eval.AddQPrefixRows(out, out, out, outRows); err != nil {
			return fmt.Errorf("Fast power %d: double: %w", n, err)
		}
	}
	if balanced {
		if !out.Scale.InDelta(schedule.targetScale, balancedScaleToleranceBits) {
			return fmt.Errorf("Fast power %d: balanced scale discrepancy: actual=%v target=%v", n, &out.Scale.Value, &schedule.targetScale.Value)
		}
		out.Scale = schedule.targetScale
	} else if !postProduct {
		outRows, rowErr := pb.workspace.rowsAt(out.Level())
		if rowErr != nil {
			return fmt.Errorf("Fast power %d rescale Q-prefix: %w", n, rowErr)
		}
		if err := pb.eval.RescaleQPrefixRows(out, outRows, out); err != nil {
			return fmt.Errorf("Fast power %d: rescale: %w", n, err)
		}
	}

	if pb.basis == bignum.Chebyshev {
		c := a - b
		if c < 0 {
			c = -c
		}
		if c == 0 {
			outRows, rowErr := pb.workspace.rowsAt(out.Level())
			if rowErr != nil {
				return fmt.Errorf("Fast power %d subtract-one Q-prefix: %w", n, rowErr)
			}
			if err = pb.eval.AddScalarQPrefixRows(out, -1, outRows, out); err != nil {
				return fmt.Errorf("Fast power %d: subtract one: %w", n, err)
			}
		} else {
			if err = pb.genPower(c, lazy); err != nil {
				return fmt.Errorf("Fast power %d: difference power %d: %w", n, c, err)
			}
			if err = pb.workspace.subAligned(pb.params, pb.eval, out, pb.values[c]); err != nil {
				return fmt.Errorf("Fast power %d: subtract difference power: %w", n, err)
			}
		}
	}
	if postProduct {
		outRows, rowErr := pb.workspace.rowsAt(out.Level())
		if rowErr != nil {
			return fmt.Errorf("Fast power %d post-product Q-prefix: %w", n, rowErr)
		}
		if err := pb.eval.RescaleQPrefixRows(out, outRows, out); err != nil {
			return fmt.Errorf("Fast power %d: post-product rescale: %w", n, err)
		}
	}

	pb.values[n] = out
	return nil
}

func (ws *fastPolynomialWorkspace) evaluatePlan(params ckks.Parameters, eval *fastckks.Evaluator, plan commonpolynomial.PatersonStockmeyerPolynomial, powers map[int]*rlwe.Ciphertext, guardFinalParent bool) (*rlwe.Ciphertext, error) {
	split := len(plan.Value)
	if guardFinalParent && !ws.planScaleOverride {
		return nil, errors.New("Fast polynomial guarded path requires an active plan-scale override")
	}
	guardedOperations := 0
	ws.guardedPlanCount = split
	if len(ws.babySteps) < split {
		ws.babySteps = append(ws.babySteps, make([]*fastBabyStep, split-len(ws.babySteps))...)
	}
	ws.babySteps = ws.babySteps[:split]

	for i := range ws.babySteps {
		step, err := ws.evaluateBabyStep(params, eval, plan.Value[i], powers, i, guardFinalParent && i == split-1)
		if err != nil {
			return nil, fmt.Errorf("Fast polynomial baby step %d: %w", i, err)
		}
		if step.GuardApplied {
			guardedOperations++
			ws.guardedPlanIndex = i
			ws.guardedScalarDegree = step.GuardedScalarDegree
			ws.guardedContraction = true
		}
		ws.babySteps[split-i-1] = step
		rows, err := ws.rowsAt(step.Value.Level())
		if err != nil {
			return nil, err
		}
		if err := eval.ObserveQPrefixCapacity(fmt.Sprintf("ps-baby-%d", i), step.Value, rows); err != nil {
			return nil, err
		}
	}
	if guardFinalParent && guardedOperations != 1 {
		return nil, fmt.Errorf("Fast polynomial guarded path selected %d final-parent scalar operations, want exactly one", guardedOperations)
	}
	ws.guardedOperations = guardedOperations
	giantCheckpointIndex := 0
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
				rows, err := ws.rowsAt(odd.Value.Level())
				if err != nil {
					return nil, err
				}
				if err := eval.ObserveQPrefixCapacity(fmt.Sprintf("ps-giant-%d", giantCheckpointIndex), odd.Value, rows); err != nil {
					return nil, err
				}
				giantCheckpointIndex++
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
	rows, err := ws.rowsAt(result.Level())
	if err != nil {
		return nil, fmt.Errorf("Fast polynomial final Q-prefix: %w", err)
	}
	if result.Degree() == 2 {
		if err := eval.RelinearizeQPrefixRows(result, result, rows); err != nil {
			return nil, fmt.Errorf("Fast polynomial final relinearization: %w", err)
		}
	}
	if err := eval.ObserveQPrefixCapacity("polynomial-before-final-rescale", result, rows); err != nil {
		return nil, err
	}
	if err := eval.RescaleQPrefixRows(result, rows, result); err != nil {
		return nil, fmt.Errorf("Fast polynomial final rescale: %w", err)
	}
	resultRows, err := ws.rowsAt(result.Level())
	if err != nil {
		return nil, err
	}
	if err := eval.ObserveQPrefixCapacity("polynomial-after-final-rescale", result, resultRows); err != nil {
		return nil, err
	}
	return result, nil
}

func (ws *fastPolynomialWorkspace) evaluateBabyStep(params ckks.Parameters, eval *fastckks.Evaluator, poly commonpolynomial.Polynomial, powers map[int]*rlwe.Ciphertext, index int, guardFinalParent bool) (*fastBabyStep, error) {
	if poly.Degree() < 0 {
		return nil, errors.New("invalid empty baby-step polynomial")
	}
	level := poly.Level
	scale := poly.Scale
	rows, err := ws.rowsAt(level)
	if err != nil {
		return nil, fmt.Errorf("baby-step Q-prefix: %w", err)
	}
	out := ws.babyBuffer(params, index, 1, level)
	zeroQPrefix(out, rows)
	*out.MetaData = *powers[1].MetaData
	out.Scale = scale

	if poly.IsEven {
		if poly.Coeffs[0] == nil {
			return nil, errors.New("nil even polynomial constant coefficient")
		}
		if err := eval.AddScalarQPrefixRows(out, poly.Coeffs[0], rows, out); err != nil {
			return nil, fmt.Errorf("constant: %w", err)
		}
	}
	guardKey := -1
	if guardFinalParent {
		candidateKeys := make([]int, 0)
		for key := poly.Degree(); key > 0; key-- {
			if !(poly.IsEven || poly.IsOdd) || (key&1 == 0 && poly.IsEven) || (key&1 == 1 && poly.IsOdd) {
				coefficient := poly.Coeffs[key]
				if coefficient != nil && !coefficient.IsInt() && !fastPolynomialCoefficientZero(coefficient) {
					candidateKeys = append(candidateKeys, key)
					guardKey = key
				}
			}
		}
		if guardKey != 2 {
			return nil, fmt.Errorf("Fast polynomial guarded path selected final scalar degree %d, want 2 (candidates=%v even=%t odd=%t)", guardKey, candidateKeys, poly.IsEven, poly.IsOdd)
		}
	}
	guardApplied := false
	for key := poly.Degree(); key > 0; key-- {
		if !(poly.IsEven || poly.IsOdd) || (key&1 == 0 && poly.IsEven) || (key&1 == 1 && poly.IsOdd) {
			if poly.Coeffs[key] == nil {
				continue
			}
			if powers[key] == nil {
				return nil, fmt.Errorf("missing Fast power %d", key)
			}
			var err error
			if guardFinalParent && key == guardKey {
				err = eval.MulThenAddOneBitScalarGuardQPrefixRows(powers[key], poly.Coeffs[key], rows, out)
				guardApplied = true
			} else {
				err = eval.MulThenAddQPrefixRows(powers[key], poly.Coeffs[key], rows, out)
			}
			if err != nil {
				return nil, fmt.Errorf("coefficient %d: %w", key, err)
			}
		}
	}
	if guardFinalParent && !guardApplied {
		return nil, errors.New("Fast polynomial guarded path did not execute its selected scalar operation")
	}
	return &fastBabyStep{Degree: poly.Degree(), Value: out, GuardApplied: guardApplied, GuardedScalarDegree: guardKey}, nil
}

func fastPolynomialCoefficientZero(coefficient *bignum.Complex) bool {
	return coefficient[0].Sign() == 0 && coefficient[1].Sign() == 0
}

func (ws *fastPolynomialWorkspace) evaluateMonomial(params ckks.Parameters, eval *fastckks.Evaluator, a, b, xpow *rlwe.Ciphertext) error {
	if xpow == nil {
		return errors.New("missing Fast giant-step power")
	}
	rows, err := ws.rowsAt(min(min(a.Level(), b.Level()), xpow.Level()))
	if err != nil {
		return err
	}
	if b.Degree() == 2 {
		if err := eval.RelinearizeQPrefixRows(b, b, rows); err != nil {
			return fmt.Errorf("relinearize: %w", err)
		}
	}
	if err := eval.RescaleQPrefixRows(b, rows, b); err != nil {
		return fmt.Errorf("rescale: %w", err)
	}
	rows, err = ws.rowsAt(min(b.Level(), xpow.Level()))
	if err != nil {
		return err
	}
	if err := eval.MulElementQPrefixRows(b, xpow.El(), rows, b); err != nil {
		return fmt.Errorf("multiply: %w", err)
	}
	toleranceBits := float64(rlwe.ScalePrecision - 12)
	if ws.planScaleOverride {
		toleranceBits = 32
	}
	if !a.Scale.InDelta(b.Scale, toleranceBits) {
		return fmt.Errorf("scale discrepancy: (rescale(b) * X^n).Scale = %v != a.Scale = %v", &b.Scale.Value, &a.Scale.Value)
	}
	if ws.planScaleOverride && !a.Scale.Equal(b.Scale) {
		b.Scale = a.Scale
		rows, err = ws.rowsAt(min(a.Level(), b.Level()))
		if err != nil {
			return err
		}
		return eval.AddQPrefixRows(b, a, b, rows)
	}
	if err := ws.addAligned(params, eval, a, b); err != nil {
		return fmt.Errorf("add: %w", err)
	}
	return nil
}

func copyQPrefix(params ckks.Parameters, src, dst *rlwe.Ciphertext, rows int) {
	_ = copyQPrefixAtLevel(params, src, dst, src.Level(), rows)
}

func copyQPrefixAtLevel(params ckks.Parameters, src, dst *rlwe.Ciphertext, level, rows int) error {
	if src == nil || dst == nil {
		return errors.New("Fast polynomial Q-prefix copy operands cannot be nil")
	}
	if level < 0 || level > src.Level() {
		return fmt.Errorf("Fast polynomial maintained copy level %d is outside source level %d", level, src.Level())
	}
	fastckks.Resize(dst, src.Degree(), level, params.N())
	*dst.MetaData = *src.MetaData
	dst.IsNTT = src.IsNTT
	dst.IsMontgomery = src.IsMontgomery
	dst.Scale = src.Scale
	width, err := fastckks.QPrefixWidth(level)
	if err != nil {
		return err
	}
	if rows < 1 || rows > width {
		return fmt.Errorf("Fast polynomial Q-prefix rows %d exceed Level %d width %d", rows, level, width)
	}
	for d := range src.Value {
		for limb := 0; limb < rows; limb++ {
			copy(dst.Value[d].Coeffs[limb], src.Value[d].Coeffs[limb])
		}
	}
	return nil
}

func (ws *fastPolynomialWorkspace) balancedCopy(params ckks.Parameters, src *rlwe.Ciphertext, level, rows int, left bool) (*rlwe.Ciphertext, error) {
	if left {
		ws.balancedLeft = ws.ensureCiphertext(params, ws.balancedLeft, src.Degree(), level)
		if err := copyQPrefixAtLevel(params, src, ws.balancedLeft, level, rows); err != nil {
			return nil, err
		}
		return ws.balancedLeft, nil
	}
	ws.balancedRight = ws.ensureCiphertext(params, ws.balancedRight, src.Degree(), level)
	if err := copyQPrefixAtLevel(params, src, ws.balancedRight, level, rows); err != nil {
		return nil, err
	}
	return ws.balancedRight, nil
}

func (ws *fastPolynomialWorkspace) subAligned(params ckks.Parameters, eval *fastckks.Evaluator, out, sub *rlwe.Ciphertext) error {
	rows, err := ws.rowsAt(min(out.Level(), sub.Level()))
	if err != nil {
		return err
	}
	if out.Scale.Equal(sub.Scale) {
		return eval.SubQPrefixRows(out, sub, out, rows)
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
		if err := eval.MulIntegerQPrefixRows(sub, ratio, rows, ws.scaleScratch); err != nil {
			return err
		}
		ws.scaleScratch.Scale = out.Scale
		return eval.SubQPrefixRows(out, ws.scaleScratch, out, rows)
	}

	ws.scaleScratch = ws.ensureCiphertext(params, ws.scaleScratch, out.Degree(), utils.Min(out.Level(), sub.Level()))
	ws.scaleScratch.IsNTT = out.IsNTT
	ws.scaleScratch.IsMontgomery = out.IsMontgomery
	*ws.scaleScratch.MetaData = *out.MetaData
	ratio := sub.Scale.Div(out.Scale).BigInt()
	if ratio.Sign() <= 0 {
		return errors.New("invalid Fast polynomial scale alignment ratio")
	}
	if err := eval.MulIntegerQPrefixRows(out, ratio, rows, ws.scaleScratch); err != nil {
		return err
	}
	ws.scaleScratch.Scale = sub.Scale
	if err := eval.SubQPrefixRows(ws.scaleScratch, sub, ws.scaleScratch, rows); err != nil {
		return err
	}
	copyQPrefixElement(params, ws.scaleScratch, out, rows)
	return nil
}

func (ws *fastPolynomialWorkspace) addAligned(params ckks.Parameters, eval *fastckks.Evaluator, a, b *rlwe.Ciphertext) error {
	rows, err := ws.rowsAt(min(a.Level(), b.Level()))
	if err != nil {
		return err
	}
	if a.Scale.Equal(b.Scale) {
		return eval.AddQPrefixRows(b, a, b, rows)
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
		if err := eval.MulIntegerQPrefixRows(a, ratio, rows, ws.scaleScratch); err != nil {
			return err
		}
		ws.scaleScratch.Scale = b.Scale
		return eval.AddQPrefixRows(b, ws.scaleScratch, b, rows)
	}

	ws.scaleScratch = ws.ensureCiphertext(params, ws.scaleScratch, b.Degree(), utils.Min(a.Level(), b.Level()))
	ws.scaleScratch.IsNTT = b.IsNTT
	ws.scaleScratch.IsMontgomery = b.IsMontgomery
	*ws.scaleScratch.MetaData = *b.MetaData
	ratio := a.Scale.Div(b.Scale).BigInt()
	if ratio.Sign() <= 0 {
		return errors.New("invalid Fast polynomial scale alignment ratio")
	}
	if err := eval.MulIntegerQPrefixRows(b, ratio, rows, ws.scaleScratch); err != nil {
		return err
	}
	ws.scaleScratch.Scale = a.Scale
	if err := eval.AddQPrefixRows(ws.scaleScratch, a, ws.scaleScratch, rows); err != nil {
		return err
	}
	copyQPrefixElement(params, ws.scaleScratch, b, rows)
	return nil
}

func copyQPrefixElement(params ckks.Parameters, src, dst *rlwe.Ciphertext, rows int) {
	fastckks.Resize(dst, src.Degree(), src.Level(), params.N())
	*dst.MetaData = *src.MetaData
	dst.IsNTT = src.IsNTT
	dst.IsMontgomery = src.IsMontgomery
	dst.Scale = src.Scale
	for d := range src.Value {
		for limb := 0; limb < rows; limb++ {
			copy(dst.Value[d].Coeffs[limb], src.Value[d].Coeffs[limb])
		}
	}
}

// cloneQPrefixResult creates an independently-owned Fast ciphertext
// without reading or materializing rows above the selected Q-prefix.
func cloneQPrefixResult(params ckks.Parameters, src *rlwe.Ciphertext, rows int) *rlwe.Ciphertext {
	dst := fastckks.NewCiphertext(params, src.Degree(), src.Level())
	copyQPrefixElement(params, src, dst, rows)
	return dst
}

func zeroQPrefix(ct *rlwe.Ciphertext, rows int) {
	for d := range ct.Value {
		for limb := 0; limb < rows; limb++ {
			ring.ZeroVec(ct.Value[d].Coeffs[limb])
		}
	}
}
