package polynomial

import (
	"errors"
	"fmt"
	"sort"

	commonpolynomial "github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// DiagnosticPower is an independently-owned snapshot of one generated Fast
// Chebyshev power. It exposes only the maintained q0/q1 ciphertext state for
// explicit experiment diagnostics; it is not part of the evaluation hot path.
type DiagnosticPower struct {
	N          int
	Ciphertext *rlwe.Ciphertext
}

// DiagnosticPlanBlock records one simulated Paterson-Stockmeyer block.
type DiagnosticPlanBlock struct {
	Degree int
	Level  int
	Scale  rlwe.Scale
}

// DiagnosticPlan records the public planner metadata used by Fast evaluation.
type DiagnosticPlan struct {
	Degree int
	Base   int
	Level  int
	Scale  rlwe.Scale
	Blocks []DiagnosticPlanBlock
}

// DiagnosticGeneratePowers runs the same Fast power-generation schedule used
// by Evaluate and returns maintained-state snapshots plus planner metadata.
// It is an explicit reference/diagnostic boundary and never materializes
// dormant q2...qL rows or invokes Standard/full-Q arithmetic.
func (eval *FastEvaluator) DiagnosticGeneratePowers(input *rlwe.Ciphertext, p bignum.Polynomial, targetScale rlwe.Scale) ([]DiagnosticPower, DiagnosticPlan, error) {
	if eval == nil || eval.Evaluator == nil {
		return nil, DiagnosticPlan{}, errors.New("Fast polynomial evaluator cannot be nil")
	}
	if input == nil || input.MetaData == nil {
		return nil, DiagnosticPlan{}, errors.New("Fast polynomial diagnostic input cannot be nil")
	}
	if len(p.Coeffs) == 0 || p.Basis != bignum.Chebyshev {
		return nil, DiagnosticPlan{}, errors.New("Fast polynomial diagnostic requires a non-empty Chebyshev polynomial")
	}
	commonPoly := commonpolynomial.NewPolynomial(p)
	levelsConsumed := eval.Parameters.LevelsConsumedPerRescaling()
	if input.Level() < levelsConsumed*commonPoly.Depth() {
		return nil, DiagnosticPlan{}, fmt.Errorf("%d levels < %d log(d) -> cannot diagnose poly", input.Level(), levelsConsumed*commonPoly.Depth())
	}

	ws := &eval.workspace
	ws.reset(eval.Parameters, input)
	if err := ws.generatePowers(eval.Parameters, eval.Evaluator, p, commonPoly); err != nil {
		return nil, DiagnosticPlan{}, err
	}

	keys := make([]int, 0, len(ws.powers))
	for n := range ws.powers {
		keys = append(keys, n)
	}
	sort.Ints(keys)
	powers := make([]DiagnosticPower, 0, len(keys))
	for _, n := range keys {
		powers = append(powers, DiagnosticPower{N: n, Ciphertext: cloneMaintainedResult(eval.Parameters, ws.powers[n])})
	}

	sim := simEvaluator{params: eval.Parameters, levelsConsumedPerRescaling: levelsConsumed}
	plan := commonPoly.PatersonStockmeyerPolynomial(eval.Evaluator.GetRLWEParameters(), input.Level(), input.Scale, targetScale, sim)
	diagnosticPlan := DiagnosticPlan{Degree: plan.Degree, Base: plan.Base, Level: plan.Level, Scale: plan.Scale, Blocks: make([]DiagnosticPlanBlock, len(plan.Value))}
	for i, block := range plan.Value {
		diagnosticPlan.Blocks[i] = DiagnosticPlanBlock{Degree: block.Degree(), Level: block.Level, Scale: block.Scale}
	}
	return powers, diagnosticPlan, nil
}
