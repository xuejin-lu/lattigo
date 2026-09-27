package fast

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// QPrefixCapacitySnapshot is compact exact coefficient-capacity evidence for
// one named ciphertext boundary. MaxAbs contains one decimal integer per
// ciphertext component; PrefixProduct is the actual product q0*...*q(rows-1).
type QPrefixCapacitySnapshot struct {
	Name          string
	Level         int
	Rows          int
	Scale         string
	Degree        int
	MaxAbs        []string
	PrefixProduct string
	StrictFit     bool
}

// SetQPrefixCapacityObserver installs an optional diagnostic observer. A nil
// observer disables all coefficient reconstruction and leaves the evaluation
// hot path unchanged. Evaluators are single-stream and not concurrency-safe.
func (eval *Evaluator) SetQPrefixCapacityObserver(observer func(QPrefixCapacitySnapshot) error) {
	if eval != nil {
		eval.qPrefixCapacityObserver = observer
	}
}

// ObserveQPrefixCapacity records a named exact centered-coefficient bound
// when an observer is installed, and rejects the first strict-capacity failure.
func (eval *Evaluator) ObserveQPrefixCapacity(name string, ct *rlwe.Ciphertext, rows int) error {
	if eval == nil || eval.qPrefixCapacityObserver == nil {
		return nil
	}
	if name == "" || ct == nil || ct.MetaData == nil {
		return errors.New("Q-prefix capacity checkpoint requires a name and ciphertext")
	}
	if ct.N() != eval.Parameters.N() || ct.Level() < 0 {
		return errors.New("Q-prefix capacity checkpoint ciphertext does not match evaluator parameters")
	}
	if err := eval.validateExplicitRows(ct.Level(), rows); err != nil {
		return err
	}
	ringQ := eval.Parameters.RingQ()
	scratch := &eval.rescaleScratch
	if err := scratch.validateWidth(rows); err != nil {
		return err
	}
	maxima := make([]uint192, len(ct.Value))
	for component := range ct.Value {
		if err := prefixToCoefficientRows(ringQ, ct.Value[component], rows, ct.IsNTT, ct.IsMontgomery, scratch.coeff); err != nil {
			return fmt.Errorf("Q-prefix capacity checkpoint %s component %d: %w", name, component, err)
		}
		for coefficient := 0; coefficient < ringQ.N(); coefficient++ {
			var residues [MaxQPrefixWidth]uint64
			for row := 0; row < rows; row++ {
				residues[row] = scratch.coeff.Coeffs[row][coefficient]
			}
			value := reconstructQPrefix(rows, residues, scratch)
			magnitude, _ := centeredQPrefix(value, scratch.modulus[rows-1], scratch.half[rows-1])
			if cmp192(magnitude, maxima[component]) > 0 {
				maxima[component] = magnitude
			}
		}
	}

	prefixProduct := uint192ToBigInt(scratch.modulus[rows-1])
	snapshot := QPrefixCapacitySnapshot{
		Name:          name,
		Level:         ct.Level(),
		Rows:          rows,
		Scale:         ct.Scale.Value.Text('g', -1),
		Degree:        ct.Degree(),
		MaxAbs:        make([]string, len(maxima)),
		PrefixProduct: prefixProduct.String(),
		StrictFit:     true,
	}
	for component, maximum := range maxima {
		bound := uint192ToBigInt(maximum)
		snapshot.MaxAbs[component] = bound.String()
		if !QPrefixCapacitySatisfied(bound, prefixProduct) {
			snapshot.StrictFit = false
		}
	}
	if err := eval.qPrefixCapacityObserver(snapshot); err != nil {
		return err
	}
	if !snapshot.StrictFit {
		for component, maximum := range maxima {
			bound := uint192ToBigInt(maximum)
			if !QPrefixCapacitySatisfied(bound, prefixProduct) {
				return &QPrefixCapacityError{Level: ct.Level(), Component: component, Bound: bound, PrefixProduct: new(big.Int).Set(prefixProduct)}
			}
		}
	}
	return nil
}
