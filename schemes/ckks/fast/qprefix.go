package fast

import (
	"errors"
	"fmt"
	"math/big"
)

// MaxQPrefixWidth is the initial Q-prefix v2 engineering cap. It is not a
// claim that four residues are sufficient for every possible CKKS operation.
const MaxQPrefixWidth = 4

// QPrefixWidth returns the number of logical-Q residues maintained at level.
// The logical Level remains authoritative; this policy only caps the active
// prefix at q0..q3.
func QPrefixWidth(level int) (int, error) {
	if level < 0 {
		return 0, fmt.Errorf("Q-prefix level must be non-negative: %d", level)
	}
	if level >= MaxQPrefixWidth-1 {
		return MaxQPrefixWidth, nil
	}
	return level + 1, nil
}

func qPrefixWidthOrPanic(level int) int {
	width, err := QPrefixWidth(level)
	if err != nil {
		panic(err)
	}
	return width
}

// QPrefixProduct returns the exact product of the actual q values in the
// maintained prefix for level. It does not infer capacity from prime bit sizes.
func QPrefixProduct(q []uint64, level int) (*big.Int, error) {
	width, err := QPrefixWidth(level)
	if err != nil {
		return nil, err
	}
	if level >= len(q) || len(q) < width {
		return nil, fmt.Errorf("Q-prefix level %d requires q values through q%d; got %d moduli", level, level, len(q))
	}

	product := big.NewInt(1)
	for i := 0; i < width; i++ {
		if q[i] == 0 {
			return nil, fmt.Errorf("Q-prefix modulus q%d must be positive", i)
		}
		product.Mul(product, new(big.Int).SetUint64(q[i]))
	}
	return product, nil
}

// QPrefixCapacityError identifies the first component whose proven absolute
// coefficient bound does not fit the strict centered-capacity contract.
// Bound and PrefixProduct are independent copies owned by the error.
type QPrefixCapacityError struct {
	Level         int
	Component     int
	Bound         *big.Int
	PrefixProduct *big.Int
}

func (err *QPrefixCapacityError) Error() string {
	if err == nil {
		return "Q-prefix capacity exceeded"
	}
	return fmt.Sprintf("Q-prefix capacity exceeded at Level %d component %d: strict 2B < S_Q failed (2B=%s, S_Q=%s)",
		err.Level, err.Component, doubledBoundString(err.Bound), bigIntString(err.PrefixProduct))
}

func doubledBoundString(bound *big.Int) string {
	if bound == nil {
		return "<nil>"
	}
	return new(big.Int).Lsh(new(big.Int).Set(bound), 1).String()
}

func bigIntString(value *big.Int) string {
	if value == nil {
		return "<nil>"
	}
	return value.String()
}

// QPrefixCapacitySatisfied reports whether a non-negative component bound is
// strictly inside the centered uniqueness interval: 2*bound < prefixProduct.
func QPrefixCapacitySatisfied(bound, prefixProduct *big.Int) bool {
	if bound == nil || prefixProduct == nil || bound.Sign() < 0 || prefixProduct.Sign() <= 0 {
		return false
	}
	twiceBound := new(big.Int).Lsh(new(big.Int).Set(bound), 1)
	return twiceBound.Cmp(prefixProduct) < 0
}

// CheckQPrefixCapacity validates each component bound against the actual
// maintained prefix product. It is a pure preflight: neither q nor any bound
// is modified, so callers can reject before committing operation outputs.
func CheckQPrefixCapacity(level int, q []uint64, componentBounds []*big.Int) error {
	prefixProduct, err := QPrefixProduct(q, level)
	if err != nil {
		return err
	}
	for component, bound := range componentBounds {
		if bound == nil {
			return fmt.Errorf("Q-prefix component %d bound cannot be nil", component)
		}
		if bound.Sign() < 0 {
			return fmt.Errorf("Q-prefix component %d bound cannot be negative", component)
		}
		if !QPrefixCapacitySatisfied(bound, prefixProduct) {
			return &QPrefixCapacityError{
				Level:         level,
				Component:     component,
				Bound:         new(big.Int).Set(bound),
				PrefixProduct: new(big.Int).Set(prefixProduct),
			}
		}
	}
	return nil
}

// CommitQPrefixBounds validates all candidate component bounds before
// replacing dst. On any error, including a capacity error, *dst is unchanged.
// The committed values are deep copies and do not alias componentBounds.
func CommitQPrefixBounds(dst *[]*big.Int, level int, q []uint64, componentBounds []*big.Int) error {
	if dst == nil {
		return errors.New("Q-prefix bounds destination cannot be nil")
	}
	if err := CheckQPrefixCapacity(level, q, componentBounds); err != nil {
		return err
	}
	staged := make([]*big.Int, len(componentBounds))
	for i, bound := range componentBounds {
		staged[i] = new(big.Int).Set(bound)
	}
	*dst = staged
	return nil
}
