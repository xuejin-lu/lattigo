// Package fastcore contains import-neutral CKKS ring kernels shared by the
// ordinary public evaluator's bounded Fast path and the explicit Q-prefix
// evaluator. It must not depend on either parent ckks or child fast packages.
package fastcore

import (
	"errors"
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// PolynomialPair describes one input/output polynomial pair for a shared
// automorphism. ApplyRows validates every pair before mutating any output.
type PolynomialPair struct {
	Input  ring.Poly
	Output ring.Poly
}

// AutomorphismCore is the narrow dependency shared by CKKS evaluator facades.
type AutomorphismCore interface {
	ApplyRows(ringQ *ring.Ring, galEl uint64, isNTT bool, rows int, pairs ...PolynomialPair) error
}

// AutomorphismWorkspace owns per-evaluator permutation cache and row scratch.
// Evaluators sharing one workspace must remain single-stream, as before.
type AutomorphismWorkspace struct {
	n          int
	maxRows    int
	indexCache map[uint64][]uint64
	scratch    [][]uint64
}

// NewAutomorphismWorkspace allocates evaluator-owned storage for at most
// maxRows residues. Invalid dimensions leave an unusable workspace that fails
// closed when ApplyRows is called.
func NewAutomorphismWorkspace(n, maxRows int) *AutomorphismWorkspace {
	w := &AutomorphismWorkspace{n: n, maxRows: maxRows}
	if n < 1 || maxRows < 1 {
		return w
	}
	w.indexCache = make(map[uint64][]uint64)
	w.scratch = make([][]uint64, maxRows)
	for row := range w.scratch {
		w.scratch[row] = make([]uint64, n)
	}
	return w
}

// ScratchRows reports the allocated workspace width for internal regression
// tests; it does not imply that every allocated ciphertext row is authoritative.
func (w *AutomorphismWorkspace) ScratchRows() int {
	if w == nil {
		return 0
	}
	return len(w.scratch)
}

// CachedGaloisElementCount reports the number of evaluator-local NTT index
// tables retained by this workspace.
func (w *AutomorphismWorkspace) CachedGaloisElementCount() int {
	if w == nil {
		return 0
	}
	return len(w.indexCache)
}

// ApplyRows applies the existing Fast coefficient/NTT permutation to exactly
// rows of every pair. Rows above the requested prefix are not read or written.
// The explicit row bound is logical Level+1 rather than the Fast Q-prefix cap;
// callers enforce their own representation policy before calling this core.
func (w *AutomorphismWorkspace) ApplyRows(ringQ *ring.Ring, galEl uint64, isNTT bool, rows int, pairs ...PolynomialPair) error {
	if w == nil || w.n < 1 || w.maxRows < 1 {
		return errors.New("Fast automorphism workspace is not initialized")
	}
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if ringQ.Type() != ring.Standard {
		return fmt.Errorf("Fast automorphism requires the Standard ring, got %s", ringQ.Type())
	}
	if ringQ.N() != w.n {
		return fmt.Errorf("automorphism workspace dimension %d does not match ring dimension %d", w.n, ringQ.N())
	}
	if rows < 1 || rows > ringQ.Level()+1 || rows > w.maxRows || rows > len(ringQ.SubRings) {
		return fmt.Errorf("automorphism row count %d must be in [1,%d] and fit workspace capacity %d", rows, ringQ.Level()+1, w.maxRows)
	}
	if len(pairs) == 0 {
		return errors.New("Fast automorphism requires at least one polynomial pair")
	}
	for row := 0; row < rows; row++ {
		if ringQ.SubRings[row] == nil || len(w.scratch) <= row || len(w.scratch[row]) != w.n {
			return fmt.Errorf("Fast automorphism q%d ring or scratch is unavailable", row)
		}
	}
	for pairIndex, pair := range pairs {
		for name, poly := range map[string]ring.Poly{"input": pair.Input, "output": pair.Output} {
			if len(poly.Coeffs) < rows || poly.Level() < rows-1 {
				return fmt.Errorf("polynomial pair %d %s has insufficient logical rows for requested row count %d", pairIndex, name, rows)
			}
			for row := 0; row < rows; row++ {
				if len(poly.Coeffs[row]) != w.n {
					return fmt.Errorf("polynomial pair %d %s q%d backing length %d does not match N=%d", pairIndex, name, row, len(poly.Coeffs[row]), w.n)
				}
			}
		}
	}

	var index []uint64
	if isNTT {
		var ok bool
		index, ok = w.indexCache[galEl]
		if !ok {
			var err error
			index, err = ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galEl)
			if err != nil {
				return fmt.Errorf("compute NTT automorphism index: %w", err)
			}
			w.indexCache[galEl] = index
		}
	}

	// All fallible validation and NTT index construction precede writes. The
	// shared row scratch also makes each individual input/output pair alias-safe.
	if isNTT {
		for _, pair := range pairs {
			for row := 0; row < rows; row++ {
				tmp := w.scratch[row]
				for j, src := range index {
					tmp[j] = pair.Input.Coeffs[row][src]
				}
				copy(pair.Output.Coeffs[row], tmp)
			}
		}
		return nil
	}

	mask := uint64(ringQ.N() - 1)
	logN := uint(bits.Len64(mask))
	for _, pair := range pairs {
		for row := 0; row < rows; row++ {
			modulus := ringQ.SubRings[row].Modulus
			tmp := w.scratch[row]
			for i, value := range pair.Input.Coeffs[row] {
				raw := uint64(i) * galEl
				position := raw & mask
				if (raw>>logN)&1 == 0 {
					tmp[position] = value
				} else {
					// Preserve ring.Ring.Automorphism's representation, including
					// modulus for an input zero.
					tmp[position] = modulus - value
				}
			}
			copy(pair.Output.Coeffs[row], tmp)
		}
	}
	return nil
}
