// Package compact implements storage checks shared by Ring and RingQP.
// It knows neither CKKS policy nor the mathematical provenance of residues.
package compact

import (
	"fmt"
	"unsafe"
)

// Width validates the requested view and returns its contiguous available prefix.
// Empty and nil rows are equivalent. Rows outside the view are not certified.
func Width(coeffs [][]uint64, level, n int) int {
	if level < 0 || len(coeffs) <= level {
		panic(fmt.Sprintf("ring compact: missing logical rows for level %d", level))
	}
	width := 0
	missing := false
	for _, row := range coeffs[:level+1] {
		switch len(row) {
		case 0:
			missing = true
		default:
			if len(row) != n || missing {
				panic("ring compact: malformed row or noncontiguous available prefix")
			}
			width++
		}
	}
	return width
}

// Overlap compares backing ranges, including shifted slices.
func Overlap(a, b []uint64) bool {
	if len(a) == 0 || len(b) == 0 {
		return false
	}
	x, y := uintptr(unsafe.Pointer(unsafe.SliceData(a))), uintptr(unsafe.Pointer(unsafe.SliceData(b)))
	return x < y+uintptr(len(b))*8 && y < x+uintptr(len(a))*8
}

// CheckAliases rejects cross-row/output overlap before writes. Same-row exact
// aliases are supported by elementwise kernels; Copy additionally permits memmove.
func CheckAliases(out [][]uint64, level int, exact, memmove bool, inputs ...[][]uint64) {
	for _, input := range inputs {
		if headersOverlap(out, input) && unsafe.SliceData(out) != unsafe.SliceData(input) {
			panic("ring compact: partially overlapping row headers")
		}
	}
	for i, dst := range out[:level+1] {
		for j := 0; j < i; j++ {
			if Overlap(dst, out[j]) {
				panic("ring compact: overlapping output rows")
			}
		}
		for _, input := range inputs {
			for j, src := range input[:level+1] {
				if Overlap(dst, src) && !(i == j && (memmove || exact && unsafe.SliceData(dst) == unsafe.SliceData(src))) {
					panic("ring compact: unsupported overlapping input/output")
				}
			}
		}
	}
}

func headersOverlap(a, b [][]uint64) bool {
	if len(a) == 0 || len(b) == 0 {
		return false
	}
	x, y := uintptr(unsafe.Pointer(unsafe.SliceData(a))), uintptr(unsafe.Pointer(unsafe.SliceData(b)))
	size := unsafe.Sizeof([]uint64{})
	return x < y+uintptr(len(b))*size && y < x+uintptr(len(a))*size
}

// Disjoint prevents a Q write from mutating P storage (or conversely).
func Disjoint(a, b [][]uint64) {
	if headersOverlap(a, b) {
		panic("ring compact: overlapping Q/P row headers")
	}
	for _, x := range a {
		for _, y := range b {
			if Overlap(x, y) {
				panic("ring compact: overlapping Q/P backing")
			}
		}
	}
}

// Plan checks all shapes/aliases before preparation or arithmetic can write.
// Destination availability constrains accumulation, but never an overwrite.
func Plan(n, level int, out [][]uint64, accumulate, exactAlias bool, inputs ...[][]uint64) int {
	return plan(n, level, out, accumulate, exactAlias, false, inputs...)
}

// PlanTransform also preserves the existing NTT convention allowing full
// oversized workspace rows. Compact rows still have exactly the Ring's N.
func PlanTransform(n, level int, out [][]uint64, inputs ...[][]uint64) int {
	return plan(n, level, out, false, true, true, inputs...)
}

func viewWidth(coeffs [][]uint64, level, n int) int {
	if len(coeffs) != 0 && len(coeffs[0]) > n {
		if Width(coeffs, level, len(coeffs[0])) != level+1 {
			panic("ring compact: oversized compact rows are unsupported")
		}
		return level + 1
	}
	return Width(coeffs, level, n)
}

func plan(n, level int, out [][]uint64, accumulate, exactAlias, transform bool, inputs ...[][]uint64) int {
	rows := level + 1
	outRows := viewWidth(out, level, n)
	for _, input := range inputs {
		rows = min(rows, viewWidth(input, level, n))
	}
	if accumulate {
		rows = min(rows, outRows)
	}
	CheckAliases(out, level, exactAlias, false, inputs...)
	if rows < level+1 && len(out) != level+1 {
		panic("ring compact: compact output requires matching logical level")
	}
	if rows < level+1 {
		if len(out[0]) > n {
			panic("ring compact: oversized workspace cannot become compact")
		}
		for _, poly := range inputs {
			if len(poly[0]) > n {
				panic("ring compact: oversized workspace cannot become compact")
			}
		}
	}
	if !transform && rows > 0 && len(inputs) > 0 {
		// Existing elementwise kernels iterate over the first input's length.
		needed := len(inputs[0][0])
		for _, input := range inputs {
			if len(input[0]) < needed {
				panic("ring compact: input row is shorter than elementwise kernel")
			}
		}
		if len(out[0]) != 0 && len(out[0]) < needed || len(out[0]) == 0 && needed > n {
			panic("ring compact: output row is shorter than elementwise kernel")
		}
	}
	return rows
}

// Permutation validates every index before the unsafe kernel, including alias
// with the output, which could otherwise change later indices during the loop.
func Permutation(n int, index []uint64, outputs ...[][]uint64) {
	if len(index) != n {
		panic("ring compact: permutation index length mismatch")
	}
	seen := make([]bool, n)
	for _, x := range index {
		if x >= uint64(n) || seen[x] {
			panic("ring compact: invalid permutation index")
		}
		seen[x] = true
	}
	for _, out := range outputs {
		for _, row := range out {
			if Overlap(index, row) {
				panic("ring compact: permutation index aliases output")
			}
		}
	}
}

// Prepare writes row headers only after the caller has validated all operands.
// A value copy sharing the outer slice is an in-place alias, not independent storage.
func Prepare(out [][]uint64, level, n, rows int) {
	for i := 0; i <= level; i++ {
		if i < rows {
			if len(out[i]) == 0 {
				out[i] = make([]uint64, n)
			}
		} else {
			out[i] = nil
		}
	}
}
