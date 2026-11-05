package ring

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

func compactFixture(n, level, rows int) Poly {
	p := NewPoly(n, level)
	for i := range p.Coeffs {
		if i >= rows {
			p.Coeffs[i] = nil
		} else {
			for j := range p.Coeffs[i] {
				p.Coeffs[i][j] = uint64(11*i + j + 1)
			}
		}
	}
	return p
}

func TestCompactCopyLifecycle(t *testing.T) {
	src := compactFixture(16, 5, 2)
	src.Coeffs[4] = []uint64{} // nil and empty must agree.
	dst := compactFixture(16, 5, 6)
	shallow := dst
	oldHigh := append([]uint64(nil), shallow.Coeffs[4]...)
	dst.CopyLvl(5, src)
	for i := 0; i < 2; i++ {
		require.Equal(t, src.Coeffs[i], dst.Coeffs[i])
	}
	for i := 2; i < 6; i++ {
		require.Empty(t, dst.Coeffs[i])
	}
	// Changing destination availability must not destroy another header owner.
	require.Equal(t, oldHigh, shallow.Coeffs[4])
	clone := dst.CopyNew()
	dst.Zero()
	for _, row := range dst.Coeffs[:2] {
		require.Equal(t, make([]uint64, 16), row)
	}
	for _, row := range dst.Coeffs[2:] {
		require.Empty(t, row)
	}
	require.Equal(t, src.Coeffs[0], clone.Coeffs[0])
	dst.Copy(src)
	dst.CopyLvl(5, dst)
	full := compactFixture(16, 5, 6)
	dst.Copy(full)
	require.Equal(t, full, dst)
	dst.Coeffs[0][0]++
	require.NotEqual(t, full.Coeffs[0][0], dst.Coeffs[0][0])
	clone.Resize(8)
	for _, row := range clone.Coeffs[2:] {
		require.Empty(t, row)
	}
	full.Resize(8)
	for _, row := range full.Coeffs[6:] {
		require.Len(t, row, 16)
	}
}

func TestCompactCopyRejectsBeforeMutationAndSupportsMemmove(t *testing.T) {
	for _, malformed := range []int{1, 2} {
		src := compactFixture(16, 3, 4)
		if malformed == 1 {
			src.Coeffs[2] = nil
		} else {
			src.Coeffs[2] = src.Coeffs[2][:15]
		}
		dst := compactFixture(16, 4, 5)
		before := dst.CopyNew()
		require.Panics(t, func() { dst.Copy(src) })
		require.Equal(t, *before, dst)
		invalidBefore := src.Coeffs.CopyNew()
		for i := range invalidBefore {
			if src.Coeffs[i] == nil {
				invalidBefore[i] = nil
			}
		}
		require.Panics(t, func() { src.Zero() })
		require.Panics(t, func() { src.CopyNew() })
		for i, row := range invalidBefore {
			require.Equal(t, row, src.Coeffs[i])
		}
	}
	backing := make([]uint64, 17)
	for i := range backing {
		backing[i] = uint64(i)
	}
	a := Poly{Coeffs: [][]uint64{backing[:16]}}
	b := Poly{Coeffs: [][]uint64{backing[1:]}}
	b.CopyLvl(0, a)
	for i, v := range b.Coeffs[0] {
		require.Equal(t, uint64(i), v)
	}
	a = compactFixture(16, 2, 3)
	b = compactFixture(16, 2, 3)
	b.Coeffs[1] = a.Coeffs[0]
	before := b.CopyNew()
	require.Panics(t, func() { b.CopyLvl(2, a) })
	require.Equal(t, *before, b)
}

type compactTestPool struct{ rows []*[]uint64 }

func (p *compactTestPool) Get() *[]uint64 {
	if len(p.rows) == 0 {
		row := make([]uint64, 16)
		return &row
	}
	last := p.rows[len(p.rows)-1]
	p.rows = p.rows[:len(p.rows)-1]
	return last
}
func (p *compactTestPool) Put(row *[]uint64) { p.rows = append(p.rows, row) }

func TestCompactPoolOwnershipAndReuse(t *testing.T) {
	pool := &compactTestPool{}
	var _ structs.BufferPool[*[]uint64] = pool
	for iteration := 0; iteration < 8; iteration++ {
		p := NewPolyFromUintPool(pool, 4)
		for _, row := range p.Coeffs {
			require.Len(t, row, 16)
		}
		p.CopyLvl(4, compactFixture(16, 4, 2))
		retainedHeaders := p.Coeffs
		RecyclePolyInUintPool(pool, p)
		require.Empty(t, p.Coeffs)
		// Recycling must not leave pool pointers aimed at mutable old header slots.
		retainedHeaders[0] = nil
		for _, row := range pool.rows {
			require.Len(t, *row, 16)
		}
		p = NewPolyFromUintPool(pool, 4)
		p.Zero() // pool storage is not implicitly initialized input.
		for _, row := range p.Coeffs {
			require.Equal(t, make([]uint64, 16), row)
		}
		RecyclePolyInUintPool(pool, p)
	}
	p := compactFixture(16, 2, 3)
	p.Coeffs[1] = p.Coeffs[0]
	before := len(pool.rows)
	require.Panics(t, func() { RecyclePolyInUintPool(pool, &p) })
	require.Len(t, pool.rows, before)
}
