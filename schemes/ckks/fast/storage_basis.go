package fast

import (
	"errors"
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// This fixed private storage basis is independent of every CKKS logical Q chain.
const (
	fastStoragePrime0 uint64 = 1152921504606584833
	fastStoragePrime1 uint64 = 1152921504598720513
	fastStoragePrime2 uint64 = 1152921504592429057
)

func fastStoragePrimes() [3]uint64 {
	return [3]uint64{fastStoragePrime0, fastStoragePrime1, fastStoragePrime2}
}

// storageUint192 is fixed-width arithmetic for storage products and values.
type storageUint192 struct{ lo, mid, hi uint64 }

// storageInteger is a signed-magnitude value represented with 192 bits.
type storageInteger struct {
	magnitude storageUint192
	negative  bool
}

// fastStorageBasis holds three private NTT subrings. Width is always explicit
// and is never inferred from a CKKS logical level.
type fastStorageBasis struct {
	logN       int
	subrings   [3]*ring.SubRing
	products   [4]storageUint192
	halves     [4]storageUint192
	inverse01  uint64
	inverse012 uint64
}

func newFastStorageBasis(logN int) (basis fastStorageBasis, err error) {
	if logN < 3 || logN > 16 {
		return basis, fmt.Errorf("unsupported storage basis LogN: %d", logN)
	}
	basis.logN = logN
	degree := 1 << uint(logN)
	primes := fastStoragePrimes()
	for i, q := range primes {
		if !ring.IsPrime(q) || q >= 1<<60 || q%(1<<18) != 1 {
			return fastStorageBasis{}, fmt.Errorf("invalid storage prime at index %d", i)
		}
		for j := 0; j < i; j++ {
			if q == primes[j] {
				return fastStorageBasis{}, errors.New("storage primes must be distinct")
			}
		}
	}
	storageRing, err := ring.NewRing(degree, primes[:])
	if err != nil {
		return fastStorageBasis{}, fmt.Errorf("storage ring: %w", err)
	}
	copy(basis.subrings[:], storageRing.SubRings)

	basis.products[0] = storageUint192{lo: 1}
	for width := 1; width <= len(primes); width++ {
		basis.products[width] = storageMul64(basis.products[width-1], primes[width-1])
		basis.halves[width] = storageHalf(basis.products[width])
	}
	var ok bool
	basis.inverse01, ok = storageInverseMod(primes[0]%primes[1], primes[1])
	if !ok {
		return fastStorageBasis{}, errors.New("storage primes are not pairwise coprime")
	}
	product01Mod2 := storageMod192(basis.products[2], primes[2])
	basis.inverse012, ok = storageInverseMod(product01Mod2, primes[2])
	if !ok {
		return fastStorageBasis{}, errors.New("storage basis products are not pairwise coprime")
	}
	return basis, nil
}

func (basis fastStorageBasis) subring(index int) (*ring.SubRing, error) {
	if index < 0 || index >= len(basis.subrings) || basis.subrings[index] == nil {
		return nil, fmt.Errorf("invalid storage prime index %d", index)
	}
	return basis.subrings[index], nil
}

func (basis fastStorageBasis) product(width int) (*big.Int, error) {
	if err := validateStorageWidth(width); err != nil {
		return nil, err
	}
	return storageToBig(basis.products[width]), nil
}

func (basis fastStorageBasis) productBitLen(width int) (int, error) {
	if err := validateStorageWidth(width); err != nil {
		return 0, err
	}
	return storageToBig(basis.products[width]).BitLen(), nil
}

func (basis fastStorageBasis) prime(index int) (uint64, error) {
	if index < 0 || index >= len(fastStoragePrimes()) {
		return 0, fmt.Errorf("invalid storage prime index %d", index)
	}
	return fastStoragePrimes()[index], nil
}

// centeredCapacity is floor((S-1)/2), the maximum magnitude satisfying 2|X|<S.
func (basis fastStorageBasis) centeredCapacity(width int) (*big.Int, error) {
	if err := validateStorageWidth(width); err != nil {
		return nil, err
	}
	return new(big.Int).Rsh(new(big.Int).Sub(storageToBig(basis.products[width]), big.NewInt(1)), 1), nil
}

func (basis fastStorageBasis) hasUniqueCenteredRepresentation(value *big.Int, width int) (bool, error) {
	if err := validateStorageWidth(width); err != nil {
		return false, err
	}
	if value == nil {
		return false, errors.New("value cannot be nil")
	}
	doubled := new(big.Int).Lsh(new(big.Int).Abs(new(big.Int).Set(value)), 1)
	return doubled.Cmp(storageToBig(basis.products[width])) < 0, nil
}

func (basis fastStorageBasis) encode(value *big.Int, width int) ([3]uint64, error) {
	var residues [3]uint64
	if err := validateStorageWidth(width); err != nil {
		return residues, err
	}
	if value == nil {
		return residues, errors.New("value cannot be nil")
	}
	magnitude := new(big.Int).Abs(new(big.Int).Set(value))
	if magnitude.BitLen() > 192 {
		return residues, errors.New("value exceeds the fixed 192-bit storage representation")
	}
	integer := storageInteger{magnitude: storageFromBig(magnitude), negative: value.Sign() < 0}
	return basis.encodeFixed(integer, width)
}

func (basis fastStorageBasis) encodeFixed(value storageInteger, width int) ([3]uint64, error) {
	var residues [3]uint64
	if err := validateStorageWidth(width); err != nil {
		return residues, err
	}
	primes := fastStoragePrimes()
	for i := 0; i < width; i++ {
		residue := storageMod192(value.magnitude, primes[i])
		if value.negative && residue != 0 {
			residue = primes[i] - residue
		}
		residues[i] = residue
	}
	return residues, nil
}

func (basis fastStorageBasis) decode(residues [3]uint64, width int) (*big.Int, error) {
	value, err := basis.decodeFixed(residues, width)
	if err != nil {
		return nil, err
	}
	result := storageToBig(value.magnitude)
	if value.negative {
		result.Neg(result)
	}
	return result, nil
}

func (basis fastStorageBasis) decodeFixed(residues [3]uint64, width int) (storageInteger, error) {
	if err := validateStorageWidth(width); err != nil {
		return storageInteger{}, err
	}
	primes := fastStoragePrimes()
	for i := 0; i < width; i++ {
		if residues[i] >= primes[i] {
			return storageInteger{}, fmt.Errorf("residue %d is not canonical", i)
		}
	}
	value := storageUint192{lo: residues[0]}
	if width >= 2 {
		q0, q1 := primes[0], primes[1]
		delta := storageSubMod(residues[1], value.lo%q1, q1)
		t := storageMulMod(delta, basis.inverse01, q1)
		value = storageAdd(value, storageMul64(storageUint192{lo: q0}, t))
	}
	if width == 3 {
		q2 := primes[2]
		delta := storageSubMod(residues[2], storageMod192(value, q2), q2)
		t := storageMulMod(delta, basis.inverse012, q2)
		value = storageAdd(value, storageMul64(basis.products[2], t))
	}
	if storageCmp(value, basis.halves[width]) > 0 {
		return storageInteger{magnitude: storageSub(basis.products[width], value), negative: true}, nil
	}
	return storageInteger{magnitude: value}, nil
}

func validateStorageWidth(width int) error {
	if width < 1 || width > 3 {
		return fmt.Errorf("storage width must be 1, 2, or 3: got %d", width)
	}
	return nil
}

func storageAdd(a, b storageUint192) storageUint192 {
	lo, carry := bits.Add64(a.lo, b.lo, 0)
	mid, carry := bits.Add64(a.mid, b.mid, carry)
	hi, _ := bits.Add64(a.hi, b.hi, carry)
	return storageUint192{lo: lo, mid: mid, hi: hi}
}

func storageSub(a, b storageUint192) storageUint192 {
	lo, borrow := bits.Sub64(a.lo, b.lo, 0)
	mid, borrow := bits.Sub64(a.mid, b.mid, borrow)
	hi, _ := bits.Sub64(a.hi, b.hi, borrow)
	return storageUint192{lo: lo, mid: mid, hi: hi}
}

func storageCmp(a, b storageUint192) int {
	if a.hi != b.hi {
		if a.hi < b.hi {
			return -1
		}
		return 1
	}
	if a.mid != b.mid {
		if a.mid < b.mid {
			return -1
		}
		return 1
	}
	if a.lo < b.lo {
		return -1
	}
	if a.lo > b.lo {
		return 1
	}
	return 0
}

func storageHalf(value storageUint192) storageUint192 {
	return storageUint192{lo: (value.lo >> 1) | (value.mid << 63), mid: (value.mid >> 1) | (value.hi << 63), hi: value.hi >> 1}
}

func storageMul64(value storageUint192, scalar uint64) storageUint192 {
	hi0, lo0 := bits.Mul64(value.lo, scalar)
	hi1, lo1 := bits.Mul64(value.mid, scalar)
	hi2, lo2 := bits.Mul64(value.hi, scalar)
	mid, carry := bits.Add64(hi0, lo1, 0)
	hi, carry := bits.Add64(hi1, lo2, carry)
	hi, _ = bits.Add64(hi, hi2, carry)
	return storageUint192{lo: lo0, mid: mid, hi: hi}
}

func storageMod192(value storageUint192, modulus uint64) uint64 {
	_, remainder := bits.Div64(0, value.hi%modulus, modulus)
	_, remainder = bits.Div64(remainder, value.mid, modulus)
	_, remainder = bits.Div64(remainder, value.lo, modulus)
	return remainder
}

func storageSubMod(a, b, modulus uint64) uint64 {
	if a >= b {
		return a - b
	}
	return modulus - (b - a)
}

func storageMulMod(a, b, modulus uint64) uint64 {
	hi, lo := bits.Mul64(a, b)
	_, remainder := bits.Div64(hi, lo, modulus)
	return remainder
}

func storageInverseMod(a, modulus uint64) (uint64, bool) {
	inverse := new(big.Int).ModInverse(new(big.Int).SetUint64(a), new(big.Int).SetUint64(modulus))
	if inverse == nil {
		return 0, false
	}
	return inverse.Uint64(), true
}

func storageToBig(value storageUint192) *big.Int {
	result := new(big.Int).SetUint64(value.hi)
	result.Lsh(result, 64).Or(result, new(big.Int).SetUint64(value.mid))
	result.Lsh(result, 64).Or(result, new(big.Int).SetUint64(value.lo))
	return result
}

func storageFromBig(value *big.Int) storageUint192 {
	copy := new(big.Int).Set(value)
	mask := new(big.Int).SetUint64(^uint64(0))
	lo := new(big.Int).And(copy, mask).Uint64()
	copy.Rsh(copy, 64)
	mid := new(big.Int).And(copy, mask).Uint64()
	hi := new(big.Int).Rsh(copy, 64).Uint64()
	return storageUint192{lo: lo, mid: mid, hi: hi}
}
