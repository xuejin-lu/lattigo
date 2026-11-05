// This test is compiled beside byte-identical original CNN512 Go files in an
// isolated temporary directory. It changes no CNN source and runs no Bootstrap.
package main

import (
	"math"
	"testing"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestCP2ActualCNNCompactQPChain(t *testing.T) {
	logQ := []int{42}
	for i := 0; i < 15; i++ {
		logQ = append(logQ, 32)
	}
	for i := 0; i < 11; i++ {
		logQ = append(logQ, 42)
	}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 8, LogQ: logQ, LogP: []int{42, 42, 42, 42, 42}, LogDefaultScale: 32,
	})
	if err != nil {
		t.Fatal(err)
	}
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	slots := params.MaxSlots()
	values, mask, bias := make([]float64, slots), make([]float64, slots), make([]float64, slots)
	for i := range values {
		values[i] = 0.2 * math.Sin(float64(i)/9)
		mask[i] = 0.25
		bias[i] = 0.01
	}
	pt := ckks.NewPlaintext(params, 15)
	if err = encoder.Encode(values, pt); err != nil {
		t.Fatal(err)
	}
	input, err := rlwe.NewEncryptor(params, kgen.GenPublicKeyNew(sk)).EncryptNew(pt)
	if err != nil {
		t.Fatal(err)
	}
	exec := NewConvExecutor(params, sk, ConvExecConfig{Slots: slots, KeyLevelQ: 15, KeyChunk: 1, Workers: 1})
	first, err := exec.ExecMaskedRot("cp2.front", []*rlwe.Ciphertext{input}, 1, []mrTerm{{Src: 0, Out: 0, R: 1, G: 0, Mask: mask}})
	if err != nil {
		t.Fatal(err)
	}
	if err = exec.Square(first); err != nil {
		t.Fatal(err)
	}
	if err = exec.AddVecs(first, [][]float64{bias}); err != nil {
		t.Fatal(err)
	}
	want := make([]float64, slots)
	for i := range want {
		v := 0.25 * values[(i+1)%slots]
		want[i] = v*v + 0.01
	}
	check := func(t *testing.T, name string, ct *rlwe.Ciphertext, expected []float64, level int, scale rlwe.Scale) {
		t.Helper()
		if ct.Level() != level || !ct.Scale.Equal(scale) {
			t.Fatalf("%s: incorrect Level/Scale", name)
		}
		for _, poly := range ct.Value {
			for row := 4; row <= level; row++ {
				if len(poly.Coeffs[row]) != 0 {
					t.Fatalf("%s: dormant q%d was materialized", name, row)
				}
			}
		}
		decoded := make([]float64, slots)
		if err := encoder.Decode(decryptor.DecryptNew(ct), decoded); err != nil {
			t.Fatal(err)
		}
		maxError := 0.0
		for i, v := range decoded {
			if math.IsNaN(v) || math.IsInf(v, 0) {
				t.Fatalf("%s: nonfinite slot %d", name, i)
			}
			maxError = math.Max(maxError, math.Abs(v-expected[i]))
		}
		t.Logf("%s: max_error=%.12g level=%d", name, maxError, level)
		if maxError > 1e-6 {
			t.Fatalf("%s: numerical oracle failed: %g", name, maxError)
		}
	}
	// Square consumes logical q14; the original frontend does not relabel its
	// scale to Delta here. Check the actual public multiplication/rescale rule.
	frontScale := params.DefaultScale().Mul(params.DefaultScale()).Div(rlwe.NewScale(params.Q()[14]))
	check(t, "front conv-square-rescale-bias", first[0], want, 13, frontScale)
	before := first[0].CopyNew()
	cases := []struct {
		name      string
		rotations [][2]int
	}{
		{"r0", [][2]int{{0, 0}}},
		{"baby", [][2]int{{1, 0}}},
		{"giant", [][2]int{{0, 5}}},
		{"bsgs-cross-evaluator", [][2]int{{1, 5}, {0, 9}, {2, 0}, {2, 5}}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			terms := make([]mrTerm, len(tc.rotations))
			expected := make([]float64, slots)
			for k, rg := range tc.rotations {
				terms[k] = mrTerm{Src: 0, Out: 0, R: rg[0], G: rg[1], Mask: mask}
				for i := range expected {
					expected[i] += 0.25 * want[(i+rg[0]+rg[1])%slots]
				}
			}
			out, err := exec.ExecMaskedRot("cp2."+tc.name, first, 1, terms)
			if err != nil {
				t.Fatalf("numerical acceptance requires success, got: %v", err)
			}
			check(t, tc.name, out[0], expected, 12, params.DefaultScale())
			if !first[0].Equal(before) {
				t.Fatal("convolution changed its source")
			}
		})
	}
}
