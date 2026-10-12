package ckks

import (
	"crypto/sha256"
	"encoding/binary"
	"fmt"
	"math/big"
	"sync"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// fastQPCompatibilityState is local to one ordinary CKKS evaluator. The
// compatibility lane is armed only after a complete high-level Q ciphertext
// has been checked against the authoritative Q-prefix CRT representative.
// Once armed, BTS's full-Q custom BSGS intermediates may use the existing
// full-RNS ModDown and Rescale implementations after backing and c1 checks.
type fastQPCompatibilityState struct {
	mu       sync.Mutex
	armed    bool
	pending  error
	trusted  map[*rlwe.Ciphertext]fastQPFingerprint
	counters fastQPCompatibilityCounters
}

type fastQPFingerprint struct {
	level int
	hash  [sha256.Size]byte
}

const maxTrustedFastQPCiphertexts = 4

type fastQPCompatibilityCounters struct {
	decomposeElisions          uint64
	applyEvaluationKeyElisions uint64
	gadgetProductElisions      uint64
	hoistedProductElisions     uint64
	keySwitchElisions          uint64
	automorphisms              uint64
	fullQRescales              uint64
}

func newFastQPCompatibilityState() *fastQPCompatibilityState {
	return &fastQPCompatibilityState{trusted: make(map[*rlwe.Ciphertext]fastQPFingerprint)}
}

// Automorphism is the ordinary CKKS public-API boundary for Fast zero-secret
// ciphertexts. The maintained Q-prefix path reuses the shared Fast
// automorphism core. When a caller (such as BTS C2S) supplies a complete
// high-level Q representation, the full-Q compatibility path is allowed only
// after proving its rows are the centered lift of the authoritative prefix;
// it performs the Galois permutation but never a key switch because c1 was
// checked zero. Standard execution keeps the embedded RLWE implementation.
func (eval Evaluator) Automorphism(ctIn *rlwe.Ciphertext, galEl uint64, opOut *rlwe.Ciphertext) error {
	if !eval.usesFastCKKSZeroSecretSimulation() {
		return eval.Evaluator.Automorphism(ctIn, galEl, opOut)
	}
	if err := eval.fastQPCompat.take(); err != nil {
		return fmt.Errorf("cannot Fast CKKS Automorphism: previous compatibility operation failed: %w", err)
	}
	level, _, err := eval.validateFastQPCompatCiphertext("input", ctIn, false)
	if err != nil {
		return err
	}
	outLevel, err := fastCKKSZeroSecretCiphertextLevel(opOut)
	if err != nil {
		return err
	}
	if opOut.MetaData == nil {
		return fmt.Errorf("cannot Fast CKKS Automorphism: output metadata cannot be nil")
	}
	if level != outLevel {
		return fmt.Errorf("cannot Fast CKKS Automorphism: input level %d and output level %d must match", level, outLevel)
	}
	if ctIn.IsNTT != opOut.IsNTT || ctIn.IsMontgomery != opOut.IsMontgomery {
		return fmt.Errorf("cannot Fast CKKS Automorphism: input and output domains/representations must match")
	}
	if ctIn.IsMontgomery && !ctIn.IsNTT {
		return fmt.Errorf("cannot Fast CKKS Automorphism: Montgomery representation requires NTT domain")
	}

	// Full-Q rows are used only when both sides have complete storage. Otherwise
	// the shared Fast core consumes and writes only the maintained prefix.
	if level+1 > fastcore.MaxQPrefixWidth &&
		fullFastPolyRows(ctIn.Value[0], level, ctIn.N()) &&
		fullFastPolyRows(opOut.Value[0], level, opOut.N()) &&
		fullFastPolyRows(opOut.Value[1], level, opOut.N()) {
		if err = eval.verifyFastFullQCiphertext(ctIn); err != nil {
			return fmt.Errorf("cannot Fast CKKS Automorphism: complete Q input is not authorized: %w", err)
		}
		params := eval.GetParameters()
		ringQ := params.RingQ().AtLevel(level)
		staged := NewCiphertext(*params, 1, level)
		staged.IsNTT = ctIn.IsNTT
		staged.IsMontgomery = ctIn.IsMontgomery
		if ctIn.IsNTT {
			ringQ.AutomorphismNTT(ctIn.Value[0], galEl, staged.Value[0])
		} else {
			ringQ.Automorphism(ctIn.Value[0], galEl, staged.Value[0])
		}
		staged.Value[1].Zero()
		*staged.MetaData = *ctIn.MetaData
		opOut.Copy(staged)
		fingerprint := fastQFingerprint(opOut, level)
		eval.fastQPCompat.trust(opOut, fingerprint)
		eval.fastQPCompat.note(func(c *fastQPCompatibilityCounters) { c.automorphisms++ })
		return nil
	}

	if err = eval.rotateFastCKKSZeroSecret(ctIn, galEl, opOut); err != nil {
		return fmt.Errorf("cannot Fast CKKS Automorphism: %w", err)
	}
	eval.fastQPCompat.note(func(c *fastQPCompatibilityCounters) { c.automorphisms++ })
	return nil
}

func (state *fastQPCompatibilityState) record(err error) {
	if state == nil || err == nil {
		return
	}
	state.mu.Lock()
	defer state.mu.Unlock()
	if state.pending == nil {
		state.pending = err
	}
}

func (state *fastQPCompatibilityState) peek() error {
	if state == nil {
		return nil
	}
	state.mu.Lock()
	defer state.mu.Unlock()
	return state.pending
}

func (state *fastQPCompatibilityState) take() error {
	if state == nil {
		return nil
	}
	state.mu.Lock()
	defer state.mu.Unlock()
	err := state.pending
	state.pending = nil
	return err
}

func (state *fastQPCompatibilityState) isArmed() bool {
	if state == nil {
		return false
	}
	state.mu.Lock()
	defer state.mu.Unlock()
	return state.armed
}

func (state *fastQPCompatibilityState) trust(ct *rlwe.Ciphertext, fingerprint fastQPFingerprint) {
	state.mu.Lock()
	if _, exists := state.trusted[ct]; !exists && len(state.trusted) >= maxTrustedFastQPCiphertexts {
		for previous := range state.trusted {
			delete(state.trusted, previous)
			break
		}
	}
	state.armed = true
	state.trusted[ct] = fingerprint
	state.mu.Unlock()
}

func (state *fastQPCompatibilityState) note(update func(*fastQPCompatibilityCounters)) {
	state.mu.Lock()
	update(&state.counters)
	state.mu.Unlock()
}

func (state *fastQPCompatibilityState) counterSnapshot() fastQPCompatibilityCounters {
	state.mu.Lock()
	defer state.mu.Unlock()
	return state.counters
}

func (state *fastQPCompatibilityState) trustedFingerprint(ct *rlwe.Ciphertext, fingerprint fastQPFingerprint) bool {
	state.mu.Lock()
	defer state.mu.Unlock()
	known, ok := state.trusted[ct]
	return ok && known == fingerprint
}

// ApplyEvaluationKey preserves the ordinary RLWE API for Fast CKKS's
// zero-secret mode. It elides key switching only after verifying that every
// authoritative c1 Q-prefix residue is zero. Ring-degree changes are outside
// this bounded BTS compatibility lane.
func (eval Evaluator) ApplyEvaluationKey(ctIn *rlwe.Ciphertext, evk *rlwe.EvaluationKey, opOut *rlwe.Ciphertext) error {
	if !eval.usesFastCKKSZeroSecretSimulation() {
		return eval.Evaluator.ApplyEvaluationKey(ctIn, evk, opOut)
	}
	if err := eval.fastQPCompat.take(); err != nil {
		return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: previous compatibility operation failed: %w", err)
	}
	if evk == nil || evk.Layout != rlwe.KeyLayoutFast || evk.GadgetCiphertext.Layout != rlwe.KeyLayoutFast {
		return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: a Fast evaluation key is required")
	}
	keyDegree, keyLevelQ, _, err := fastGadgetCiphertextShape(&evk.GadgetCiphertext, eval.GetParameters())
	if err != nil || keyDegree != 1 {
		return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: malformed degree-one Fast evaluation key")
	}
	levelIn, _, err := eval.validateFastQPCompatCiphertext("input", ctIn, false)
	if err != nil {
		return err
	}
	if opOut == nil || opOut.MetaData == nil || len(opOut.Value) != 2 {
		return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: output must have degree-one storage and metadata")
	}
	for component := range opOut.Value {
		if len(opOut.Value[component].Coeffs) == 0 || len(opOut.Value[component].Coeffs[0]) == 0 {
			return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: output component %d has no backing", component)
		}
	}
	if ctIn.N() != opOut.N() || ctIn.N() != eval.GetParameters().N() {
		return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: ring-degree conversion is not supported by the Fast Q-prefix path")
	}
	level := utils.Min(levelIn, opOut.Level())
	if level > keyLevelQ {
		return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: input level %d exceeds Fast evaluation-key level %d", level, keyLevelQ)
	}
	rows, err := fastcore.QPrefixWidth(level)
	if err != nil {
		return err
	}
	if err = fastcore.ValidatePrefixRows(eval.GetParameters().RingQ(), level, rows, ctIn.Value[0], ctIn.Value[1]); err != nil {
		return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: %w", err)
	}
	if !fastQPolyIsZero(ctIn.Value[1], rows) {
		return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: input c1 is nonzero on the authoritative Q-prefix")
	}
	if opOut == ctIn && level == levelIn {
		eval.fastQPCompat.note(func(c *fastQPCompatibilityCounters) { c.applyEvaluationKeyElisions++ })
		return nil
	}
	if err = fastcore.ValidatePrefixBacking(eval.GetParameters().RingQ(), rows, opOut.Value[0], opOut.Value[1]); err != nil {
		return fmt.Errorf("cannot Fast CKKS ApplyEvaluationKey: output: %w", err)
	}
	fastcore.ResizeCompactCiphertext(opOut, 1, level, eval.GetParameters().N())
	for row := 0; row < rows; row++ {
		copy(opOut.Value[0].Coeffs[row], ctIn.Value[0].Coeffs[row])
		clear(opOut.Value[1].Coeffs[row])
	}
	*opOut.MetaData = *ctIn.MetaData
	eval.fastQPCompat.note(func(c *fastQPCompatibilityCounters) { c.applyEvaluationKeyElisions++ })
	return nil
}

// CheckAndGetGaloisKey is the Fast CKKS boundary for the in-memory key set.
// RLWE's generic getter intentionally rejects Fast keys; this method exposes
// one only to the CKKS zero-secret compatibility methods, which still verify
// its Fast layout and the zero c1/decomposition before eliding key switching.
func (eval Evaluator) CheckAndGetGaloisKey(galEl uint64) (*rlwe.GaloisKey, error) {
	if !eval.usesFastCKKSZeroSecretSimulation() {
		return eval.Evaluator.CheckAndGetGaloisKey(galEl)
	}
	var gk *rlwe.GaloisKey
	switch keySet := eval.Evaluator.EvaluationKeySet.(type) {
	case *rlwe.MemEvaluationKeySet:
		gk = keySet.GaloisKeys[galEl]
	case rlwe.MemEvaluationKeySet:
		gk = keySet.GaloisKeys[galEl]
	default:
		return nil, fmt.Errorf("Fast CKKS zero-secret compatibility requires an in-memory Galois key set")
	}
	if gk == nil {
		return nil, fmt.Errorf("GaloisKey[%d] is missing", galEl)
	}
	if gk.Layout != rlwe.KeyLayoutFast || gk.GadgetCiphertext.Layout != rlwe.KeyLayoutFast {
		return nil, fmt.Errorf("Fast CKKS zero-secret compatibility requires a Fast Galois key")
	}
	degree, _, _, err := fastGadgetCiphertextShape(&gk.GadgetCiphertext, eval.GetParameters())
	if err != nil || degree != 1 {
		return nil, fmt.Errorf("Fast CKKS zero-secret compatibility requires a well-formed degree-one Galois key")
	}
	_ = eval.AutomorphismIndex(galEl)
	return gk, nil
}

// DecomposeNTT has no error return in the RLWE API. For Fast zero-secret
// inputs, a verified zero c1 decomposes exactly to zero, so the expensive
// decomposition is elided. Errors are retained and surfaced by the next
// error-returning boundary (or Trace); output buffers remain untouched.
func (eval Evaluator) DecomposeNTT(levelQ, levelP, nbPi int, c2 ring.Poly, c2IsNTT bool, decompQP []ringqp.Poly) {
	if !eval.usesFastCKKSZeroSecretSimulation() {
		eval.Evaluator.DecomposeNTT(levelQ, levelP, nbPi, c2, c2IsNTT, decompQP)
		return
	}
	state := eval.fastQPCompat
	if state.peek() != nil {
		return
	}
	if err := eval.validateZeroDecompositionInput(levelQ, levelP, c2, decompQP); err != nil {
		state.record(err)
		return
	}
	zeroQPDecomposition(levelQ, levelP, decompQP)
	state.note(func(c *fastQPCompatibilityCounters) { c.decomposeElisions++ })
}

// GadgetProductLazy returns the exact zero QP product only when the source
// polynomial is proven zero on its authoritative Q-prefix. It does not run a
// Standard gadget product with Fast-layout key material.
func (eval Evaluator) GadgetProductLazy(levelQ int, cx ring.Poly, gadgetCt *rlwe.GadgetCiphertext, ctQP *rlwe.Element[ringqp.Poly]) error {
	if !eval.usesFastCKKSZeroSecretSimulation() {
		return eval.Evaluator.GadgetProductLazy(levelQ, cx, gadgetCt, ctQP)
	}
	if err := eval.fastQPCompat.take(); err != nil {
		return fmt.Errorf("cannot Fast CKKS GadgetProductLazy: previous compatibility operation failed: %w", err)
	}
	keyDegree, keyLevelQ, levelP, err := fastGadgetCiphertextShape(gadgetCt, eval.GetParameters())
	if err != nil || keyDegree == 0 {
		return fmt.Errorf("cannot Fast CKKS GadgetProductLazy: malformed Fast gadget ciphertext")
	}
	if levelQ < 0 || levelQ > eval.GetParameters().MaxLevel() || len(cx.Coeffs) < levelQ+1 || cx.N() != eval.GetParameters().N() {
		return fmt.Errorf("cannot Fast CKKS GadgetProductLazy: input polynomial has invalid Q level or backing")
	}
	if levelQ > keyLevelQ {
		return fmt.Errorf("cannot Fast CKKS GadgetProductLazy: input level %d exceeds Fast gadget level %d", levelQ, keyLevelQ)
	}
	rows, err := fastcore.QPrefixWidth(levelQ)
	if err != nil {
		return err
	}
	if err = fastcore.ValidatePrefixRows(eval.GetParameters().RingQ(), levelQ, rows, cx); err != nil {
		return fmt.Errorf("cannot Fast CKKS GadgetProductLazy: %w", err)
	}
	if !fastQPolyIsZero(cx, rows) {
		return fmt.Errorf("cannot Fast CKKS GadgetProductLazy: source polynomial is nonzero on the authoritative Q-prefix")
	}
	if err = validateFastQPOutput(eval.GetParameters(), ctQP, levelQ, levelP); err != nil {
		return fmt.Errorf("cannot Fast CKKS GadgetProductLazy: %w", err)
	}
	zeroFastQPOutput(ctQP, levelQ, levelP)
	eval.fastQPCompat.note(func(c *fastQPCompatibilityCounters) { c.gadgetProductElisions++ })
	return nil
}

// GadgetProductHoistedLazy skips the product only for a fully shaped, all-zero
// QP decomposition. This is the exact output required by BTS's zero-secret
// key-switch inputs.
func (eval Evaluator) GadgetProductHoistedLazy(levelQ int, decompQP []ringqp.Poly, gadgetCt *rlwe.GadgetCiphertext, ctQP *rlwe.Element[ringqp.Poly]) error {
	if !eval.usesFastCKKSZeroSecretSimulation() {
		return eval.Evaluator.GadgetProductHoistedLazy(levelQ, decompQP, gadgetCt, ctQP)
	}
	if err := eval.fastQPCompat.take(); err != nil {
		return fmt.Errorf("cannot Fast CKKS GadgetProductHoistedLazy: previous compatibility operation failed: %w", err)
	}
	keyDegree, keyLevelQ, levelP, err := fastGadgetCiphertextShape(gadgetCt, eval.GetParameters())
	if err != nil || keyDegree == 0 {
		return fmt.Errorf("cannot Fast CKKS GadgetProductHoistedLazy: malformed Fast gadget ciphertext")
	}
	if levelQ > keyLevelQ {
		return fmt.Errorf("cannot Fast CKKS GadgetProductHoistedLazy: input level %d exceeds Fast gadget level %d", levelQ, keyLevelQ)
	}
	if err := validateFastQPDecomposition(eval.GetParameters(), levelQ, levelP, decompQP, true); err != nil {
		return fmt.Errorf("cannot Fast CKKS GadgetProductHoistedLazy: %w", err)
	}
	if !fastQPDecompositionIsZero(decompQP, levelQ, levelP) {
		return fmt.Errorf("cannot Fast CKKS GadgetProductHoistedLazy: decomposition is nonzero")
	}
	if err := validateFastQPOutput(eval.GetParameters(), ctQP, levelQ, levelP); err != nil {
		return fmt.Errorf("cannot Fast CKKS GadgetProductHoistedLazy: %w", err)
	}
	zeroFastQPOutput(ctQP, levelQ, levelP)
	eval.fastQPCompat.note(func(c *fastQPCompatibilityCounters) { c.hoistedProductElisions++ })
	return nil
}

// GadgetProductHoisted is the legacy void-returning Q-only interface used by
// BTS ModUp. A malformed/nonzero decomposition is recorded and surfaced by
// the following Trace call; valid zero input writes an exact zero ciphertext.
func (eval Evaluator) GadgetProductHoisted(levelQ int, decompQP []ringqp.Poly, gadgetCt *rlwe.GadgetCiphertext, ct *rlwe.Ciphertext) {
	if !eval.usesFastCKKSZeroSecretSimulation() {
		eval.Evaluator.GadgetProductHoisted(levelQ, decompQP, gadgetCt, ct)
		return
	}
	state := eval.fastQPCompat
	if state.peek() != nil {
		return
	}
	keyDegree, keyLevelQ, levelP, shapeErr := fastGadgetCiphertextShape(gadgetCt, eval.GetParameters())
	if shapeErr != nil || keyDegree == 0 {
		state.record(fmt.Errorf("cannot Fast CKKS GadgetProductHoisted: malformed Fast gadget ciphertext"))
		return
	}
	if levelQ > keyLevelQ {
		state.record(fmt.Errorf("cannot Fast CKKS GadgetProductHoisted: input level %d exceeds Fast gadget level %d", levelQ, keyLevelQ))
		return
	}
	if err := validateFastQPDecomposition(eval.GetParameters(), levelQ, levelP, decompQP, true); err != nil {
		state.record(fmt.Errorf("cannot Fast CKKS GadgetProductHoisted: %w", err))
		return
	}
	if !fastQPDecompositionIsZero(decompQP, levelQ, levelP) {
		state.record(fmt.Errorf("cannot Fast CKKS GadgetProductHoisted: decomposition is nonzero"))
		return
	}
	level, rows, err := eval.validateFastQPCompatCiphertext("output", ct, false)
	if err != nil || levelQ > level {
		if err == nil {
			err = fmt.Errorf("output level %d is below requested level %d", level, levelQ)
		}
		state.record(fmt.Errorf("cannot Fast CKKS GadgetProductHoisted: %w", err))
		return
	}
	if err = validateFastPolyRows("output c0", ct.Value[0], levelQ, ct.N(), false); err == nil {
		err = validateFastPolyRows("output c1", ct.Value[1], levelQ, ct.N(), false)
	}
	if err != nil {
		state.record(fmt.Errorf("cannot Fast CKKS GadgetProductHoisted: %w", err))
		return
	}
	writeRows := min(levelQ+1, rows)
	if fullFastPolyRows(ct.Value[0], levelQ, ct.N()) && fullFastPolyRows(ct.Value[1], levelQ, ct.N()) {
		writeRows = levelQ + 1
	}
	for row := 0; row < writeRows; row++ {
		clear(ct.Value[0].Coeffs[row])
		clear(ct.Value[1].Coeffs[row])
	}
	state.note(func(c *fastQPCompatibilityCounters) { c.hoistedProductElisions++ })
}

// AutomorphismHoistedLazy retains the actual Galois automorphism while
// eliding key switching only after both the ciphertext c1 and supplied
// decomposition are verified zero. Its result is P-scaled in Q and zero in P,
// so the existing ModDownQPtoQNTT recovers the automorphed c0 exactly.
func (eval Evaluator) AutomorphismHoistedLazy(levelQ int, ctIn *rlwe.Ciphertext, decompQP []ringqp.Poly, galEl uint64, ctQP *rlwe.Element[ringqp.Poly]) error {
	if !eval.usesFastCKKSZeroSecretSimulation() {
		return eval.Evaluator.AutomorphismHoistedLazy(levelQ, ctIn, decompQP, galEl, ctQP)
	}
	if err := eval.fastQPCompat.take(); err != nil {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: previous compatibility operation failed: %w", err)
	}
	level, _, err := eval.validateFastQPCompatCiphertext("input", ctIn, false)
	if err != nil {
		return err
	}
	if levelQ < 0 || levelQ > level {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: requested level %d exceeds input level %d", levelQ, level)
	}
	if galEl == 1 {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: identity automorphism does not require a hoisted key")
	}
	var gk *rlwe.GaloisKey
	if gk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: %w", err)
	}
	if gk.Layout != rlwe.KeyLayoutFast || gk.GadgetCiphertext.Layout != rlwe.KeyLayoutFast {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: a Fast Galois key is required")
	}
	keyDegree, keyLevelQ, levelP, shapeErr := fastGadgetCiphertextShape(&gk.GadgetCiphertext, eval.GetParameters())
	if shapeErr != nil || keyDegree == 0 {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: malformed Fast Galois key")
	}
	if levelQ > keyLevelQ {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: input level %d exceeds Galois-key level %d", levelQ, keyLevelQ)
	}
	if err = validateFastQPOutput(eval.GetParameters(), ctQP, levelQ, levelP); err != nil {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: %w", err)
	}
	if ctIn.IsNTT != ctQP.IsNTT || ctIn.IsMontgomery != ctQP.IsMontgomery {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: input and QP output domains must match")
	}
	if err = validateFastQPDecomposition(eval.GetParameters(), levelQ, levelP, decompQP, true); err != nil {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: %w", err)
	}
	if !fastQPDecompositionIsZero(decompQP, levelQ, levelP) {
		return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: c1 decomposition is nonzero")
	}
	if levelQ+1 > fastcore.MaxQPrefixWidth {
		if err = eval.verifyFastFullQCiphertext(ctIn); err != nil {
			return fmt.Errorf("cannot Fast CKKS AutomorphismHoistedLazy: full-Q input is not authorized: %w", err)
		}
	}
	params := eval.GetParameters()
	ringQ := params.RingQ().AtLevel(levelQ)
	scaled := ringQ.NewPoly()
	if levelP >= 0 {
		ringQ.MulScalarBigint(ctIn.Value[0], params.RingP().ModulusAtLevel[levelP], scaled)
	} else {
		scaled.CopyLvl(levelQ, ctIn.Value[0])
	}
	if ctIn.IsNTT {
		ringQ.AutomorphismNTTWithIndex(scaled, eval.AutomorphismIndex(galEl), ctQP.Value[0].Q)
	} else {
		ringQ.Automorphism(scaled, galEl, ctQP.Value[0].Q)
	}
	for row := 0; row <= levelP; row++ {
		clear(ctQP.Value[0].P.Coeffs[row])
	}
	for row := 0; row <= levelQ; row++ {
		clear(ctQP.Value[1].Q.Coeffs[row])
	}
	for row := 0; row <= levelP; row++ {
		clear(ctQP.Value[1].P.Coeffs[row])
	}
	*ctQP.MetaData = *ctIn.MetaData
	eval.fastQPCompat.note(func(c *fastQPCompatibilityCounters) {
		c.keySwitchElisions++
		c.automorphisms++
	})
	return nil
}

// ModDownQPtoQNTT validates the complete QP representation before invoking
// the shared basis extender. The underlying API cannot return an error, so a
// malformed call is deferred to the next error-returning evaluator boundary.
func (eval Evaluator) ModDownQPtoQNTT(levelQ, levelP int, p1Q, p1P, p2Q ring.Poly) {
	if !eval.usesFastCKKSZeroSecretSimulation() {
		eval.Evaluator.ModDownQPtoQNTT(levelQ, levelP, p1Q, p1P, p2Q)
		return
	}
	if eval.fastQPCompat.peek() != nil {
		return
	}
	if levelQ < 0 || levelQ > eval.GetParameters().MaxLevel() || levelP < 0 || levelP > eval.GetParameters().MaxLevelP() {
		eval.fastQPCompat.record(fmt.Errorf("cannot Fast CKKS ModDownQPtoQNTT: invalid Q/P level (%d,%d)", levelQ, levelP))
		return
	}
	if levelQ+1 > fastcore.MaxQPrefixWidth && !eval.fastQPCompat.isArmed() {
		eval.fastQPCompat.record(fmt.Errorf("cannot Fast CKKS ModDownQPtoQNTT: high-level dormant Q rows are not authorized"))
		return
	}
	if err := validateFastPolyRows("ModDown Q input", p1Q, levelQ, eval.GetParameters().N(), true); err != nil {
		eval.fastQPCompat.record(err)
		return
	}
	if err := validateFastPolyRows("ModDown P input", p1P, levelP, eval.GetParameters().N(), true); err != nil {
		eval.fastQPCompat.record(err)
		return
	}
	if !fastQPolyIsZero(p1P, levelP+1) {
		eval.fastQPCompat.record(fmt.Errorf("cannot Fast CKKS ModDownQPtoQNTT: P input is nonzero in the zero-secret lazy-QP lane"))
		return
	}
	if err := eval.verifyFastPScaledFullQPoly(p1Q, levelQ, levelP); err != nil {
		eval.fastQPCompat.record(fmt.Errorf("cannot Fast CKKS ModDownQPtoQNTT: P-scaled Q input is not authorized by the Q-prefix representation: %w", err))
		return
	}
	if err := validateFastPolyRows("ModDown Q output", p2Q, levelQ, eval.GetParameters().N(), true); err != nil {
		eval.fastQPCompat.record(err)
		return
	}
	eval.Evaluator.ModDownQPtoQNTT(levelQ, levelP, p1Q, p1P, p2Q)
}

// Trace surfaces errors recorded by the legacy void-returning low-level
// decomposition/key-switch adapters before executing further work.
func (eval Evaluator) Trace(ctIn *rlwe.Ciphertext, logN int, opOut *rlwe.Ciphertext) error {
	if eval.usesFastCKKSZeroSecretSimulation() {
		if err := eval.fastQPCompat.take(); err != nil {
			return fmt.Errorf("cannot Fast CKKS Trace: previous compatibility operation failed: %w", err)
		}
	}
	return eval.Evaluator.Trace(ctIn, logN, opOut)
}

func (eval Evaluator) validateFastQPCompatCiphertext(name string, ct *rlwe.Ciphertext, requireFull bool) (level, rows int, err error) {
	if ct == nil || ct.MetaData == nil || len(ct.Value) != 2 {
		return 0, 0, fmt.Errorf("cannot Fast CKKS QP compatibility: %s must be a degree-one ciphertext with metadata", name)
	}
	if len(ct.Value[0].Coeffs) == 0 || len(ct.Value[0].Coeffs[0]) == 0 {
		return 0, 0, fmt.Errorf("cannot Fast CKKS QP compatibility: %s has no c0 backing", name)
	}
	level = len(ct.Value[0].Coeffs) - 1
	if level < 0 || level > eval.GetParameters().MaxLevel() || ct.N() != eval.GetParameters().N() {
		return 0, 0, fmt.Errorf("cannot Fast CKKS QP compatibility: %s has an invalid level or ring degree", name)
	}
	for component := range ct.Value {
		if len(ct.Value[component].Coeffs) != level+1 {
			return 0, 0, fmt.Errorf("cannot Fast CKKS QP compatibility: %s component %d has inconsistent logical Q rows", name, component)
		}
	}
	rows, err = fastcore.QPrefixWidth(level)
	if err != nil {
		return 0, 0, err
	}
	if err = fastcore.ValidatePrefixRows(eval.GetParameters().RingQ(), level, rows, ct.Value...); err != nil {
		return 0, 0, fmt.Errorf("cannot Fast CKKS QP compatibility: %s: %w", name, err)
	}
	if !fastQPolyIsZero(ct.Value[1], rows) {
		return 0, 0, fmt.Errorf("cannot Fast CKKS QP compatibility: %s c1 is nonzero on authoritative q-prefix rows", name)
	}
	if requireFull {
		for component := range ct.Value {
			if err = validateFastPolyRows(fmt.Sprintf("%s component %d", name, component), ct.Value[component], level, ct.N(), true); err != nil {
				return 0, 0, err
			}
		}
		if !fastQPolyIsZero(ct.Value[1], level+1) {
			return 0, 0, fmt.Errorf("cannot Fast CKKS QP compatibility: %s c1 is nonzero in a full-Q row", name)
		}
	}
	return level, rows, nil
}

func (eval Evaluator) validateZeroDecompositionInput(levelQ, levelP int, c2 ring.Poly, decompQP []ringqp.Poly) error {
	if levelQ < 0 || levelQ > eval.GetParameters().MaxLevel() || levelP < -1 || levelP > eval.GetParameters().MaxLevelP() {
		return fmt.Errorf("cannot Fast CKKS DecomposeNTT: invalid Q/P level (%d,%d)", levelQ, levelP)
	}
	if len(c2.Coeffs) < levelQ+1 || len(c2.Coeffs) == 0 || len(c2.Coeffs[0]) != eval.GetParameters().N() {
		return fmt.Errorf("cannot Fast CKKS DecomposeNTT: input polynomial has insufficient dimensions or Q rows")
	}
	rows, err := fastcore.QPrefixWidth(levelQ)
	if err != nil {
		return err
	}
	if err = fastcore.ValidatePrefixRows(eval.GetParameters().RingQ(), levelQ, rows, c2); err != nil {
		return fmt.Errorf("cannot Fast CKKS DecomposeNTT: %w", err)
	}
	if !fastQPolyIsZero(c2, rows) {
		return fmt.Errorf("cannot Fast CKKS DecomposeNTT: c1 is nonzero on authoritative q-prefix rows")
	}
	return validateFastQPDecomposition(eval.GetParameters(), levelQ, levelP, decompQP, false)
}

func (eval Evaluator) verifyFastFullQCiphertext(ct *rlwe.Ciphertext) error {
	// Only c0 enters the complete-Q lazy transform. c1 is intentionally a
	// compact zero polynomial, so validate its authoritative prefix without
	// reading or requiring dormant rows.
	level, _, err := eval.validateFastQPCompatCiphertext("full-Q input", ct, false)
	if err != nil {
		return err
	}
	if err = validateFastPolyRows("full-Q c0", ct.Value[0], level, ct.N(), true); err != nil {
		return err
	}
	fingerprint := fastQFingerprint(ct, level)
	if eval.fastQPCompat.trustedFingerprint(ct, fingerprint) {
		return nil
	}
	if err = eval.verifyFastFullQPoly(ct.Value[0], level, ct.IsNTT, ct.IsMontgomery); err != nil {
		return err
	}
	eval.fastQPCompat.trust(ct, fingerprint)
	return nil
}

func (eval Evaluator) verifyFastFullQPoly(poly ring.Poly, level int, isNTT, isMontgomery bool) error {
	if err := validateFastPolyRows("full-Q polynomial", poly, level, eval.GetParameters().N(), true); err != nil {
		return err
	}
	ringQ := eval.GetParameters().RingQ()
	fullRing := ringQ.AtLevel(level)
	normal := ring.NewPoly(eval.GetParameters().N(), level)
	if isMontgomery {
		fullRing.IMForm(poly, normal)
	} else {
		normal.CopyLvl(level, poly)
	}
	coefficients := ring.NewPoly(eval.GetParameters().N(), level)
	if isNTT {
		fullRing.INTT(normal, coefficients)
	} else {
		coefficients.CopyLvl(level, normal)
	}
	return eval.verifyFastFullQCoefficients(coefficients, level)
}

// verifyFastPScaledFullQPoly checks a BTS lazy-QP Q component without
// incorrectly requiring P*x itself to fit in the Q-prefix. It removes the
// known P factor modulo each q_i, then proves that x's dormant rows match its
// centered Q-prefix representative (2B < S_Q). No dormant row is reconstructed.
func (eval Evaluator) verifyFastPScaledFullQPoly(poly ring.Poly, levelQ, levelP int) error {
	params := eval.GetParameters()
	if levelP < 0 || levelP > params.MaxLevelP() {
		return fmt.Errorf("invalid P level %d", levelP)
	}
	if err := validateFastPolyRows("P-scaled full-Q polynomial", poly, levelQ, params.N(), true); err != nil {
		return err
	}
	ringQ := params.RingQ()
	coefficients := ring.NewPoly(params.N(), levelQ)
	ringQ.AtLevel(levelQ).INTT(poly, coefficients)

	pFactor := params.RingP().ModulusAtLevel[levelP]
	qProduct := ringQ.ModulusAtLevel[levelQ]
	pInverse := new(big.Int).ModInverse(pFactor, qProduct)
	if pInverse == nil {
		return fmt.Errorf("P product at level %d is not invertible modulo the active Q product", levelP)
	}
	recovered := ring.NewPoly(params.N(), levelQ)
	ringQ.AtLevel(levelQ).MulScalarBigint(coefficients, pInverse, recovered)
	return eval.verifyFastFullQCoefficients(recovered, levelQ)
}

func (eval Evaluator) verifyFastFullQCoefficients(coefficients ring.Poly, level int) error {
	rows, err := fastcore.QPrefixWidth(level)
	if err != nil {
		return err
	}
	if err = validateFastPolyRows("full-Q coefficient polynomial", coefficients, level, eval.GetParameters().N(), true); err != nil {
		return err
	}
	ringQ := eval.GetParameters().RingQ()
	prefixRing := ringQ.AtLevel(rows - 1)
	for row := 0; row <= level; row++ {
		modulus := ringQ.SubRings[row].Modulus
		for coefficient, value := range coefficients.Coeffs[row] {
			if value >= modulus {
				return fmt.Errorf("full-Q polynomial q%d coefficient %d is not a canonical residue", row, coefficient)
			}
		}
	}
	if level+1 <= rows {
		return nil
	}
	centered := make([]*big.Int, eval.GetParameters().N())
	for i := range centered {
		centered[i] = new(big.Int)
	}
	prefixRing.PolyToBigintCentered(coefficients, 1, centered)
	capacity := new(big.Int).Set(ringQ.ModulusAtLevel[rows-1])
	capacityHalf := new(big.Int).Rsh(new(big.Int).Set(capacity), 1)
	moduli := make([]*big.Int, level+1)
	for row := rows; row <= level; row++ {
		moduli[row] = new(big.Int).SetUint64(ringQ.SubRings[row].Modulus)
	}
	tmp, abs := new(big.Int), new(big.Int)
	for i, value := range centered {
		if abs.Abs(value).Cmp(capacityHalf) > 0 {
			return fmt.Errorf("Q-prefix CRT representative exceeds S_Q/2 at coefficient %d", i)
		}
		for row := rows; row <= level; row++ {
			tmp.Mod(value, moduli[row])
			if coefficients.Coeffs[row][i] != tmp.Uint64() {
				return fmt.Errorf("dormant q%d row disagrees with the bounded Q-prefix CRT representative at coefficient %d", row, i)
			}
		}
	}
	return nil
}

func (eval Evaluator) fastQPRescale(op0, opOut *rlwe.Ciphertext) error {
	if err := eval.fastQPCompat.take(); err != nil {
		return fmt.Errorf("cannot Fast CKKS Rescale: previous compatibility operation failed: %w", err)
	}
	level, _, err := eval.validateFastQPCompatCiphertext("input", op0, true)
	if err != nil {
		return err
	}
	if err = eval.verifyFastFullQCiphertext(op0); err != nil {
		return fmt.Errorf("cannot Fast CKKS Rescale: %w", err)
	}
	if opOut == nil || opOut.MetaData == nil || len(opOut.Value) != 2 {
		return fmt.Errorf("cannot Fast CKKS Rescale: output must have degree-one storage and metadata")
	}
	if len(opOut.Value[0].Coeffs) == 0 || len(opOut.Value[0].Coeffs[0]) != eval.GetParameters().N() {
		return fmt.Errorf("cannot Fast CKKS Rescale: output ring degree does not match parameters")
	}
	if !op0.IsNTT {
		return fmt.Errorf("cannot Fast CKKS Rescale: full-Q compatibility Rescale requires NTT-domain input")
	}
	if op0.IsMontgomery {
		return fmt.Errorf("cannot Fast CKKS Rescale: full-Q compatibility Rescale does not accept Montgomery-form input")
	}
	nbRescales := eval.GetParameters().LevelsConsumedPerRescaling()
	if level < nbRescales {
		return fmt.Errorf("cannot Fast CKKS Rescale: input Ciphertext level is too low")
	}
	targetLevel := level - nbRescales
	stage := rlwe.NewCiphertext(eval.GetParameters(), 1, targetLevel)
	*stage.MetaData = *op0.MetaData
	ringQ := eval.GetParameters().RingQ().AtLevel(level)
	for i := 0; i < nbRescales; i++ {
		stage.Scale = stage.Scale.Div(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
	}
	buffer := ringQ.NewPoly()
	for component := range stage.Value {
		ringQ.DivRoundByLastModulusManyNTT(nbRescales, op0.Value[component], buffer, stage.Value[component])
	}
	if err = eval.verifyFastFullQPoly(stage.Value[0], targetLevel, stage.IsNTT, stage.IsMontgomery); err != nil {
		return fmt.Errorf("cannot Fast CKKS Rescale: output Q rows failed the bounded-prefix check: %w", err)
	}
	if !fastQPolyIsZero(stage.Value[1], targetLevel+1) {
		return fmt.Errorf("cannot Fast CKKS Rescale: output c1 is nonzero")
	}
	opOut.Resize(1, targetLevel)
	for component := range stage.Value {
		for row := 0; row <= targetLevel; row++ {
			if len(opOut.Value[component].Coeffs[row]) != eval.GetParameters().N() {
				opOut.Value[component].Coeffs[row] = make([]uint64, eval.GetParameters().N())
			}
			copy(opOut.Value[component].Coeffs[row], stage.Value[component].Coeffs[row])
		}
	}
	*opOut.MetaData = *stage.MetaData
	fingerprint := fastQFingerprint(opOut, targetLevel)
	eval.fastQPCompat.trust(opOut, fingerprint)
	eval.fastQPCompat.note(func(c *fastQPCompatibilityCounters) { c.fullQRescales++ })
	return nil
}

// fastGadgetCiphertextShape reads only the first QP vector needed by the
// GadgetCiphertext metadata accessors, after validating every backing slice
// those accessors would otherwise index. The zero-product shortcuts do not
// consume key coefficients, but still require a well-formed Fast key header.
func fastGadgetCiphertextShape(gadgetCt *rlwe.GadgetCiphertext, params *Parameters) (degree, levelQ, levelP int, err error) {
	if gadgetCt == nil || params == nil || gadgetCt.Layout != rlwe.KeyLayoutFast || len(gadgetCt.Value) == 0 || len(gadgetCt.Value[0]) == 0 {
		return 0, 0, 0, fmt.Errorf("missing Fast gadget ciphertext storage")
	}
	vector := gadgetCt.Value[0][0]
	if len(vector) != 2 {
		return 0, 0, 0, fmt.Errorf("Fast gadget ciphertext must contain two polynomial components")
	}
	levelQ = len(vector[0].Q.Coeffs) - 1
	if levelQ < 0 || levelQ > params.MaxLevel() {
		return 0, 0, 0, fmt.Errorf("Fast gadget ciphertext has invalid Q rows")
	}
	if err = validateFastPolyRows("Fast gadget Q component 0", vector[0].Q, levelQ, params.N(), true); err != nil {
		return 0, 0, 0, err
	}
	if err = validateFastPolyRows("Fast gadget Q component 1", vector[1].Q, levelQ, params.N(), true); err != nil {
		return 0, 0, 0, err
	}
	levelP = len(vector[0].P.Coeffs) - 1
	if len(vector[0].P.Coeffs) == 0 {
		levelP = -1
	} else {
		if levelP > params.MaxLevelP() {
			return 0, 0, 0, fmt.Errorf("Fast gadget ciphertext has invalid P rows")
		}
		if err = validateFastPolyRows("Fast gadget P component 0", vector[0].P, levelP, params.N(), true); err != nil {
			return 0, 0, 0, err
		}
		if err = validateFastPolyRows("Fast gadget P component 1", vector[1].P, levelP, params.N(), true); err != nil {
			return 0, 0, 0, err
		}
	}
	return len(vector) - 1, levelQ, levelP, nil
}

func (eval Evaluator) validateFastQPDecomposition(levelQ, levelP int, decompQP []ringqp.Poly, requireP bool) error {
	return validateFastQPDecomposition(eval.GetParameters(), levelQ, levelP, decompQP, requireP)
}

func validateFastQPDecomposition(params *Parameters, levelQ, levelP int, decompQP []ringqp.Poly, requireP bool) error {
	if params == nil || levelQ < 0 || levelQ > params.MaxLevel() || levelP < -1 || levelP > params.MaxLevelP() {
		return fmt.Errorf("invalid decomposition Q/P level (%d,%d)", levelQ, levelP)
	}
	beta := params.BaseRNSDecompositionVectorSize(levelQ, levelP)
	if len(decompQP) != beta {
		return fmt.Errorf("decomposition has %d digits, expected %d", len(decompQP), beta)
	}
	for i := range decompQP {
		if err := validateFastPolyRows(fmt.Sprintf("decomposition digit %d Q", i), decompQP[i].Q, levelQ, params.N(), true); err != nil {
			return err
		}
		if levelP >= 0 {
			if err := validateFastPolyRows(fmt.Sprintf("decomposition digit %d P", i), decompQP[i].P, levelP, params.N(), true); err != nil {
				return err
			}
		} else if requireP && len(decompQP[i].P.Coeffs) != 0 {
			return fmt.Errorf("decomposition digit %d unexpectedly contains P rows", i)
		}
	}
	return nil
}

func validateFastQPOutput(params *Parameters, output *rlwe.Element[ringqp.Poly], levelQ, levelP int) error {
	if output == nil || output.MetaData == nil || len(output.Value) != 2 {
		return fmt.Errorf("QP output must contain two polynomials and metadata")
	}
	for component := range output.Value {
		if err := validateFastPolyRows(fmt.Sprintf("QP output component %d Q", component), output.Value[component].Q, levelQ, params.N(), true); err != nil {
			return err
		}
		if levelP >= 0 {
			if err := validateFastPolyRows(fmt.Sprintf("QP output component %d P", component), output.Value[component].P, levelP, params.N(), true); err != nil {
				return err
			}
		}
	}
	return nil
}

func validateFastPolyRows(name string, poly ring.Poly, level, n int, full bool) error {
	if level < 0 || len(poly.Coeffs) < level+1 || poly.N() != n {
		return fmt.Errorf("%s has insufficient logical rows or wrong ring degree", name)
	}
	rows := min(level+1, fastcore.MaxQPrefixWidth)
	if full {
		rows = level + 1
	}
	for row := 0; row < rows; row++ {
		if len(poly.Coeffs[row]) != n {
			return fmt.Errorf("%s q/p row %d has backing length %d, expected %d", name, row, len(poly.Coeffs[row]), n)
		}
	}
	return nil
}

func fullFastPolyRows(poly ring.Poly, level, n int) bool {
	if len(poly.Coeffs) < level+1 || poly.N() != n {
		return false
	}
	for row := 0; row <= level; row++ {
		if len(poly.Coeffs[row]) != n {
			return false
		}
	}
	return true
}

func fastQPolyIsZero(poly ring.Poly, rows int) bool {
	if len(poly.Coeffs) < rows {
		return false
	}
	for row := 0; row < rows; row++ {
		for _, value := range poly.Coeffs[row] {
			if value != 0 {
				return false
			}
		}
	}
	return true
}

func fastQPDecompositionIsZero(decompQP []ringqp.Poly, levelQ, levelP int) bool {
	for i := range decompQP {
		if !fastQPolyIsZero(decompQP[i].Q, levelQ+1) {
			return false
		}
		if levelP >= 0 && !fastQPolyIsZero(decompQP[i].P, levelP+1) {
			return false
		}
	}
	return true
}

func zeroQPDecomposition(levelQ, levelP int, decompQP []ringqp.Poly) {
	for i := range decompQP {
		for row := 0; row <= levelQ; row++ {
			clear(decompQP[i].Q.Coeffs[row])
		}
		if levelP >= 0 {
			for row := 0; row <= levelP; row++ {
				clear(decompQP[i].P.Coeffs[row])
			}
		}
	}
}

func zeroFastQPOutput(output *rlwe.Element[ringqp.Poly], levelQ, levelP int) {
	for i := range output.Value {
		for row := 0; row <= levelQ; row++ {
			clear(output.Value[i].Q.Coeffs[row])
		}
		if levelP >= 0 {
			for row := 0; row <= levelP; row++ {
				clear(output.Value[i].P.Coeffs[row])
			}
		}
	}
}

func fastQFingerprint(ct *rlwe.Ciphertext, level int) fastQPFingerprint {
	hash := sha256.New()
	var encoded [8]byte
	flags := uint64(level) << 2
	if ct.IsNTT {
		flags |= 1
	}
	if ct.IsMontgomery {
		flags |= 2
	}
	binary.LittleEndian.PutUint64(encoded[:], flags)
	_, _ = hash.Write(encoded[:])
	for row := 0; row <= level; row++ {
		coefficients := ct.Value[0].Coeffs[row]
		const wordsPerChunk = 1024
		chunk := make([]byte, wordsPerChunk*8)
		for start := 0; start < len(coefficients); start += wordsPerChunk {
			end := min(start+wordsPerChunk, len(coefficients))
			for i := start; i < end; i++ {
				binary.LittleEndian.PutUint64(chunk[(i-start)*8:], coefficients[i])
			}
			_, _ = hash.Write(chunk[:(end-start)*8])
		}
	}
	var sum [sha256.Size]byte
	copy(sum[:], hash.Sum(nil))
	return fastQPFingerprint{level: level, hash: sum}
}
