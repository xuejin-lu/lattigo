package fastdiag

import (
	"math"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

type Scope string

const (
	Stage   Scope = "stage"
	Power   Scope = "power"
	Rescale Scope = "rescale"
)

type Fields struct {
	LevelIn      *int     `json:"level_in,omitempty"`
	LevelOut     *int     `json:"level_out,omitempty"`
	RowsIn       *int     `json:"rows_in,omitempty"`
	RowsOut      *int     `json:"rows_out,omitempty"`
	DegreeIn     *int     `json:"degree_in,omitempty"`
	DegreeOut    *int     `json:"degree_out,omitempty"`
	ScaleLog2In  *float64 `json:"scale_log2_in,omitempty"`
	ScaleLog2Out *float64 `json:"scale_log2_out,omitempty"`
	InPlace      *bool    `json:"in_place,omitempty"`
	Power        *int     `json:"power,omitempty"`
	SplitA       *int     `json:"split_a,omitempty"`
	SplitB       *int     `json:"split_b,omitempty"`
	Component    string   `json:"component,omitempty"`
	Count        *int     `json:"count,omitempty"`
}

type Event struct {
	Scope          Scope  `json:"scope"`
	Name           string `json:"name"`
	ElapsedNS      int64  `json:"elapsed_ns"`
	Sequence       uint64 `json:"sequence"`
	ParentSequence uint64 `json:"parent_sequence,omitempty"`
	Fields
}

type Span struct {
	sequence uint64
	parent   uint64
	scope    Scope
	name     string
	fields   Fields
	started  time.Time
	active   bool
}

type Timer struct{ started time.Time }

func Input(ct *rlwe.Ciphertext, rows int) Fields {
	if ct == nil || ct.MetaData == nil {
		return Fields{}
	}
	level, degree := ct.Level(), ct.Degree()
	fields := Fields{LevelIn: &level, DegreeIn: &degree}
	if rows >= 0 {
		fields.RowsIn = &rows
	}
	if logScale, ok := scaleLog2(ct.Scale); ok {
		fields.ScaleLog2In = &logScale
	}
	return fields
}

func Output(ct *rlwe.Ciphertext, rows int) Fields {
	if ct == nil || ct.MetaData == nil {
		return Fields{}
	}
	level, degree := ct.Level(), ct.Degree()
	fields := Fields{LevelOut: &level, DegreeOut: &degree}
	if rows >= 0 {
		fields.RowsOut = &rows
	}
	if logScale, ok := scaleLog2(ct.Scale); ok {
		fields.ScaleLog2Out = &logScale
	}
	return fields
}

func InPlace(value bool) Fields         { return Fields{InPlace: &value} }
func GeneratedPower(value int) Fields   { return Fields{Power: &value} }
func PowerSplit(a, b int) Fields        { return Fields{SplitA: &a, SplitB: &b} }
func ComponentName(value string) Fields { return Fields{Component: value} }
func RepetitionCount(value int) Fields  { return Fields{Count: &value} }

func (fields Fields) Merge(other Fields) Fields {
	if other.LevelIn != nil {
		fields.LevelIn = other.LevelIn
	}
	if other.LevelOut != nil {
		fields.LevelOut = other.LevelOut
	}
	if other.RowsIn != nil {
		fields.RowsIn = other.RowsIn
	}
	if other.RowsOut != nil {
		fields.RowsOut = other.RowsOut
	}
	if other.DegreeIn != nil {
		fields.DegreeIn = other.DegreeIn
	}
	if other.DegreeOut != nil {
		fields.DegreeOut = other.DegreeOut
	}
	if other.ScaleLog2In != nil {
		fields.ScaleLog2In = other.ScaleLog2In
	}
	if other.ScaleLog2Out != nil {
		fields.ScaleLog2Out = other.ScaleLog2Out
	}
	if other.InPlace != nil {
		fields.InPlace = other.InPlace
	}
	if other.Power != nil {
		fields.Power = other.Power
	}
	if other.SplitA != nil {
		fields.SplitA = other.SplitA
	}
	if other.SplitB != nil {
		fields.SplitB = other.SplitB
	}
	if other.Component != "" {
		fields.Component = other.Component
	}
	if other.Count != nil {
		fields.Count = other.Count
	}
	return fields
}

func scaleLog2(scale rlwe.Scale) (float64, bool) {
	if scale.Value.Sign() <= 0 {
		return 0, false
	}
	var mantissa big.Float
	exponent := scale.Value.MantExp(&mantissa)
	value, _ := mantissa.Float64()
	return math.Log2(value) + float64(exponent), true
}
