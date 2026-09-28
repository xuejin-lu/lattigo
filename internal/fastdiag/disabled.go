//go:build !fastdiag

package fastdiag

const Enabled = false

func Selected(Scope) bool                         { return false }
func ConfigureCSV(string) error                   { return nil }
func Reset()                                      {}
func Begin(Scope, string, uint64, Fields) Span    { return Span{} }
func (Span) Sequence() uint64                     { return 0 }
func (Span) End(Fields)                           {}
func Events() []Event                             { return nil }
func StartTimer() Timer                           { return Timer{} }
func (Timer) ElapsedNS() int64                    { return 0 }
func Record(Scope, string, uint64, Fields, int64) {}
