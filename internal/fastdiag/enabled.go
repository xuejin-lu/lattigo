//go:build fastdiag

package fastdiag

import (
	"errors"
	"sort"
	"strings"
	"sync"
	"sync/atomic"
	"time"
)

const Enabled = true

func StartTimer() Timer              { return Timer{started: time.Now()} }
func (timer Timer) ElapsedNS() int64 { return time.Since(timer.started).Nanoseconds() }

var collector struct {
	mu     sync.Mutex
	mask   atomic.Uint32
	next   atomic.Uint64
	events []Event
}

func ConfigureCSV(value string) error {
	var mask uint32
	for _, raw := range strings.Split(value, ",") {
		name := strings.TrimSpace(raw)
		if name == "" {
			continue
		}
		if name == "all" {
			mask = 7
			continue
		}
		switch Scope(name) {
		case Stage:
			mask |= 1
		case Power:
			mask |= 2
		case Rescale:
			mask |= 4
		default:
			return errors.New("unknown fastdiag scope " + name)
		}
	}
	collector.mask.Store(mask)
	return nil
}

func Selected(scope Scope) bool {
	var bit uint32
	switch scope {
	case Stage:
		bit = 1
	case Power:
		bit = 2
	case Rescale:
		bit = 4
	}
	return bit != 0 && collector.mask.Load()&bit != 0
}

func Reset() {
	collector.mu.Lock()
	collector.events = nil
	collector.next.Store(0)
	collector.mu.Unlock()
}

func Begin(scope Scope, name string, parent uint64, fields Fields) Span {
	if !Selected(scope) {
		return Span{}
	}
	return Span{sequence: collector.next.Add(1), parent: parent, scope: scope, name: name, fields: fields, started: time.Now(), active: true}
}

func (span Span) Sequence() uint64 { return span.sequence }

func (span Span) End(fields Fields) {
	if !span.active {
		return
	}
	event := Event{
		Scope: span.scope, Name: span.name, ElapsedNS: time.Since(span.started).Nanoseconds(),
		Sequence: span.sequence, ParentSequence: span.parent, Fields: span.fields.Merge(fields),
	}
	collector.mu.Lock()
	collector.events = append(collector.events, event)
	collector.mu.Unlock()
}

func Record(scope Scope, name string, parent uint64, fields Fields, elapsedNS int64) {
	if !Selected(scope) {
		return
	}
	event := Event{Scope: scope, Name: name, ElapsedNS: elapsedNS, Sequence: collector.next.Add(1), ParentSequence: parent, Fields: fields}
	collector.mu.Lock()
	collector.events = append(collector.events, event)
	collector.mu.Unlock()
}

func Events() []Event {
	collector.mu.Lock()
	defer collector.mu.Unlock()
	out := append([]Event(nil), collector.events...)
	sort.Slice(out, func(i, j int) bool { return out[i].Sequence < out[j].Sequence })
	return out
}
