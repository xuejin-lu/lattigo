//go:build fastdiag

package fastdiag

import (
	"testing"

	"github.com/stretchr/testify/require"
)

func TestScopeParserAndUnknownScope(t *testing.T) {
	require.NoError(t, ConfigureCSV("stage,power,rescale"))
	require.NoError(t, ConfigureCSV("all"))
	require.ErrorContains(t, ConfigureCSV("stage,typo"), "unknown fastdiag scope typo")
}

func TestNestedEventsHaveStableParentAndSequence(t *testing.T) {
	require.NoError(t, ConfigureCSV("stage,power"))
	Reset()
	parent := Begin(Stage, "bootstrap", 0, Fields{})
	child := Begin(Power, "generated", parent.Sequence(), Fields{})
	child.End(Fields{})
	parent.End(Fields{})
	events := Events()
	require.Len(t, events, 2)
	require.Equal(t, uint64(1), events[0].Sequence)
	require.Equal(t, uint64(2), events[1].Sequence)
	require.Equal(t, parent.Sequence(), events[1].ParentSequence)
}

func TestUnselectedScopeEmitsNoEvent(t *testing.T) {
	require.NoError(t, ConfigureCSV("stage"))
	Reset()
	require.True(t, Selected(Stage))
	require.False(t, Selected(Power))
	require.False(t, Selected(Rescale))
	span := Begin(Power, "ignored", 0, Fields{})
	span.End(Fields{})
	require.Empty(t, Events())
}
