# Current Task

Status: WAITING_FOR_PRIMARY_TASK

FAST-INTEGRATION-001 is accepted at:
`31efadc693559217b48d3e76a2e3655b9e6cd14d`

Accepted production integration:
- `FastEvaluator.modUpBasis` now uses the private-F width-3 bridge;
- compact logical export materializes only maintained q rows;
- canonical q0 midpoint semantics are authoritative;
- downstream Trace/DFT/EvalMod/packing/production Rescale remain on the existing path.

Performance follow-up is required:
- current LogN13 bridge is about 1.53x slower than the Standard basis-raise benchmark;
- allocations are about 87,263/op versus 46/op;
- do not treat the current bridge as the final optimized hot path.

Do not begin the next optimization or downstream private-F migration without a new Primary specification.

Wait for the next Primary task.
