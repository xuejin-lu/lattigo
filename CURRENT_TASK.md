# Current Task

Task: FAST-INTEGRATION-002
Status: NEEDS_BENCHMARK_EVIDENCE

Candidate implementation:
`57ffb88744c82778c0a9392ecab394e19f712a3d`

Web source review found no blocking correctness defect.

Before acceptance, run and report the required LogN13 performance evidence on this candidate:
- `BenchmarkFastModUpBasisLogN13`
- `BenchmarkStandardModUpBasisLogN13`

Run each at least three times and report:
- ns/op
- B/op
- allocs/op

Required gates:
- Fast allocs/op <= 873;
- Fast B/op materially below ~3.35 MB/op;
- Fast ns/op below ~4.50 ms/op.

Do not change code unless a gate fails and one bounded repair pass is needed under the existing FAST-INTEGRATION-002 spec.

Then report `READY_FOR_WEB_REVIEW` if all gates pass, otherwise `NEEDS_WEB_REVIEW`.
