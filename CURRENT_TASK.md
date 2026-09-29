# Current Task

Task: QPREFIX-PERF-MEASURE-LOGN16-001
Status: MEASUREMENT_ONLY

Working branch remains: `fast-qprefix`
Required marker: `FAST-QPREFIX-ACTIVE-001`

The prior numerical diagnosis is temporarily suspended.

Authoritative Primary spec:
`specs/QPREFIX-PERF-MEASURE-LOGN16-001-STANDARD-VS-FAST.md`

Purpose:
Use temporary detached worktrees to measure genuine Standard `5dbffbde...` against current Q-prefix `82601ea...` at LogN16.

Do not:
- modify Secondary production code;
- checkout Standard main in this authoritative worktree;
- commit Secondary;
- resume numerical diagnosis during this task.

Use the Primary LogN16 config, warmup 1, repetitions 7, staged timing.
