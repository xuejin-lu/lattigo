# Current Task

This Secondary repository does not independently define the active project task.

For orchestrated `fast-qprefix` work, the authoritative task is in the Primary repository `xuejin-lu/heart-lattigo-bootstrap`:

- `CURRENT_TASK.md`
- the active specification referenced by that file.

This file must not introduce task scope, architecture, semantics, or an independent `READY` / `COMPLETE` status.

Safely synchronize and read the Primary active task before acting. If Primary is blocked, completed, or does not authorize Secondary work, do not start a Secondary task.
