# Current Task

Status: WAITING_FOR_PRIMARY_TASK

FAST-STORAGE-008 is accepted at:
`1aafc442595da9af41fdc97b03b91ccb2431fe3a`

Classification:
`FAST_STORAGE_008_WIDTH2_VALID_NOT_COMPETITIVE`

Accepted evidence:
- explicit 3->2 contraction is correct under strict capacity proof;
- actual LogN13 C2S fits width 2 for all four factors;
- F2/F3 authoritative lifts and logical maintained rows match exactly;
- F2 is about 30% faster than F3;
- F2 full four-factor chain is still about 1.93x slower than existing Logical Fast.

Production Bootstrap/C2S remains unchanged.

A final bounded C2S performance-root-cause experiment is appropriate before closing this direction:
- scratch/workspace parity for private-F LinearTransform;
- per-stage timing breakdown for LinearTransform, Rescale, and restore;
- no production routing changes.

Wait for the next Primary task.
