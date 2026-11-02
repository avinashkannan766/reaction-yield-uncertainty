# Profile: uspto_yields

Raw rows **60000** -> clean rows **55644**

Cleaning steps (what was removed and why):

- missing yield/components/product: removed 0, remaining 60000
- yield out of range: removed 118, remaining 59882 (scale=fraction; 73 within ±5.0 clipped to [0.0,100.0])
- unparsable SMILES: removed 0, remaining 59882 (11 distinct molecules kept via unsanitized fallback (metal complexes etc.))
- element balance: removed 0, remaining 59882 (0 reactions have product elements absent from inputs (flag))
- duplicates: removed 4238, remaining 55644 (merged by mean; 5 groups with replicate std > 20.0 (kept, see yield_std))

## Yield

n = 55644; mean 64.0, std 26.0, median 68.0
- exactly 0: 0.0%   |   >= 90: 20.0%   |   >= 99: 5.6%
- deciles: 0, 25, 39, 50, 60, 68, 76, 83, 90, 96, 100
- (a piled-up distribution means a constant predictor is already hard to beat, and the benchmark's `mean baseline` interval width will be small)

## Components

components per reaction: min 2, median 5, max 18
- unique agent: 2907
- unique reactant: 56465
- unique products: 54795 (reactions per product: median 1)

Role signatures (most common):

- 7702 reactions: {'agent': 3, 'reactant': 2}
- 6464 reactions: {'agent': 4, 'reactant': 2}
- 6394 reactions: {'agent': 2, 'reactant': 2}
- 5284 reactions: {'agent': 1, 'reactant': 2}
- 4309 reactions: {'agent': 5, 'reactant': 2}
- 3782 reactions: {'reactant': 2}

## Group structure (can OOD splits be built?)

| group_by | n groups | largest group | groups with >=20 rows | median group size | group split | group-val split |
|---|---|---|---|---|---|---|
| generic_scaffold | 13825 | 12.3% | 220 | 1 | yes | yes |
| scaffold | 26939 | 8.5% | 187 | 1 | yes | yes |
| role:agent | 18865 | 7.6% | 302 | 1 | yes | yes |

Read: a group split needs a few dozen groups, none dominating (largest group well under 20 % of rows); `group-val` needs more groups still.
