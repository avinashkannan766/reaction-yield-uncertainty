# Profile: uspto_yields

Raw rows **50000** -> clean rows **46929**

Cleaning steps (what was removed and why):

- missing yield/components/product: removed 0, remaining 50000
- yield out of range: removed 96, remaining 49904 (scale=fraction; 56 within ±5.0 clipped to [0.0,100.0])
- unparsable SMILES: removed 0, remaining 49904 (11 distinct molecules kept via unsanitized fallback (metal complexes etc.))
- element balance: removed 0, remaining 49904 (0 reactions have product elements absent from inputs (flag))
- duplicates: removed 2975, remaining 46929 (merged by mean; 3 groups with replicate std > 20.0 (kept, see yield_std))

## Yield

n = 46929; mean 64.0, std 26.0, median 68.0
- exactly 0: 0.0%   |   >= 90: 19.8%   |   >= 99: 5.6%
- deciles: 0, 25, 39, 50, 60, 68, 76, 83, 89, 96, 100
- (a piled-up distribution means a constant predictor is already hard to beat, and the benchmark's `mean baseline` interval width will be small)

## Components

components per reaction: min 2, median 5, max 18
- unique agent: 2598
- unique reactant: 48875
- unique products: 46309 (reactions per product: median 1)

Role signatures (most common):

- 6528 reactions: {'agent': 3, 'reactant': 2}
- 5457 reactions: {'agent': 4, 'reactant': 2}
- 5365 reactions: {'agent': 2, 'reactant': 2}
- 4452 reactions: {'agent': 1, 'reactant': 2}
- 3667 reactions: {'agent': 5, 'reactant': 2}
- 3181 reactions: {'reactant': 2}

## Group structure (can OOD splits be built?)

| group_by | n groups | largest group | groups with >=20 rows | median group size | group split | group-val split |
|---|---|---|---|---|---|---|
| cluster:100 | 100 | 3.1% | 93 | 374 | yes | yes |
| cluster:300 | 300 | 1.6% | 248 | 80 | yes | yes |

Read: a group split needs a few dozen groups, none dominating (largest group well under 20 % of rows); `group-val` needs more groups still.
