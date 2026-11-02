# Profile: uspto_yields

Raw rows **20000** -> clean rows **19401**

Cleaning steps (what was removed and why):

- missing yield/components/product: removed 0, remaining 20000
- yield out of range: removed 43, remaining 19957 (scale=fraction; 23 within ±5.0 clipped to [0.0,100.0])
- unparsable SMILES: removed 0, remaining 19957 (8 distinct molecules kept via unsanitized fallback (metal complexes etc.))
- element balance: removed 0, remaining 19957 (0 reactions have product elements absent from inputs (flag))
- duplicates: removed 556, remaining 19401 (merged by mean; 1 groups with replicate std > 20.0 (kept, see yield_std))

## Yield

n = 19401; mean 64.0, std 25.9, median 68.0
- exactly 0: 0.0%   |   >= 90: 19.6%   |   >= 99: 5.7%
- deciles: 0, 25, 39, 50, 60, 68, 76, 83, 89, 96, 100
- (a piled-up distribution means a constant predictor is already hard to beat, and the benchmark's `mean baseline` interval width will be small)

## Components

components per reaction: min 2, median 5, max 18
- unique agent: 1560
- unique reactant: 22902
- unique products: 19286 (reactions per product: median 1)

Role signatures (most common):

- 2696 reactions: {'agent': 3, 'reactant': 2}
- 2225 reactions: {'agent': 4, 'reactant': 2}
- 2203 reactions: {'agent': 2, 'reactant': 2}
- 1858 reactions: {'agent': 1, 'reactant': 2}
- 1498 reactions: {'agent': 5, 'reactant': 2}
- 1283 reactions: {'reactant': 2}

## Group structure (can OOD splits be built?)

| group_by | n groups | largest group | groups with >=20 rows | median group size | group split | group-val split |
|---|---|---|---|---|---|---|
| scaffold | 11294 | 8.6% | 44 | 1 | yes | yes |
| product | 19286 | 0.0% | 0 | 1 | yes | yes |
| role:agent | 8338 | 7.4% | 83 | 1 | yes | yes |

Read: a group split needs a few dozen groups, none dominating (largest group well under 20 % of rows); `group-val` needs more groups still.
