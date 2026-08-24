# Acceptance benchmarks

Two measurement-only scripts, promoted from the PoC's benchmark harness. Run
them from the repo root with the bench environment, which adds
DiscreteGlobalGrids and Geomorphometry as hard dependencies:

```
julia --project=examples/GeomorphometryNG.jl/bench examples/GeomorphometryNG.jl/bench/accept.jl 12
julia --project=examples/GeomorphometryNG.jl/bench examples/GeomorphometryNG.jl/bench/accept_rect.jl
```

- **`accept.jl <level> [full|sub]`** — the cell-grid matrix on the tutorial's
  IGeo7 region: TPI against DGG's raw streaming pass, `steepest_slope` on-demand
  vs `precompute`d, and `flowaccumulation` against Geomorphometry, with the four
  ACCEPT gates and a value-agreement line for each. Default level is 13.
- **`accept_rect.jl`** — the rectilinear regression: two 1000×1000 Float64 DEMs,
  one smooth (`steepest_slope`, `flow_direction`, PlaneFit, TPI) and one random
  (`flowaccumulation` vs Geomorphometry, directions compared code for code).
- **`setup.jl`** — `tutorial_elevation(level)`, the synthetic GLO-30 stand-in
  regridded to IGeo7. Regridded values are cached in `cache/*.jls` (ignored by
  git); `cache/` ships with levels 11 and 12, and any other level — including 13,
  whose cache is 62 MiB — is regenerated on first use, which takes a few minutes.

Reference numbers, 2026-08-24, Julia 1.12.6, 8 threads (`scratchpad/v3/accept_L13.log`,
`accept_rect.log`), measured on the PoC core these scripts were promoted from:

| level 13 (16,181,892 cells) | time | alloc |
| --- | --- | --- |
| TPI — raw `DGG.mapneighbors(Values())` | 503.6 ms | 61.7 MiB |
| TPI — `topographic_position_index` | 496.2 ms (0.99× baseline) | 61.7 MiB |
| `steepest_slope`, on-demand geometry | 3.182 s | 123.5 MiB (8.0 B/cell = output floor) |
| `steepest_slope`, `precompute`d | 1.242 s (build 2.173 s) | 123.5 MiB (build 370.4 MiB) |
| `flowaccumulation` — `GM.flowaccumulation(D8())` | 9.641 s | 2084.9 MiB |
| `flowaccumulation` — `flowaccumulation` | 3.117 s (0.32×) | 1994.9 MiB, 100.0000% identical direction codes |

Rectilinear, 1000×1000: `steepest_slope` 2.28 ms / 7.65 MiB; `flowaccumulation`
139.7 ms / 49.4 MiB vs Geomorphometry's 516.5 ms / 123.6 MiB (3.70×), directions
identical and `max rel acc diff = 0`.
