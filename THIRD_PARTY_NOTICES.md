# Third-party notices

pgphase is distributed under the MIT License (see `LICENSE`). It contains or
builds against the following third-party code, under the MIT License unless
stated otherwise.

## Code derived from longcallD

The BAM candidate collection, read-haplotype assignment, noisy-region MSA and
chunk stitching in `src/collect_var.cpp`, `src/collect_phase.cpp`,
`src/collect_phase_noisy.cpp`, `src/align.cpp` and `src/bam_digar.cpp` are a C++
port of longcallD (<https://github.com/yangao07/longcallD>).
Functions that keep their upstream names cite the upstream source location.

    MIT License
    Copyright (c) 2024- Dana-Farber Cancer Institute

## Files vendored in `src/`

| File | Origin | Copyright |
|---|---|---|
| `khash.h`, `kvec.h` | klib (<https://github.com/attractivechaos/klib>) | Copyright (c) 2008-2011 Attractive Chaos |
| `cgranges.c`, `cgranges.h` | cgranges (<https://github.com/lh3/cgranges>) | Copyright (c) 2019 Dana-Farber Cancer Institute |
| `kalloc.c`, `kalloc.h`, `kdq.h`, `sdust.c`, `sdust.h` | minimap2 (<https://github.com/lh3/minimap2>) | Copyright (c) 2018- Dana-Farber Cancer Institute |

## Git submodules in `third_party/`

| Path | Project | Licence file |
|---|---|---|
| `WFA2-lib` | WFA2-lib, Copyright (c) 2017 Santiago Marco-Sola | `third_party/WFA2-lib/LICENSE` |
| `abPOA` | abPOA, Copyright (c) 2020 Yan Gao | `third_party/abPOA/LICENSE` |
| `edlib` | edlib, Copyright (c) 2014 Martin Šošić | `third_party/edlib/LICENSE` |
| `minimap2` | minimap2, Copyright (c) 2018- Dana-Farber Cancer Institute (evaluation only) | `third_party/minimap2/LICENSE.txt` |
| `catch2` | Catch2 (tests only), Boost Software License 1.0 | upstream |

Each submodule's licence travels with its source.

## pgcat (`tools/pgcat`)

pgcat is built by cargo from these crates (versions pinned in
`tools/pgcat/Cargo.lock`), each with its licence in its source:

| Crate | Project | Licence |
|---|---|---|
| `gbz`, `simple-sds` | Jouni Sirén (<https://github.com/jltsiren/gbwt-rs>, <https://github.com/jltsiren/simple-sds>) | MIT |
| `rayon` | <https://github.com/rayon-rs/rayon> | MIT or Apache-2.0 |
| `flate2` | <https://github.com/rust-lang/flate2-rs> | MIT or Apache-2.0 |
| `zstd`, `zstd-safe` | <https://github.com/gyscos/zstd-rs> | MIT |
| `zstd-sys` | the same, bundling zstd (<https://github.com/facebook/zstd>) | BSD-3-Clause |

pgphase itself links the system zstd library (BSD-3-Clause) to read `.pgsites`
catalogs.
