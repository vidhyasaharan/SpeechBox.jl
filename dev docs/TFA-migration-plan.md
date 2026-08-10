# Migration plan: SpeechBox.jl → TimeFrequencyAnalysis.jl

**Goal.** Move the generic time-frequency core of SpeechBox.jl into a new standalone package,
[TimeFrequencyAnalysis.jl](https://github.com/vidhyasaharan/TimeFrequencyAnalysis.jl) (TFA), so it can be shared by
SpeechBox and future packages. SpeechBox becomes a speech-analysis layer that depends on TFA and re-exports its API,
so **existing SpeechBox user code keeps working unchanged**.

**Quality bar.** TFA must be *easier to maintain, inspect and understand* than SpeechBox. Concretely, code only moves
together with:

- a docstring on **every** function and struct, including internal helpers;
- a Documenter page that includes the docstring, organised by topic;
- code comments that explain intent (not restating the code);
- generic, non-speech naming and wording (a waveform is a waveform, not a speech signal);
- tests in TFA's own suite using **synthetic fixtures** (no WAV files, no WAV test dependency);
- one topic per source file, with a module file that reads as a table of contents.

No behaviour changes ride along with a move. Every commit in both repositories leaves the full test suite green.

---

## Package boundaries

| | TimeFrequencyAnalysis.jl | SpeechBox.jl (keeps) |
|---|---|---|
| Types | `waveform` (renamed from `speech_waveform`), `framed_signal`, `spectrum`, `timefreq`, `comp` | `pitch_timefreq`, `RAPT_*`, `GMM`, `Gaussian`, `Categorical`, `filter_coefs`(→ TFA in stage 2) |
| Analyses | `dft`, `magspec`, `specgram`, `periodogram`; correlations in stage 2 | `melfcc`, `vad`, `lpc`, `pitch` (spectral comb + RAPT), distortion measures |
| Utilities | framing, windows, frequency grids, `frqindex`, `time2nsamples`, signal generators, element-wise ops (`log`, `amp2db`, …) | `preemphasis`, peak finding (`findpeaks` family), stats/clustering layer |
| Dependencies | DSP, FFTW, LinearAlgebra, Random | TFA, Reexport, DSP, FFTW, LinearAlgebra, RecipesBase |

Back-compat: SpeechBox defines `const speech_waveform = waveform` and re-exports it, plus
`@reexport using TimeFrequencyAnalysis`, so the public SpeechBox API is a superset of what it was.

The stats/clustering layer (kmeans, GMM, statistics utilities) is out of scope for this migration; it may become a
third package later.

---

## Stages

### Stage 1 — minimal core (types, spectral analyses, utilities)  ✅ this stage first

The smallest self-consistent slice: the four data types, framing (which `specgram`/`periodogram` require), the
spectral analyses, and their supporting utilities.

**Moves to TFA** (source file in TFA → origin in SpeechBox):

| TFA file | Contents | From |
|---|---|---|
| `types.jl` | `waveform`, `framed_signal`, `spectrum`, `timefreq`, `comp` | signal_objects.jl (minus `pitch_timefreq`), SpeechBox.jl |
| `framing.jl` | `enframe(!)`, `view_frame`, `extract_frame`, `number_signal_frames`, `frame_energy`, padding helpers | speechframing.jl, internal_utilities.jl |
| `windows.jl` | `window` | utilities.jl |
| `spectral_analyses.jl` | `dft`, `magspec`, `specgram`, `periodogram`, `cexp`, `cexp_proj_matrix` | spectralanalyses.jl, internal_utilities.jl |
| `frequency_grids.jl` | `logfreq_array`, `linfreq_array`, `frqindex`, `findclosest` | internal_utilities.jl |
| `sampling.jl` | `time2nsamples`, `resample` | internal_utilities.jl, utilities.jl |
| `signal_generators.jl` | `white_noise`, `ar_process`, `impulse_train` | internal_utilities.jl |
| `element_ops.jl` | `log`/`log10`/`amp2db`/`pow2db` methods for `spectrum`/`timefreq` — rewritten as **explicit method definitions** (the `eval`-based generator in SpeechBox is deleted) | internal_utilities.jl |

**SpeechBox changes:**

- `Project.toml`: add TimeFrequencyAnalysis (with `[sources]` url) and Reexport; drop Random; version → 0.4.0.
- Module file: `@reexport using TimeFrequencyAnalysis`, `const speech_waveform = waveform`,
  `using TimeFrequencyAnalysis: findclosest` (the one internal SpeechBox still needs — kmeans++/GMM sampling).
- Delete moved code; `pitch_timefreq` moves into a new `pitch_objects.jl`; `utilities.jl` keeps only `preemphasis`;
  `internal_utilities.jl` keeps only the peak-finding family (`Δ`, `findpeaks`, `findpeaks_sorted`,
  `remove_nearest_peak!`).
- Plot recipes stay in SpeechBox for this stage (they attach to TFA types via dispatch; moved in stage 2).
- Tests: types/framing/window/grids/spectral/element-op/Float32-core test sets move to TFA and gain synthetic
  fixtures; SpeechBox keeps speech, correlation and clustering tests (which exercise the TFA core end-to-end).
- Docs: signal_objects.md, speechframing.md, spectralanalyses.md move to TFA (adapted); SpeechBox pages link to the
  TFA documentation.

### Stage 2 — analysis toolkit and plotting

- `correlations.jl` (`xcorr(!)`, `acorr(!)`, `acf`, `nacf`) → TFA. RAPT keeps using `nacf` via the re-export.
- `dsp_utilities.jl` (`filter_coefs`, `freq2θ`, `freq2z`, `H`, `Hmag`, `filter_resp`, `filter_magresp`) → TFA;
  SpeechBox's LPC imports what it needs.
- Plot recipes for `waveform`/`spectrum`/`timefreq` + `generate_ticks` → TFA (TFA gains the RecipesBase dependency);
  SpeechBox keeps the pitch-overlay recipes and imports `generate_ticks`.
- TFA version → 0.2.0; SpeechBox compat bump.

### Stage 3 — documentation polish, CI hardening, registration

- Complete TFA documentation site (tutorial-style index, examples with plots).
- Add a *downstream* job to TFA's CI that runs SpeechBox's test suite against each TFA PR
  (dev both into one environment, `Pkg.test("SpeechBox")` — the pattern verified below).
- Decide on registration: register TFA in the General registry first, then SpeechBox (deps must be registered
  before dependents), then remove the `[sources]` section. Until then, consumers use the two-step add below.
- Revisit marking semi-public helpers with `public` (needs Julia ≥ 1.11; TFA currently supports 1.10).

---

## Mechanics (verified empirically on Julia 1.12.6, 2026-08-10)

These were tested with throwaway packages before committing to the design:

1. **`[sources]` is not transitive.** `Pkg.add`ing SpeechBox by URL does *not* resolve unregistered TFA.
   Consumers must add TFA first:

   ```julia
   pkg> add https://github.com/vidhyasaharan/TimeFrequencyAnalysis.jl.git
   pkg> add https://github.com/unsw-edu-au/SpeechBox.jl.git
   ```

2. **`[sources]` works when the package itself is the active project** — SpeechBox CI resolves TFA from the URL
   with no workflow changes, and a fresh clone of SpeechBox instantiates.

3. **Local co-development:** keep both checkouts in `~/.julia/dev` and dev them into a shared environment:

   ```julia
   pkg> dev ~/.julia/dev/TimeFrequencyAnalysis ~/.julia/dev/SpeechBox.jl
   ```

   From that environment, `pkg> test SpeechBox` uses the **local TFA working tree** (verified: uncommitted TFA edits
   are picked up), and with Revise, edits to TFA are live in a SpeechBox session.

4. **Do not `Pkg.develop` inside the SpeechBox project itself** — it conflicts with the committed `[sources]` url
   ("`path` and `url` are conflicting specifications"). Always dev from the outer environment.

**Push order:** TFA must be pushed to GitHub *before* SpeechBox, otherwise SpeechBox CI and fresh clones resolve the
stale TFA `main`. Same rule for any TFA-breaking change: land TFA first.

**Versioning:** TFA follows semver from 0.1.0; every breaking change bumps the minor version even while unregistered
(compat is enforced against the cloned `Project.toml`), and SpeechBox's `[compat]` entry is updated in the same
change.

---

## Status log

| Date | Stage | Notes |
|---|---|---|
| 2026-08-10 | Pre-migration | SpeechBox fixed first: CI matrix, dead `timefreq` constructor, specgram frequency axis, `A_FACT`, and the `Float` alias replaced by parametric types generic over `AbstractFloat` (so the core moves over already generic). |
| 2026-08-10 | Stage 1 | **Done.** TFA v0.1.0 populated (8 source files, full docstrings, 97,531-test suite on synthetic fixtures, Documenter site with 6 pages, all building clean). SpeechBox v0.4.0 switched to the TFA dependency with `@reexport` + `speech_waveform` alias; suite green (5,108 tests) against the dev'ed TFA; docs build clean; benchmarks unchanged. Remember: **push TFA `main` before pushing SpeechBox** — SpeechBox CI resolves TFA from its GitHub url. |
| 2026-08-10 | Post-stage 1 | **TFA repo moved** to [github.com/vidhyasaharan/TimeFrequencyAnalysis.jl](https://github.com/vidhyasaharan/TimeFrequencyAnalysis.jl) (public; the old `unsw-edu-au` url is dead — no redirect). Both `[sources]` urls, the README/docs install instructions and links, and this plan updated; the `TFA_READ_TOKEN` auth steps deleted from CI.yml and Documentation.yml (public repo needs no token). SpeechBox itself stays at `unsw-edu-au/SpeechBox.jl`. **Flaky tests fixed:** the SpeechBox suite was unseeded, and two statistical GMM assertions (GMM_tests.jl:109 — exact match between two independent k-means runs; :185 — posterior argmax of sampled points) failed on ~7% of runs, including the first post-migration CI run on master. `test/setup.jl` now does `Random.seed!(2026)` (seed verified across the stochastic test files); `white_noise`/`ar_process` stay entropy-seeded by design, and the lpc tests that use them assert loose tolerances only. |
