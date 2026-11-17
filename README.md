# Cal-Allocation

**Seed-reproducible sample-to-batch allocation for large-cohort
biomolecular detection with multi-feature balance and cohort integrity.**

Cal-Allocation partitions a batch of samples into detection sets (batches,
e.g. 35 samples per set) so that every set mirrors the feature composition of
the whole cohort — minimizing batch effects *at the experimental design
stage* rather than correcting them downstream. It was developed for
large-cohort metabolomics/proteomics studies on mass spectrometry platforms.

## Key features

1. **Multi-feature balancing.** Any number of categorical features (sex, age
   group, region, disease status, ...) are balanced simultaneously through a
   weighted, level-normalized objective
   `d(s) = Σ_f w_f · (1/L_f) · Σ_l |p̂(s,f,l) − p(f,l)|`. The primary
   criterion is the mean across sets, `mean_d = mean_s d(s)`; the dispersion
   `D = SD(d(1), …, d(S))` is reported as a secondary criterion.
2. **Seed-reproducible allocation route.** Samples with the same feature
   profile follow the same allocation path; randomness only orders
   interchangeable cohorts *within* a combination. With a fixed seed the
   result is reproducible bit-for-bit, and criterion values are highly stable
   across seeds (unlike pure shuffling methods).
3. **Local optimality under a stated neighbourhood.** The exchange phase
   terminates when no single exchange of two eligible equal-sized cohorts can
   further decrease the chosen objective. This is a local, neighbourhood-level
   guarantee; global optimality is not claimed.
4. **Cohort and subject integrity.** Samples with the same feature profile
   and the same `cohortID` are always placed in the same set. For longitudinal
   or repeated-measures studies, set `subject_block = TRUE` to make *all*
   samples from the same subject indivisible across visits, even when their
   feature profiles differ — a hard constraint that conventional
   randomization methods do not support.
5. **Exact set capacities by default.** With `allow_overshoot = FALSE` (the
   default), the algorithm stops with an error if an exact packing cannot be
   found, so every reported set size is a hard constraint rather than a
   post-hoc claim. Set `allow_overshoot = TRUE` only when a minimal-overflow
   fallback is explicitly acceptable.

## Algorithm overview

1. **Combo generation** — enumerate all feature-level combinations, discard
   empty ones, sort by size. Features with >10 distinct values are binned
   beforehand.
2. **Comb unit construction** — greedily build a set of combos covering every
   level of every feature ("full level coverage").
3. **Rotation** — each empty set draws one cohort from every combo of the
   unit in turn; exhausted combos are replenished. Placed samples are marked
   *essential* and protected from later exchanges.
4. **Greedy fill** — remaining cohorts are placed where they minimize the
   set deviation `d(s)`.
5. **Exchange optimization** — pairs of sets exchange equal-sized cohorts
   whenever the chosen criterion decreases (`objective_metric = "mean"` in the
   manuscript; `"sd"` retains the original dispersion objective), repeated
   until no improving exchange remains in that neighbourhood.
6. **Optional subject blocking** — with `subject_block = TRUE`, the algorithm
   first treats every subject as one indivisible unit (using subject size as
   the combo feature), then refines the partition with
   “one subject ↔ several subjects” subset exchanges that preserve each set's
   exact sample count. The result is a hard 0-split subject constraint.

## Installation

The core algorithm is pure base R (no package dependencies).

```r
source("R/cal_allocation.R")

result <- cal_allocation(
  data   = my_samples,          # data.frame: sampleID, cohortID, feature columns
  weight = list(disease = 1.0, sex = 0.5, center = 0.2, age_bin = 0.1),
  set_size = 35,                # samples per detection set
  objective_metric = "mean",    # primary criterion in the manuscript
  seed   = 20260813             # fixed seed -> bit-for-bit reproducible
)
head(result$assignment)         # original columns + key + new_set
```

For longitudinal/repeated-measures data, keep every subject in one set:

```r
result_blocked <- cal_allocation(
  data          = my_samples,
  weight        = list(disease = 1.0, sex = 0.5, center = 0.2, age_bin = 0.1),
  set.max.size  = rep(35, 156),
  cohort_col    = "cohortID",   # subject identifier
  subject_block = TRUE,         # all visits of a subject stay together
  seed          = 20260813
)
result_blocked$n_split_subjects  # must be 0
```

Input format: a data frame whose first column is a sample identifier and
which contains a cohort column (e.g. `cohortID`; repeated samples of the same
subject share one `cohortID`), followed by the feature columns. See
`demo/DEMO.DATA.xlsx` for a ready-to-run example.

The demo dataset is the simulated longitudinal cohort of the manuscript
(Cohort B): 5,500 samples from 2,827 subjects (541 with 2-10 longitudinal
visits), features `center` (C1-C6), `sex` (F/M), `disease` (D1-D9,
visit-specific) and `age_bin` (A1-A4), with weights disease 1.0 / sex 0.5 /
center 0.2 / age_bin 0.1. It is generated by `R/sim_cohort_b.R` with
`seed = 20260813`; running `cal_allocation()` on it with
`objective_metric = "mean"`, the same seed and
`set.max.size = c(rep(36, 40), rep(35, 116))` (156 sets) reproduces the
manuscript mean deviation (`mean_d = 0.019433`). With the original
`objective_metric = "sd"` the same data reproduce `D = 1.69 x 10^-5`.

## Repository layout

```
R/
  cal_allocation.R   # core five-phase algorithm
  metrics.R          # evaluation metrics (D, level RSD, coverage, chi2/Hellinger/max-d, cohort integrity)
  comparators.R      # complete / stratified / block randomization, minimization
  sim_lib.R          # simulation scenario definitions and data generator
  sim_cohort_b.R     # simulated longitudinal cohort (Cohort B) generator for the demo dataset
demo/
  DEMO.DATA.xlsx     # de-identified demo dataset (samples + weights sheets)
```

## Availability

A Shiny web application is available at
<https://calomics-best.shinyapps.io/metab-sample-allocation/>.

A manuscript describing the method is in preparation. Until it is published,
please cite this repository:

```
Chen R, et al. Cal-Allocation: deterministic and reproducible sample-to-batch
allocation for large-cohort biomolecular detection. (in preparation)
```

## License

MIT — see [LICENSE](LICENSE). The demo dataset is synthetic/de-identified
and provided for illustration only.
