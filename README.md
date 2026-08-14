# Cal-Allocation

**Deterministic, reproducible sample-to-batch allocation for large-cohort
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
   `d(s) = Σ_f w_f · (1/L_f) · Σ_l |p̂(s,f,l) − p(f,l)|`, minimizing the
   dispersion `D = SD(d(1), …, d(S))` across sets.
2. **Deterministic allocation route.** Samples with the same feature profile
   follow the same allocation path; randomness only orders interchangeable
   cohorts *within* a combination. With a fixed seed the result is
   reproducible bit-for-bit, and criterion values are highly stable across
   seeds (unlike pure shuffling methods).
3. **Quantifiable optimality.** The exchange phase terminates at a local
   optimum of `D` over the pairwise cohort-exchange neighbourhood: no single
   exchange of two eligible cohorts can further decrease `D`.
4. **Cohort integrity.** Repeated samples from the same donor/subject
   (longitudinal visits, aliquots) are treated as one indivisible unit and
   are always placed in the same set — a hard constraint that conventional
   randomization methods do not support.

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
   whenever the global criterion `D` decreases, repeated until no improving
   exchange remains.

## Installation

The core algorithm is pure base R (no package dependencies).

```r
source("R/cal_allocation.R")

result <- cal_allocation(
  data   = my_samples,          # data.frame: sampleID, cohortID, feature columns
  weight = list(center = 0.2, sex = 0.5, age = 0.1, group = 1.0),
  set_size = 35,                # samples per detection set
  seed   = 20260813             # fixed seed -> bitwise reproducible
)
head(result$assignment)         # original columns + key + new_set
```

Input format: a data frame whose first column is a sample identifier and
which contains a cohort column (e.g. `cohortID`; repeated samples of the same
subject share one `cohortID`), followed by the feature columns. See
`demo/DEMO.DATA.xlsx` for a de-identified example.

## Repository layout

```
R/
  cal_allocation.R   # core five-phase algorithm
  metrics.R          # evaluation metrics (D, level RSD, coverage, chi2/Hellinger/max-d, cohort integrity)
  comparators.R      # complete / stratified / block randomization, minimization
  sim_lib.R          # simulation scenario definitions and data generator
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
