# Methods represented by this release

The accompanying paper contains the full screening and implementation appendix.
This page defines the operations actually reproduced by the repository.

## Units and samples

Saved points and SEs use the original author-sample outcome standard deviation.
The denominator is common across estimators and is recomputed within the
upstream bootstrap. Baseline TWFE retains eligible source rows, including
always-treated observations; the approved omission of untreated reversal cells
after first adoption remains in force. CSA and BJS retain their own
identifiable targets. An aligned TWFE coordinate uses the exact union of rows
used by that modern estimator, including relevant baseline/control/nuisance
rows. Equal rows do not imply equal estimands or weights.

The source release has 351 specifications (323 scalar and 28 post event-study
averages). A dynamic comparison averages registered bins over common
structurally estimable horizons. The full joint regression is retained upstream.
Pre averages enter as separate points alongside post averages, never pooled
with them. They remain diagnostics, not treatment effects or an omnibus
pre-trend test. Comparison CSA retains its universal untreated baseline;
native adjacent-period pre contrasts are not substituted. Ordinary BJS pre
averages use an untreated-only joint lead regression, except for the preserved
mixed-reference contrasts. Aligned TWFE/BJS pre pairs can consequently be
algebraically identical. One
Study 19 (Muñoz) event window differs between sample families. Published mixed
reference bins and the explicit unrestricted fixed-effect exception for
Study 8 remain as documented in the paper.

The two added banking specifications in Colonnelli and Prem (RESTUD, 2022)
retain municipality and calendar-year fixed effects and municipality clustering.
The fourteen added public-source specifications in Cantoni and Pons (QJE, 2021)
retain their registered state-year designs, source-specific samples and controls,
and state clustering. Original native estimates remain separate from comparison
estimates where reversal removal changes the sample. Each study has a shared
positive-weight bootstrap across its specifications and estimator/sample
coordinates. The static additions have no manufactured pre averages or event
coefficient vectors.

## Covariance and paired inference

There are 999 upstream positive-weight bootstrap draws per declared block.
Independent exponential innovations were normalized to mean one over the
recorded dependence frame and shared across compatible specifications and
estimators. Fits, aggregation, and original-sample SDs were updated per draw.
Only **estimate vectors**, not observation or cluster-weight vectors, are
distributed here. Covariance is computed around draw means with divisor 998.

For a modern-minus-TWFE difference, the variance is
`V_modern + V_twfe - 2 * Cov(modern, twfe)`. Direct paired draw differences
provide a stable calculation near algebraic identities. Printed native SEs
and the joint bootstrap SEs are different quantities. No cross-block
covariance is supplied or assumed. Numerical draw completion does not
establish confidence-interval coverage.

## Concentration screen

The saved design metrics use net original-observation outcome coefficients:
`L_c = sum(a_i**2 for i in cluster c) / sum(a_i**2)` and
`G_eff = 1 / sum(L_c**2)`. This is an equal-observation-variance design
diagnostic, not an unrestricted cluster variance decomposition.

`config/screen.json` is applied anew when figures are generated:

- Exclude a pair if either estimate has `G_eff < 5` or `Lmax > 0.5`.
- Otherwise mark it fragile if either estimate or its nonidentity difference
  has `G_eff < 10` or `Lmax > 0.25`, or an earlier inference warning exists.
- An algebraically zero difference is not itself a concentration failure.
- Boundary tolerance is `1e-10`.

The concentration cutoff was selected after inspecting outliers. It is not
represented as a prespecified confirmatory screen. The repository reproduces
the inclusion rule from saved metrics; it does not rebuild those metrics from
microdata.

## Figure construction

`twfe-se` divides both estimates by that pair's saved TWFE bootstrap SE.
`outcome-sd` plots the stored estimates directly. The modern vertical coordinate
in the first presentation is not its own t statistic. The vertical gap is not
the covariance-aware paired-difference t statistic.

Circles identify static specifications; squares identify event-study averages,
using the same shape for pre and post. Filled markers require an unflagged paired absolute t statistic strictly above
1.96. Hollow markers are the other unflagged pairs. Fragile pairs use crosses.
Coordinates outside the chosen bounds are clipped independently; diamonds
identify capped pairs; capped event squares receive an additional diamond
outline. Flagged event points retain a square outline behind the cross.
Sign comparisons use uncapped estimates. Plotted
distances can therefore differ from actual gaps.

Shares give equal weight to every retained pair, including pre averages and fragile points in
the denominator. Significant-plus-reversal requires opposite nonzero point
estimate signs and a filled marker; it does not require two individually
significant coefficients. Tests are pointwise without multiplicity adjustment.
Specifications are dependent within studies and are not given equal study
weights. Colors and anonymous IDs are retained across panels and scales within
a release. Anonymous numbering follows journal publication year, with source
registry order breaking ties; earlier-year additions can renumber later studies
between releases. DOI/source-specification mappings are the persistent study
and specification identities.
`--summaries treatment` or `--summaries pre` permits separate descriptive views;
the default is both, and the export also gives annotation counts by period.
