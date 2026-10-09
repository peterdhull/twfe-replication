# Data dictionary

CSV files are UTF-8. Blank cells mean unavailable/not recorded, not zero. Read
numerical CSVs with `pandas.read_csv(..., float_precision="round_trip")` to
preserve the released double-precision values. Do not round before calculating
contrasts or covariances.

## Numerical result tables

`data/results/comparisons.csv` has 1,396 rows: one row per active specification,
modern estimator, and sample family (698 rows per family). The unique scientific
key is `(specification_id, modern_estimator, family, summary_period)`; `comparison_id` is a
convenience row identifier. All screen-excluded pairs are included.

| Field(s) | Meaning |
| --- | --- |
| `study_id` | Integer ID 1--24, also used in the anonymized legends; fixed within this release |
| `specification_id` | `Sxx-Pyyy` active specification ID, fixed within this release |
| `family` | `baseline` or `aligned_samples` |
| `specification_type`, `summary_period` | Source adapter's type; treatment summary is `static` or `post` |
| `modern_estimator` | `csa` or `bjs` |
| `twfe_estimate`, `modern_estimate` | Point estimates in original-sample outcome-SD units |
| `twfe_positive_se`, `modern_positive_se` | SEs from all 999 estimate draws, including upstream outcome-scale variation |
| `covariance_twfe_modern` | Within-pair covariance in squared outcome-SD units |
| `gap`, `t_signed` | **TWFE minus modern** and its signed paired t statistic |
| `modern_minus_twfe` | Negative of `gap` |
| `se_gap`, `abs_t` | Direct paired-draw difference SE and absolute paired t statistic |
| `paired_normal_p_unadjusted` | Two-sided pointwise normal-reference p-value; no multiplicity adjustment |
| `paired_draws` | Number of paired draws (999 for each released comparison) |
| `original_twfe_estimate` | Original-sample rerun when supplied by that adapter; not uniformly a printed number |
| `published_estimator`, `twfe_is_published_estimator` | Native estimator designation/flag when recorded; constructed TWFE is distinguished |
| `published_native_estimate`, `published_native_positive_se` | Saved native-estimator coordinate and joint-bootstrap SE, when available |
| `published_estimate_raw`, `published_se_raw` | Recorded publication anchor in its original units |
| `native_b_raw`, `native_se_hc1_raw`, `native_se_cluster_raw` | Native rerun coefficient/inference in original units, when recorded |
| `author_native_rows`, `author_native_clusters` | Reported/reconstructed native sample and dependence counts |
| `original_target_rows`, `modern_target_rows` | Treated-target observation counts where recorded |
| `retained_target_share` | Saved supported treated analysis-weight share; not necessarily the ratio of row counts |
| `exact_union_rows` | Rows in the estimator-specific aligned TWFE union, when recorded |
| `reversal_rows_removed`, `published_vs_normalized_sample_mismatch` | Recorded differences between native and comparison samples |
| `bootstrap_cluster_frame` | Number of clusters in the recorded shared frame; **no cluster identifiers are included** |
| `registered_horizons`, `common_horizons`, `retained_horizon_count` | Original and comparison event terms/bins; blank for unrecorded/static fields |
| `common_estimator_set`, `horizon_shortened` | Estimator set used to choose common structural support and its shortening flag |
| `twfe_G_eff`, `modern_G_eff`, `gap_G_eff` | Effective-cluster design-concentration diagnostics |
| `twfe_Lmax`, `modern_Lmax`, `gap_Lmax` | Maximum cluster design-weight shares |
| `gap_structural_identity` | Algebraic identity flag; difference concentration is undefined for a zero coefficient vector |
| `point_check_max_error` | Upstream outcome-linear coefficient reconstruction check |
| `screen_prior_inference_warning` | Earlier warning preserved independently of the current concentration policy |
| `bjs_only_sparse_support`, `modern_se_degenerate_or_near_zero` | Source inference-support flags retained for audit |
| `original_modern_se_degenerate_or_near_zero`, `matched_gap_se_degenerate_or_near_zero` | Aligned-comparison auxiliary warning records |
| `screen_*`, `twfe_screen_*`, `modern_screen_*`, `gap_screen_*` | Frozen screening records; figure code **recomputes** decisions and fills from metrics plus prior warnings |
| `screen_exclusion_reason` | Which estimator axis failed the hard concentration screen |
| `published_fixed_effect_basis_exception`, `ancillary_placebo`, `native_design` | Explicit adapter/design qualifiers where recorded |
| `covariance_block_id` | Which joint distribution contains this pair |
| `twfe_coordinate_id`, `modern_coordinate_id` | Foreign keys into `data/covariance/coordinates.csv` |
| `twfe_coordinate_index`, `modern_coordinate_index` | **Zero-based** draw-column/covariance-row indices |

`studies.csv` supplies the 24 IDs, journal publication years, fixed colors, and
`draw_order_tiebreak`. The last field preserves the original layering of
overlapping points; it has no statistical interpretation. `specifications.csv`
contains the 351 active specifications and their scalar/post distinction.
Study numbers are assigned chronologically by journal year, with the source
registry order breaking ties. The version 0.3.0 additions are Study 12
(Cantoni and Pons, QJE 2021) and Study 14 (Colonnelli and Prem, RESTUD 2022).
IDs can therefore change between releases when earlier-year studies are added;
DOI and source specification labels in `data/provenance/` retain the identities.

`auxiliary_specifications.csv` documents 48 additional specifications retained
inside the complete source covariance blocks. They do **not** increase the
active 351-specification sample or appear in the treatment-effect figures.

`pre_comparisons.csv` adds 88 paired pre-treatment summaries: 22 event-study
specifications in five studies, with CSA and BJS in each sample family. Its
core estimate, uncertainty, concentration, and covariance-pointer fields use
the same definitions as `comparisons.csv`; `summary_period` is always `pre`.
Concatenate the two tables for the default figures, retaining `summary_period`
in every comparison key. A pre average is an additional summary of an existing
event specification, not another underlying regression specification.
`pre_comparison_availability.csv` records all 112 prospective pre pairs from
the 28 event specifications; 24 are unavailable, never zero-filled.

Pre rows record the actual common bins, the registered bins, and
`reference_not_pure_untreated_baseline`. Numeric bin labels need not denote
single years: Study 20 retains decade bins and the mixed reference `(-10,0]`.
Study 1's `le_m2` and Study 19's `-5` are open-ended bins. Human-readable
definitions are in `pre_window_definition`, `pre_reference_definition`,
and `pre_diagnostic_type`. Treatment-target
counts/shares are not copied into pre rows. Their blank values do not mean
zero support. The pre concentration vectors are reconstructed and checked
against the saved points; post-period concentration metrics are not reused.

## Pooled pre-treatment diagnostics

`pooled_pre.csv` has 393 saved diagnostic/availability records. Of these, 138
have a complete estimate and are linked to a covariance coordinate. Rows with
`available=False` are availability records, not zero estimates. The table
preserves family, coordinate/estimator role, pre-term counts/lists, estimates,
SEs, and recorded sample information. `unavailable_reason_code` uses generic
reason codes rather than private source text. `reference_not_pure_untreated_baseline`
flags mixed-reference interpretations. Raw-unit/native fields are explicitly
prefixed `raw_` or `native_`; missing fields remain blank.
The default figures select paired comparison coordinates through
`pre_comparisons.csv`, not every native/auxiliary row of this diagnostic table.
The sixteen static additions contribute 128 explicit unavailable records;
they do not change the existing available pre estimates or create event vectors.

## Covariance blocks

`data/covariance/blocks.json` lists 25 active blocks for 24 studies. A block's
NPZ file contains exactly three numeric arrays:

| Array | Shape | Meaning |
| --- | --- | --- |
| `point` | `(K,)` | Saved point estimates, not bootstrap means |
| `draws` | `(999, K)` | Joint **estimate** draws; columns index estimator coordinates |
| `covariance` | `(K, K)` | Sample covariance of the estimate draws (`ddof=1`) |

`coordinates.csv` labels all 2,957 coordinates with neutral IDs, study and
specification IDs, sample family, estimator, summary period/event term, point,
SE, and `used_by_active_comparison`. Complete source blocks preserve supporting
native, coefficient-vector, pre-diagnostic, and archived-alignment coordinates.
Use the explicit pointers in `comparisons.csv` for active figures; do not infer
membership from all coordinates in a block. The legacy flag
`used_by_active_comparison` refers to the original scalar/post table; the new
pre table supplies its own explicit coordinate pointers without altering the
saved arrays. In particular, family `archived_alignment_A`
coordinates are not a current figure family.

Coordinates can be algebraically identical, and the largest block has more
coordinates than draws. Covariance matrices may therefore be singular.
Positive definiteness is neither required nor implied. Do not invert them
without a rank-aware, substantively justified analysis.

Covariance is available only inside the same declared block. The same draw
number across different blocks does not mean joint resampling, including for
the study with two different dependence partitions. Missing cross-block
covariance must not be filled with zero.

Example: recover a paired gap SE directly from the released estimate draws:

```python
from pathlib import Path
import numpy as np
import pandas as pd

root = Path('.')
comparisons = pd.read_csv(root / 'data/results/comparisons.csv',
                          float_precision='round_trip')
row = comparisons.iloc[0]
with np.load(root / 'data/covariance' / f'{row.covariance_block_id}.npz',
             allow_pickle=False) as block:
    delta = (block['draws'][:, int(row.modern_coordinate_index)]
             - block['draws'][:, int(row.twfe_coordinate_index)])
    print(delta.std(ddof=1), row.se_gap)
```

## Public provenance

`data/provenance/studies.csv` gives titles, compact authors, journal/year,
citations, and official DOI links. `specifications.csv` maps active and
auxiliary release IDs to source specification labels and verified descriptive
metadata where available. `coordinate_labels.csv` maps neutral covariance IDs
to descriptive source coordinate labels. These labels do not change numerical
keys or the anonymized figure presentation. Blank descriptions are not inferred.

No public provenance table contains private author links, download tokens,
local filesystem paths, observation identifiers, or cluster identifiers.
