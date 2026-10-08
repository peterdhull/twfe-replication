# Release boundaries

This is an estimate-level replication repository. It reproduces the screen,
annotation shares, and figures from frozen point estimates and design
diagnostics, and permits recalculation of SEs and covariances from estimate
draws. It does not claim independent recovery of those point estimates or
design diagnostics from observation-level source data.

## Included

- All available baseline and aligned treatment-effect estimator pairs, including
  pairs excluded by the figure screen.
- Stable-ID study/specification metadata, paired pooled pre comparisons, and
  explicit pre-comparison availability records.
- Study citations, official DOI links, recorded specification descriptions,
  source specification IDs, and descriptive covariance-coordinate labels.
- Numeric estimator draws, block-level covariance matrices, and coordinate
  mappings, with the covariance boundaries stated explicitly.
- Standalone Python code, pinned dependencies, policy configuration, tests,
  and the anonymized writeup.

## Excluded

- Individual/establishment/panel observations and prepared estimation panels.
- Raw replication archives, original study code, restricted source series,
  author correspondence, and account/access information.
- Original unit or cluster IDs, bootstrap cluster multipliers, support masks,
  raw outcome arrays, and the private filesystem/source-file crosswalk.
- Original filesystem paths, private download URLs or access tokens, and
  historical working-directory artifacts.
- The compiled manuscript PDF, withheld pending the author's final edits.
  Editable LaTeX and separate figure assets are included.

Some inputs used upstream were author-supplied and are not cleared for online
redistribution. They remain outside this repository. No claim that this release
provides complete raw-data replication is made.

The user authorized named study provenance in this repository. Default figure
labels remain anonymous, while `data/provenance/` supplies their bibliographic
mapping. This is not a blind-review repository. The separate private source-file
mapping and input access details remain outside this release.

The same numeric draw ID in different covariance blocks does not establish a
joint draw. Cross-block and between-study covariance are unavailable, not zero.
One study retains two different published dependence partitions. A future
observation-level release would require separate preparation, access, and
redistribution decisions.
