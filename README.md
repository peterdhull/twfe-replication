# Comparing TWFE, CSA, and BJS: Replication

This repository reproduces the baseline and aligned-sample figures
from "Two-Way Fixed Effects and Modern Difference-in-Differences: An Empirical Comparison," by Peter Hull. 
It also contains estimate-level bootstrap draws and
covariance blocks for further analysis. It contains **no observation-level
estimation data** and does not re-estimate models from raw observations.

## Quick start

Use Python 3.12 and the pinned environment:

```sh
python -m venv .venv
# Activate .venv using the command appropriate to your shell.
python -m pip install -r requirements.lock
python reproduce.py
python -m unittest discover -s tests -v
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1`.
On macOS/Linux, use `source .venv/bin/activate`.

The default command regenerates both sample families, both axis scales
(outcome standard deviations and TWFE standard errors), and capped/full views.
It writes plots and diagnostic tables under `outputs/`. No download, account,
API key, proprietary software, or original replication package is needed.

For the two capped figures used in the accompanying writeup:

```sh
python reproduce.py --scale twfe-se --view capped
```

The default includes both treatment summaries and pre averages. Use
`--summaries treatment` or `--summaries pre` to show just one set.

For companion figures labeled with study short titles and smaller publication
citations, use:

```sh
python reproduce.py --named
```

The named figures are written under `outputs/named/` by default. They use the
same estimates, screening, colors, chronological four-column legend, and
annotation counts as the anonymized figures. `--output` selects another output
directory for either presentation.

The Python entry point also supports `python -m didfigures` after installing
the local package with `python -m pip install --no-deps -e .`.

## What is included

| Item | Coverage |
| --- | --- |
| Implemented specifications | 351 across 24 studies |
| Treatment estimates before screening | 347 CSA and 351 BJS pairs per sample family |
| Retained treatment estimates | 251 CSA and 275 BJS pairs per family |
| Retained figures including pre averages | 273 CSA and 297 BJS pairs per family |
| Retained study/specification union | 22 studies, 285 specifications |
| Event-study treatment summaries | 28 pooled post-treatment averages |
| Additional available pre summaries | 22 across five studies, 22 CSA and 22 BJS pairs per family |
| Other treatment summaries | 323 scalar specifications before screening |
| Joint inference | 999 positive-weight bootstrap draws in each declared block |

Event studies enter as separate pre- and post-treatment averages, not once per
horizon; these are summaries of the same underlying specification. The default
figures combine the original scalar/post pairs with the available pre pairs,
applying the same screen separately to each pooled contrast. Circles identify
static specifications; squares identify both kinds of event-study averages.
The exact combined and pre-only counts are in
`reference/expected_summary.json`, and generated annotation tables break out
static, pre, and post summaries separately. The treatment-only counts for this
release are frozen separately in `reference/treatment_only_summary.json`.

Version 0.3.0 adds two banking specifications from Colonnelli and Prem's
*Corruption and Firms* (RESTUD, 2022) and fourteen public-source specifications
from Cantoni and Pons's *Strict ID Laws Don't Stop Voters* (QJE, 2021).
Both additions contain static targets. All existing pre averages are retained;
the additions create no event-study vectors or pre averages. Public numerical
tables include screen-excluded pairs, and the named provenance gives their
source table/column labels, recorded controls, fixed effects, and clustering.

All scalar estimates and their SEs are saved in the common original-sample
outcome-SD scale. The TWFE-SE presentation divides both plotted estimates by
that pair's TWFE bootstrap SE. It does not change paired inference.

## Repository layout

```text
config/               screening and figure policies
data/results/         complete estimate-level tables and anonymous metadata
data/covariance/      numeric estimate draws, covariance blocks, coordinate index
data/provenance/      study citations and source specification/coordinate labels
src/didfigures/       standalone screening, summaries, and plotting
tests/                numerical, covariance, and package checks
docs/                 data dictionary, methods, and release boundaries
paper/figures/        standalone figure assets (no manuscript source or PDF)
reference/            frozen expected summaries and release validation
outputs/              regenerated results (ignored by Git)
```

Read [the methods](docs/methods.md), [the data dictionary](docs/data_dictionary.md),
and [the release boundaries](docs/release_scope.md) before combining estimates.

The [validation record](reference/release_validation.json) documents the clean
environment checks. `release_manifest.json` records SHA-256 checksums for all
distributed files other than itself; generated outputs are not part of the
release archive. Git attributes preserve exact file bytes across platforms.

## Interpretation

These are descriptive estimator gaps, not observed causal bias. Baseline TWFE
retains its eligible sample, including always-treated observations. CSA and BJS
use their own identifiable targets. Aligned TWFE is refitted on the exact union
of observations used by the corresponding modern estimator; alignment does
not equalize treatment-effect weights or estimands.

Filled points use a paired difference test, including TWFE-modern covariance,
with a pointwise normal threshold of 1.96. Warning-marked points are retained as
crosses and cannot be filled. No multiplicity adjustment is applied. The
concentration screen is a design diagnostic, not a guarantee of valid inference.

Study and specification IDs link numerical tables to the named public
provenance in `data/provenance/`. The default figures retain their
anonymized labels. Original files, private correspondence, local source paths,
and private acquisition links are not included. Missing descriptive metadata
are left blank rather than inferred.

Study numbers follow journal publication year, with the source registry order
breaking ties. They are stable across views within this release; adding papers
in an earlier publication year can renumber later studies between releases.
Use DOI links and source specification labels when comparing releases.

## Figure assets and manuscript status

`paper/figures/` contains standalone PGF and vector PDF figure assets.
The draft manuscript LaTeX and compiled manuscript PDF are withheld from the
current release while the author makes final edits. Neither LaTeX nor Stata/R
is needed to regenerate the Python figures.

## Licensing and publication status

See [LICENSE.md](LICENSE.md). Public visibility does not assign a software or
data reuse license; that choice remains pending. The absence of observation-level data
does not transfer ownership of any underlying third-party inputs.
