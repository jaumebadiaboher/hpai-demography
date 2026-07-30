# HPAI demographic impact — tutorial & companion code

A step-by-step, annotated walkthrough of the integrated population model (IPM) behind:

> Badia-Boher, J.A., Schaub, M., Mollet, M., van Geneijgen, P., van der Jeugd, H.P.,
> Caliendo, V., & Kéry, M. (2026). Evaluating the demographic impacts of the highly
> pathogenic avian influenza panzootic. *Journal of Applied Ecology*, 63, e70234.
> [doi.org/10.1111/1365-2664.70234](https://doi.org/10.1111/1365-2664.70234)

The paper shows how to evaluate the past and future demographic impact of highly pathogenic avian influenza in long-lived species. This repo builds the Bayesian model to evaluate,
piece by piece — data, model structure, MCMC fitting, and the resilience analysis that
turns "the population crashed" into "here's how long recovery will take, and how far
short of its old size it may permanently fall."

## What's here

| File | What it is |
|---|---|
| `hpai_tutorial.qmd` | The tutorial source (Quarto). Renders to HTML with live R/JAGS code. |
| `hpai_tutorial.html` | The rendered tutorial — open this directly if you just want to read it. |
| `run_full_model.R` | Standalone script to reproduce the paper's exact MCMC run (HPAI scenario). ~6-7 hours. |
| `run_nohpai_model.R` | Companion script for the pre-outbreak (counterfactual) scenario, needed for the resilience comparison. |

**Want to just read it without downloading anything?** [View the rendered tutorial via htmlpreview](https://htmlpreview.github.io/?https://github.com/jaumebadiaboher/hpai-demographic-impact/blob/main/hpai_tutorial.html).

## Reproducing it yourself

The tutorial fetches the actual data and (optionally) the published posterior results
directly from the paper's Zenodo archive at render time — nothing is duplicated here.
Open `hpai_tutorial.qmd` in RStudio (or run `quarto render hpai_tutorial.qmd` from the
command line) and it'll walk through fetching everything it needs via the `zen4R`
package.

Requires: R, [JAGS](https://sourceforge.net/projects/mcmc-jags/), and the R packages
`jagsUI`, `zen4R`, `ggplot2`, `dplyr`, `gridExtra`, and [Quarto](https://quarto.org/).

## Data & code provenance

The underlying data and analysis code are archived openly on Zenodo under a CC-BY 4.0
license: [doi.org/10.5281/zenodo.17571734](https://doi.org/10.5281/zenodo.17571734).
This repo doesn't re-host that data — see `hpai_tutorial.qmd` for exactly how it's
fetched.

## Citation

If this tutorial or the underlying analysis is useful to you, please cite the paper above.
See `CITATION.cff` for a machine-readable citation (GitHub renders a "Cite this
repository" button from it automatically), and the Zenodo DOI above if you're
specifically reusing the code or data.

## License

CC-BY 4.0 — see `LICENSE`. Free to reuse and adapt with attribution.
