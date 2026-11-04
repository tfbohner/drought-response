# Tree drought resistance and recovery

Code for a field study of how conifers in montane forests of California respond
to drought: how much their growth drops during a drought year, how fully they
recover afterward, and how species, region, tree size and competition shape
both. Each tree's annual growth comes from tree-ring cores, detrended and matched
to a drought index (12-month SPEI) for its site.

Published as: Bohner, T. & Diez, J. (2021). Tree resistance and recovery from
drought mediated by multiple abiotic and biotic processes across a large
geographic gradient. *Science of the Total Environment* 789: 147744.

*Teresa Bohner, PhD research*

## Where to start

- [`server model scripts/response recovery models final.R`](server%20model%20scripts/response%20recovery%20models%20final.R):
  the main Bayesian models (brms). Growth reduction and recovery are modeled
  with tree-level random effects; the first models let the residual variance
  differ by species (`sigma ~ Species`) and test those differences with
  `hypothesis()`. Later models add region and site structure with varying
  slopes, then tree size, competition and their interactions with drought.
- [`Maths.Rmd`](Maths.Rmd): the hierarchical model written out in notation.
- [`server model scripts/splineModelsDescription.Rmd`](server%20model%20scripts/splineModelsDescription.Rmd):
  a comparison of where terms go in brms models (population vs. group level)
  and how that affects estimates, convergence and run time.

## Workflow

1. `DataProcessing.Rmd`: detrends the tree-ring series (dplR), computes SPEI
   from site climate data, joins rings to field measurements (diameter,
   competition), and builds the drought-response and recovery datasets.
2. `server model scripts/`: model fitting, run on a server.
3. `Figures for Drought 2 MS.Rmd`: results, figures and tables for the
   manuscript; rendered to `docs/index.html`.
4. `prediction comparison.R`: checks predicted against observed growth for
   individual trees around drought years.

## Data

**Data are not included.** Tree-ring, field and climate data, and the fitted
model objects, stay off this repository. `DataProcessing.Rmd` reads the field
data from a Google Sheet: set `GOOGLE_ACCOUNT` and `FIELD_DATA_SHEET_ID` in your
`.Renviron`. Other scripts expect `Raw Data/`, `Processed Data/` and
`Brms models/` folders next to the code; all three are git-ignored.

## Notes

This code dates from 2019–2020 and was written for a dissertation chapter, so
some scripts assume objects created earlier in the same session.
