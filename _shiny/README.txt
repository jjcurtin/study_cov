# Covariate Dashboard

## Developer

Chris Janssen (cjanssen3@wisc.edu)

## Hosting

Deployed on Shinyapps.io at:
https://christopher-janssen.shinyapps.io/covariate-dashboard/

To redeploy after changes:
```r
library(rsconnect)
rsconnect::deployApp("path/to/recovered-dashboard", appName = "covariate-dashboard")
```

## Rendering

Built with R Shiny (R 4.5.2) using `bslib` for layout and `DT` for the interactive table. The app reads `results_tibble.csv` at startup — both `app.R` and `results_tibble.csv` are required to run.

To run locally:
```r
shiny::runApp("path/to/recovered-dashboard")
```
