# plotIMC

Source code for plotIMC from Warren et al. (2020) ([v0.0.1](https://github.com/CnrLwlss/Warren_2019/releases/tag/v0.0.1)).

## Web version

The app can be run locally, but an alternative live instance can be found online here:

http://mito.ncl.ac.uk/warren_2019/

## Local installation

To install the app locally, first download and install [R](https://www.r-project.org/), then use R to install the required dependencies (particularly [Shiny](https://shiny.rstudio.com/)):

```R
install.packages(c("shiny","corrgram","data.table"))
```

## Starting the shiny app

To start the app, clone this repository onto your local machine, navigate to the by_patient directory, start an R session, load the shiny library:

```R
library(shiny)
```

then run the app:

```R
runApp()
```

plotIMC should then open in your web browser.

## Editing plotIMC

It is possible to adapt the plotIMC source code to run with a different dataset, however it does require adapting the code, and therefore requires understanding R and Shiny.  plotIMC does not have a GUI for uploading different datasets.
