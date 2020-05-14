# plotIMC

Source code for plotIMC from Warren et al. (2020).  

## Web version

The app can be run locally, but an alternative live instance can be found online here:

http://mito.ncl.ac.uk/warren_2019/

## Local installation

To install the app locally, first install R, then install the required dependencies (particularly Shiny):

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
