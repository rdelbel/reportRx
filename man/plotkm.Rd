\name{plotkm}
\alias{plotkm}
\title{Plot KM curve}
\usage{
plotkm(data, response, group = 1, units = "months", show.nrisk = T,
  CI = F, HR = F, legend.pos = "bottomleft", hr.pos = "topright",
  cex = 1, title = "")
}
\arguments{
  \item{data}{dataframe containing your data}

  \item{response}{character vector with names of columns to
  use for response}

  \item{group}{string specifiying the column name of
  stratification variable}

  \item{legend.pos}{what position you want the legend to
  be. use NA for no legend}

  \item{hr.pos}{what position do you want the hr and
  p-value to be. Use NA for no hr and p-value}

  \item{units}{string specifying what the unit of time is
  use lower case and plural}

  \item{CI}{boolean to specify if you want confidence
  intervals}

  \item{HR}{boolean to specify if you want the hazard ratio
  from a cox model}

  \item{legend}{boolean to specify if you want a legend}

  \item{title}{title of plot}

  \item{show.nrisk}{boolean indicating if you want to show
  the number at risk}

  \item{cex}{controls the size of text in the legends}
}
\description{
This function will plot a KM curve with possible
stratification. You can specifyif you want a legend or
confidence bands as well as the units of time used.
}
\examples{
require(survival)
lung$sex<-factor(lung$sex)
plotkm(lung,c("time","status"))
plotkm(lung,c("time","status"),"sex")
}
\keyword{plot}

