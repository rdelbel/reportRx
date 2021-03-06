\name{plotci}
\alias{plotci}
\title{Plot CI curve}
\usage{
  plotci(data, response, group = NULL, units = "months",
    main = "Viral Infections", CI = F, legpos = "topleft",
    xlim = NULL, outcomes = NULL)
}
\arguments{
  \item{data}{dataframe containing data}

  \item{response}{character vector or list of character
  vector. If a list it plot the '1' event for all outcomes
  on the same plot}

  \item{group}{string of the group want to stratify by}

  \item{units}{units of time}

  \item{main}{String corresponding to title}

  \item{CI}{Bool If True will plot CI and only the '1'
  event. if F will plot all events except for the final
  one}

  \item{legpos}{string indicating which position to put
  legend choies are "topright" etc}

  \item{xlim}{numeric vector corresponding to xlimits.
  Default is NULL}

  \item{outcomes}{character vector of the names of the
  different competing outcomes}
}
\description{
  Plots a CI curve. Currently not very powerful. Only plots
  a single curve
}
\keyword{print}

