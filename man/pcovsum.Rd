\name{pcovsum}
\alias{pcovsum}
\title{Print covariate summary Latex}
\usage{
pcovsum(data, covs, maincov = NULL, numobs = NULL, nonpar = F, TeX = F)
}
\arguments{
  \item{data}{dataframe containing data}

  \item{covs}{character vector with the names of columns to
  include in table}

  \item{maincov}{covariate to stratify table by}

  \item{numobs}{named list overriding the number of people
  you expect to have the covariate}

  \item{nonpar}{boolean indicating if you want
  nonparametric (kruskal wallis) or parametric (anova)
  tests for continuous covariates. Default is FALSE.}

  \item{TeX}{boolean indicating if you want to be able to
  view extra long tables in the LaTeX pdf. If TeX is T then
  the table will not convert properly to docx}
}
\description{
Returns a dataframe corresponding to a descriptive table
}
\keyword{print}

