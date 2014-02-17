\name{covsum}
\alias{covsum}
\title{Get covariate summary dataframe}
\usage{
covsum(data, covs, maincov = NULL, numobs = NULL, nonpar = F,
  markup = T, sanitize = T, nicenames = T)
}
\arguments{
  \item{data}{dataframe containing data}

  \item{covs}{character vector with the names of columns to
  include in table}

  \item{maincov}{covariate to stratify table by}

  \item{numobs}{named list overriding the number of people
  you expect to have the covariate}

  \item{nonparametric}{boollean indicating if you want
  parametric (anova) or nonparametric (Kruskal-Wallis)
  p-values for continuous variables.}

  \item{markup}{boolean indicating if you want latex
  markup}

  \item{sanitize}{boolean indicating if you want to
  sanitize all strings to not break LaTeX}

  \item{nicenames}{booling indicating if you want to
  replace . and _ in strings with a space}
}
\description{
Returns a dataframe corresponding to a descriptive table
}
\keyword{dataframe}

