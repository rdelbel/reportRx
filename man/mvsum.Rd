\name{mvsum}
\alias{mvsum}
\title{Get multivariate summary dataframe}
\usage{
  mvsum(model, data, markup = T, sanitize = T,
    nicenames = T)
}
\arguments{
  \item{model}{fitted model object}

  \item{data}{dataframe containing data}

  \item{markup}{boolean indicating if you want latex
  markup}

  \item{sanitize}{boolean indicating if you want to
  sanitize all strings to not break LaTeX}

  \item{nicenames}{booling indicating if you want to
  replace . and _ in strings with a space}
}
\description{
  Returns a dataframe corresponding to a univariate table
}
\keyword{dataframe}

