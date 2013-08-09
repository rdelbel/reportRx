\name{crrRx}
\alias{crrRx}
\title{fit crr model}
\usage{
  crrRx(f, data)
}
\arguments{
  \item{f}{formula for the model. Currently the formula
  only works by using the name of the column in a
  dataframe. It does not work by using $ or [] notation.}

  \item{data}{dataframe containing data}
}
\description{
  Wrapper function to fit fine and gray competing risk
  model using function crr from package cmprsk
}
\keyword{model}

