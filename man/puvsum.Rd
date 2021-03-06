\name{puvsum}
\alias{puvsum}
\title{Print univariate summary LaTeX table}
\usage{
  puvsum(response, covs, data, type = NULL, strata = 1,
    TeX = F)
}
\arguments{
  \item{response}{string vector with name of response}

  \item{covs}{character vector with the names of columns to
  fit univariate models to}

  \item{data}{dataframe containing data}

  \item{type}{string indicating he type of univariate model
  to fit. The function will try and guess what type you
  want based on your response. If you want to override this
  you can manually specify the type. Options in clude
  "linear", "logistic", "coxph", "crr",
  "boxcox","logistic"}

  \item{strata}{character vector of covariates to stratify
  by. Only used for coxph and crr}

  \item{TeX}{boolean indicating if you want to be able to
  view extra long tables in the LaTeX pdf. If TeX is T then
  the table will not convert properly to docx}
}
\description{
  Returns a LaTeX table of the univariate summary
}
\keyword{dataframe}

