\name{makedocx}
\alias{makedocx}
\title{Convert .TeX to .docx}
\usage{
  makedocx(dir, fname, pdwd, imwd = "")
}
\arguments{
  \item{dir}{full path of .TeX file directory}

  \item{fname}{.TeX file file name. Do not include
  extension}

  \item{pdwd}{full path to pandoc}

  \item{imwd}{full path to image magick. Only include if
  there is at least one graphic.}
}
\description{
  Covertes the knitr-compiled .TeX file to a .docx file
}
\keyword{print}

