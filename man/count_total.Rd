% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calba.R
\name{count_total}
\alias{count_total}
\title{Count Total Trees}
\usage{
count_total(gx, gy, r)
}
\arguments{
\item{gx}{A numeric vector of x-coordinates for the trees.}

\item{gy}{A numeric vector of y-coordinates for the trees.}

\item{r}{A numeric scalar for the radius parameter.}
}
\value{
A numeric vector containing the count of all trees within the radius for each focal tree.
}
\description{
This function counts the total number of trees within a given radius for each focal tree.
}
\examples{
sample_data <- data.frame(
  gx = runif(100, 0, 10),
  gy = runif(100, 0, 10)
)
count_total(
  gx = sample_data$gx,
  gy = sample_data$gy,
  r = 3
)

}
