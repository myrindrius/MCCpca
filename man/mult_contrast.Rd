% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mult_contrast.R
\name{mult_contrast}
\alias{mult_contrast}
\title{Multivariate Mean Comparison Test}
\usage{
mult_contrast(data, groups, contraste)
}
\arguments{
\item{data}{a data matrix}

\item{groups}{a vector or factor object giving the group for the
corresponding elements of data matrix.}

\item{contraste}{a contrast vector or contrast matrix where the sum of each
row is zero.}
}
\description{
This function applies a multivariate hypothesis test to compare means between groups
 according to a given vector or matrix of contrast.
}
\value{
A list of class "mult_contrast" containing:

\item{xbar.ip}{the matrix of mean vectors of wich group}
\item{delta}{the linear combination of the means}
\item{var.delta}{the variance's matrix of delta}
\item{T2}{the Hotelling t² value}
\item{wilks}{the Wilk's Lambda value}
}
\author{
Gualberto Segundo Agamez Montalvo; Bruno Marinho Estevam de Oliveira
}
