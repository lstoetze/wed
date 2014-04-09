\name{wed}
\alias{wed}
\title{
wed
}
\description{
Allows to estimate euclidean weighted distance models of spatial voting
}
\usage{
wed(model_list, rep = 0, lik.func =lik, sep=FALSE ,se.A=FALSE)
}
\arguments{
  \item{model_list}{
Should ba a list that contains the data needed for estimation.  Rat contains a vector of voting decisions,  S is a matrix of selfplacment on the k political dimensions, P of parties position on the k dimensions, X is matrix containing respondent specific controls.
}
  \item{rep}{
Allows to specifiy if optimization algotithm should be started at other randomly drawn other points than zero. 
}
  \item{sep}{
Should seperability be allowed?
}
  \item{se.A}{
Should standard errors for the A matrix be simulated?
}

}


\author{
Lukas Stoetzer
}

\examples{

data(anes)
wed(model.data)
wed(model.data,sep=TRUE)


}

