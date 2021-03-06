\name{term_similarity_from_MSigDB}
\alias{term_similarity_from_MSigDB}
\title{
Similarity between MSigDB terms
}
\description{
Similarity between MSigDB terms
}
\usage{
term_similarity_from_MSigDB(term_id, category = NULL, subcategory = NULL, ...)
}
\arguments{

  \item{term_id}{A vector of MSigDB gene set names.}
  \item{category}{E.g., 'C1', 'C2', pass to \code{\link[msigdbr]{msigdbr}}.}
  \item{subcategory}{E.g., 'CGP', 'BP', pass to \code{\link[msigdbr]{msigdbr}}.}
  \item{...}{Pass to \code{\link{term_similarity}}.}

}
\value{
A symmetric matrix.
}
\examples{
# There is no example
NULL

}
