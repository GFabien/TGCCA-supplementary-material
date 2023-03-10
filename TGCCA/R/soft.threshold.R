# The function soft.threshold() soft-thresholds a vector such that the L1-norm constraint is satisfied.
# @param x A numeric vector.
# @param sumabs A numeric constraint on x's L1 norm.
#
# @return Returns a vector resulting from the soft thresholding of \eqn{x} given sumabs
# @keywords manip
#
soft.threshold <- function(x, sumabs = 1) {
  return(soft(x, BinarySearch(x, sumabs)))
}
