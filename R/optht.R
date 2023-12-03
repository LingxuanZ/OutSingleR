#' OPTHT
#'
#' A fast implementation of simple linear regression
#' @param data A data frame or matrix.
#' @param sv A vector of singular values from the matrix
#' @param sigma True variance if known, is NA by default
#' @param debug A boolean indicating if want verbose information
#' @importFrom stats integrate
#' @return An index for the hardthreshold for the singular values
#' @export
optht <- function(data, sv, sigma = NA, debug){
    if (is.array(data)) {
        m <- min(dim(data))
        n <- max(dim(data))
        beta <- m/n
    }
    if(is.na(sigma)){
        coef <- optimal_SVHT_coef_sigma_unknown(beta)
        cutoff <- coef * median(sv)
    }else{
        # Do some stuff... but I don't think we are ever in this case
        return(-1)
    }
    greater_than_cutoff <- which(sv > cutoff)
    k <- ifelse(length(greater_than_cutoff) > 0,
                max(greater_than_cutoff) + 1,
                0)
    return(k-1)
}
library(stats)
optimal_SVHT_coef_sigma_known <- function(beta){
    w <- (8 * beta) / (beta + 1 + sqrt(beta^2 + 14 * beta + 1))
    lambda_star <- sqrt(2 *(beta + 1) + w)
    return(lambda_star)
}

optimal_SVHT_coef_sigma_unknown <- function(beta){
    coef <- optimal_SVHT_coef_sigma_known(beta)
    omega <- coef / sqrt(median_marcenko_pastur(beta))
    return(omega)
}

mar_pas <- function(x, topSpec, botSpec, beta){
    # Implement Marcenko-Pastur distribution
    check <- ((topSpec - x) * (x - botSpec)) > 0
    return(
        ifelse(check,
               sqrt((topSpec - x) * (x - botSpec)) / (beta * x) / (2 * pi),
               0)
    )
}

median_marcenko_pastur <- function(beta, debug = F){
    # Compute median of Marcenko-Pastur distribution
    lb <- (1 - sqrt(beta))^2
    ub <- (1 + sqrt(beta))^2
    change <- T
    steps <- 1
    while(change & ((ub - lb) > 0.001)){
        if(debug){
            print(steps)
            print(lb)
            print(ub)
            steps <- steps + 1
        }
        change <- F
        x <- seq(lb, ub, length.out = 10)
        y <- 1 - sapply(x, FUN = median_integrate, lb = lb, ub = ub, beta = beta)
        if(any(y < 0.5)){
            lb <- max(x[y < 0.5])
            change <- T
        }
        if(any(y > 0.5)){
            ub <- min(x[y > 0.5])
            change <- T
        }
    }
    return((ub + lb) / 2)
}

median_integrate <- function(x, lb, ub, beta){
    ret <- stats::integrate(mar_pas, lower = x, upper = ub, topSpec = ub,
                            botSpec = lb, beta = beta)
    return(ret$value)
}
