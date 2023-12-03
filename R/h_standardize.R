#' Standardize
#'
#' Standardize a data matrix according to the given axis into mean of 0 and a standard deviation of 1.
#' @param data A data matrix need to be standardize
#' @param axis The axis needed to be standardized, default NA
#' @return standardized data
#' @export
standardize = function(data, axis=NA){
    return(h_transform(data, h_standardize,axis = axis,print=FALSE))
}

h_standardize = function(data){
    ## The same as the assert, assert if data is one-dimensional
    stopifnot(length(dim(data))==1 || is.null(dim(data)))
    N = length(data)
    data = clean_zs(data)
    mu = mean(data)
    #mad = median(data)
    ##Need to check this part. gamma(172) returns to Inf
    if(N <= 171) {
        c4 = sqrt(2/(N - 1)) * gamma(N/2) / gamma((N - 1)/2)
    } else {
        # if N is too big to lead a overflow
        c4 <- 1 - 1/(4*N) - 7/(32*N^2) - 19/(128*N^3)
    }
    std = sd(data)/c4
    # std = 1.4826* mad(data,constant =1)/c4
    # std = sqrt(sum(data-mu)^2)/(N-1))
    # Need check this
    if (std == 0) {
        std = .Machine$double.eps
        #std = 0.000000000000001
    }
    return((data-mu)/std)
}
