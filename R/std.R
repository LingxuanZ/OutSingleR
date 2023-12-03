#' std
#'
#' Calculation the standardized standard deviation
#' @param data a data matrix.
#' @param axis The axis to perform the calculation.
#' @return Standardized standard deviation
#' @export
std = function(data, axis = NA){
    return(h_transform(data, std_, axis = axis))
}

std_ = function(data){
    N = length(data)
    data = clean_zs(data)
    mu = median(data)
    if(N <= 171) {
        c4 = sqrt(2/(N - 1)) * gamma(N/2) / gamma((N - 1)/2)
    } else {
        # if N is too big to lead a overflow
        c4 <- 1 - 1/(4*N) - 7/(32*N^2) - 19/(128*N^3)
    }
    std = 1.4826 * mad2(data)/c4
    if (std ==0){
        std = 1e-15}
    return(std)
}

mad2 = function(arr){
    arr = as.vector(t(arr))
    med = median(arr)
    mad = median(abs(arr-med))
    return(mad)
}
