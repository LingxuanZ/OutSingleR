#' Clean z-score
#'
#' Clean the z-score matrix through replacing the infinite values
#' @param data The z score matrix that need to be cleaned
#' @return Cleaned z score matrix
#' @export
clean_zs = function(data){
    temp = data ##Create a copy of data.
    temp[is.infinite(temp)] = 0 ##Replace the infinite with 0
    ## Replace the positive infinity with the max value of 7 and max absolute value in temp
    data[is.infinite(data) & data>0] = max(7,max(abs(temp)))
    ## Replace the negative infinity with the min value of 7 and -max absolute value in temp
    data[is.infinite(data) & data<0] = min(-7,-max(abs(temp)))
    return(data)
}
