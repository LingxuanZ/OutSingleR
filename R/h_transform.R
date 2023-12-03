#' Transform data
#'
#' Perform a user_identified transformation to the matrix.
#' @param a A matrix needed to be operated
#' @param transform_f The function need to be applied
#' @param axis The axis transform_f needed to be applied. Default: NULL
#' @param print_ Whether to print axis message, default FALSE
#' @return the size factor fore each column
#' @export
h_transform = function (a, transform_f, axis = NULL, print_= FALSE){
    a_new = array (0, dim = dim(a))

    if (is.null(dim(a)) || is.null(axis)){
        shape = dim(a)
        a = as.vector(a)
        a_new = transform_f(a)
        dim(a_new) = shape
        return(a_new)}

    else if (axis ==0){ ##Do column-wise
        for (c in 1:ncol(a)){
            if (print_){
                print(paste0("Processing row",c))}
            a_new[,c] = transform_f(a[,c])
        }
        return(a_new)
    }

    else if (axis ==1){ ##Row-wise
        for (r in 1:nrow(a)){
            if (print_){
                print(paste0("Processing row",r))
            }
            a_new[r,] = transform_f(a[r,])
        }
        return(a_new)
    }

    else{
        stop("Cannot happen")
    }
}
