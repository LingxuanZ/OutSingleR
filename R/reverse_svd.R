#' Reverse SVD
#'
#' Calculate the read count matrix after removing the noises.
#' @param data The original gene count matrix
#' @param U the U matrix after SVD decomposition of the original read count
#' @param s the singular matrix after SVD decomposition of the original read count
#' @param VT the VT matrix after SVD decomposition of the original read count
#' @param s_dimension cut off value calculated by optht function
#' @return the reconstructed gene count matrix after removing noises.
#' @export
reverse_svd <- function(data,U, s, VT, s_dimension){
    U_new = U[,1:s_dimension]
    s_new = s[1:s_dimension]
    VT_new = VT[1:s_dimension,]
    return(U_new%*% diag(s_new)%*% VT_new)
}
