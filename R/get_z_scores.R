#' get z score
#' @param data A size J*N dataframe
#' @param l_ji_other__ A size *N dataframe
#' @return A list containing Z scores and l_ji__
#' @examples
#' data = read.csv("geneCounts.tsv",sep="\t")
#' get_z_scores(data)
#' @export

get_z_scores <- function(data, l_ji_other__=NULL){
    J = nrow(data)
    N = ncol(data)

    c4 <- tryCatch({
        sqrt(2 / (N - 1)) * gamma(N / 2) / gamma((N - 1) / 2)
    }, error = function(err) {
        # Handle the OverflowError
        1 - 1 / (4 * N) - 7 / (32 * (N^2)) - 19 / (128 * (N^3))
    })

    mu__ = rowMeans(data, na.rm = TRUE)
    mu__ = matrix(mu__, ncol = 1, nrow = nrow(data))

    l_ji__ = log2((data + 1) / (mu__ + 1))
    l_j_ = rowMeans(l_ji__, na.rm = TRUE)
    l_j_ = matrix(l_j_, ncol = 1, nrow = nrow(data))
    l_j_std_ = apply(l_ji__, 1, function(x) sd(x, na.rm = TRUE) / c4)
    #  l_j_std_[l_j_std_ == 0] = 0.001
    l_j_std_[l_j_std_==0] = 0.000000000000001
    l_j_std_ = matrix(l_j_std_,ncol = 1, nrow = nrow(data))

    if (is.null(l_ji_other__)) {
        l_ji_other__ = l_ji__
    }
    z_scores = (l_ji_other__ - l_j_) / l_j_std_

    return(list(z_scores=z_scores,l_ji__=l_ji__))
}
