#' Get size factor
#'
#' Calculation of the size factor for each column (sample) as the median value
#' of the non_zero normalized counts. The original read counts for each are normalized
#' to the geometric means of the counts.
#' @param df a dataframe with the gene counts.
#' @return the size factor fore each column
#' @export

get_size_factor = function(df){
    data_ori = as.matrix(df)
    N = ncol(data_ori)
    # Filter rows where all values are non-zero
    valid_rows = rowSums(data_ori !=0) == ncol(data_ori)
    data = data_ori[valid_rows,]

    #Initialize a count_normalized matrix
    counts_normalized = matrix (0, nrow = nrow(data), ncol=ncol(data))
    #Calcualte normalized counts for each gene
    for (j in 1:nrow(data)){
        row = data[j,]
        #Remove 0 values in each row
        counts = row[row != 0]
        if (length(counts)>0){
            #Geometric mean
            denominator = mp_gmean(counts)
        }
        else{
            denominator = 0
        }

        if (denominator != 0){
            counts_normalized[j,] = row / denominator
        }
        else{
            counts_normalized[j,] = rep(0,ncol(data))
        }
    }

    size_factors = numeric(ncol(data))
    for (i in 1:ncol(data)){
        column = counts_normalized[,i]
        column = column[column != 0]
        size_factors[i] = median(column)
    }
    return(size_factors)
}

mp_gmean = function(array){
    log_sum <- sum(log(array))
    n <- length(array)
    exp(log_sum / n)
}
