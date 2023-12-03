#' OPTHT-SVD z score estimation
#'
#' A function to estimate and generate a z-score matrix after removing the noise using SVD and OPT, working as the outlier score.
#' p value will be calcualed according to the z-score and saved.
#' @param data_file the name of z_score for the read count
#' @import tools
#' @import dplyr
#' @return absolute path of the saved file
optht_svd_zs = function(data_file){
    out_basename = paste0(file_path_sans_ext(data_file),"-svd-optht")
    out_name_zs_optht = paste0(out_basename,"-zs.tsv")
    out_name_pv_optht = paste0(out_basename,"-pv.tsv")

    df = read.csv(data_file, sep = "\t")
    data_original = clean_zs(as.matrix(df))
    data = data_original

    svd_result = svd(data)
    U = svd_result$u
    s = svd_result$d
    VT = t(svd_result$v)
    s_dimension = optht(data, sv = s, sigma = NA)
    data_new = reverse_svd(data,U,s,VT,s_dimension)

    stds = std(data_new,axis = 1)[,1]
    stds = matrix(stds, ncol = 1)
    #zs_outrider__ = sweep((data-data_new),1,stds,"/")
    zs_outrider__ = (data-data_new)/as.vector(stds)
    data2 = standardize(zs_outrider__,axis =1)
    data2_df = as.data.frame(data2, row.names = rownames(df))
    colnames(data2_df) = colnames(df)
    write.table(data2_df,out_name_zs_optht, sep="\t")

    p_value = function(k){
        return(2*min(pnorm(k), (1-pnorm(k))))
    }
    p_value_matrix = apply(data2_df,c(1,2),FUN = p_value)
    adjusted_p_value = matrix(0, ncol = ncol(data), nrow = nrow(data))
    for ( i in 1:ncol(p_value_matrix)){
        adjusted_p_value[,i] = p.adjust(p_value_matrix[,i],method = "BY")
    }
    adjusted_p_value_df = as.data.frame(adjusted_p_value, row.names = rownames(df))
    colnames(adjusted_p_value_df ) = colnames(df)
    write.table(adjusted_p_value_df,out_name_pv_optht, sep="\t")

    abs_path_pv_optht = normalizePath(out_name_pv_optht)
    abs_path_zs_optht = normalizePath(out_name_zs_optht)
    return(list(abs_path_pv_optht = abs_path_pv_optht, abs_path_zs_optht = abs_path_zs_optht))
}
