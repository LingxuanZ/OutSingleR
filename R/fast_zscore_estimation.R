#' fast z score estimation
#'
#' A function to estimate and generate a z-score matrix for the read count
#' @import tools
#' @import dplyr
#' @param data_file the path of the .csv or .tsv file
#' @return absolute path of the saved file storing the z_scores data frame
fast_zscore_estimation = function(data_file){
    df = read.csv(data_file, sep="\t")
    ColName = names(df)
    if(is.character(df[,1])){
        row.names(df) = df[,1]
        df = df[,-1]
    }
    data = sweep(as.data.frame(lapply(df, as.double), stringsAsFactors = FALSE),2,get_size_factor(df),FUN="/")
    data_cleaned = data

    ZScoresResult = get_z_scores(data)
    z_scores=ZScoresResult$z_scores
    l_ji__ = ZScoresResult$l_ji__
    # print(which(abs(z_scores)>THRESHOLD))
    #ZScoresResult_2 = get_z_scores(data_cleaned,l_ji__1)
    #z_scores_2=ZScoresResult_2$z_scores
    #l_ji__2 = ZScoresResult_2$l_ji__

    fname = paste0(tools::file_path_sans_ext(data_file), "-fzse-zs.tsv")# Get the base file name without extension
    z_scores_df = as.data.frame(z_scores, row.names = rownames(df), col.names = colnames(df))
    #z_scores_df = bind_cols(rownames(z_scores_df),z_scores_df)
    #rownames(z_scores_df) = NULL
    #colnames(z_scores_df) = ColName
    write.table(z_scores_df, fname,sep="\t")
    #save_dfz_to_csv(z_scores_df, fname) # Save 2 dataframes: z_scores and their p-values  into CSV files
    abs_path = normalizePath(fname) # absolute path of the saved file
    return(abs_path)
}
