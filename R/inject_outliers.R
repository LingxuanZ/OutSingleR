#' Inject outliers
#'
#' A function to inject outliers
#' @import tools
#' @import dplyr
#' @param fname the name of the gene count file to inject the outliers
#' @param outlier_frequency the number of outliers needed to be injected in each sample
#' @param z_score the magnitude of the outlies injected
#' @param outlier_type injected, 'o', overexpression, 'u', underexpressed, 'b', both
#' @export save the read table with injected outliers
inject_outliers = function (fname, outlier_frequency, z_score, outlier_type){

    name = file_path_sans_ext(fname)
    ext = paste0(".",file_ext(fname))
    id_ = sprintf("-wo-f%d-%s-z%.2f", outlier_frequency,outlier_type,z_score)

    df = read.csv(fname, sep="\t",row.names = 1, header = TRUE)
    sf_ = matrix(get_size_factor(df),nrow=1,ncol=ncol(df))
    data__ = as.matrix(df)
    data_sf__ = sweep(data__,2,sf_,"/")
    J = nrow(data_sf__)
    N = ncol(data_sf__)

    if(N <= 171) {
        c4 = sqrt(2/(N - 1)) * gamma(N/2) / gamma((N - 1)/2)
    } else {
        # if N is too big to lead a overflow
        c4 <- 1 - 1/(4*N) - 7/(32*N^2) - 19/(128*N^3)
    }

    mu__ = rowSums(data_sf__)/ncol(data_sf__)
    mu__ = matrix(mu__, nrow=nrow(data_sf__),ncol=1)

    l_ji__ = log2(sweep((data_sf__ +1),1,(mu__+1),"/"))
    l_j_ = rowSums(l_ji__)/ncol(l_ji__)
    l_j_ = matrix(l_j_,nrow =nrow(data_sf__),ncol=1)
    l_j_std_ = apply(l_ji__,1,sd)/c4
    l_j_std_[l_j_std_ == 0] <- 1e-15
    l_j_std_ = matrix(l_j_std_, nrow = nrow(data_sf__), ncol = 1)
    zs__ = sweep(sweep(l_ji__,1,l_j_,"-"),1,l_j_std_,"/")

    svd_result <- svd(zs__)
    U <- svd_result$u
    s <- svd_result$d
    VT = t(svd_result$v)

    s_dimension = optht(zs__, sv=s)
    zs_lr_optht__ = reverse_svd(zs__,U,s,VT,s_dimension)
    noise = zs__ - zs_lr_optht__
    zs_optht__ = matrix(NA,nrow=nrow(noise),ncol=ncol(noise))
    mu_optht_ = numeric(J)
    std_optht_ = numeric(J)

    for (j in 1:J){
        data_j_ = noise[j,]
        data_j_ = clean_zs(data_j_)
        mu = mean(data_j_)
        std = sd(data_j_)/c4
        if (std ==0){
            std = 1e-15
        }
        zs_optht__[j,] = (data_j_ -mu)/std
        mu_optht_[j] = mu
        std_optht_[j] = std
    }

    print(paste0("mu_optht_ ",mean(mu_optht_)))
    print(paste0("std_optht_ ",mean(std_optht_)))

    ##Create an outlier mask
    data_with_outliers = as.double(data__)
    outlier_mask = matrix(0,nrow=nrow(data__),ncol= ncol(data__),dimnames = dimnames(data__))

    for (i in 1:N){
        ois = sample(1:J,outlier_frequency)
        for (oi in ois){
            outlier_value = z_score
            if (outlier_type == "u"){
                outlier_value = -outlier_value
            }
            if (outlier_type =="b"){
                outler_value = sample(c(-1, 1), 1) * outlier_value
            }
            outlier_mask[oi,i] = outlier_value
        }
    }

    num_elements = nrow(zs_optht__) * ncol(zs_optht__)
    zs_optht_wo__ = matrix(rnorm(num_elements), nrow = nrow(zs_optht__), ncol = ncol(zs_optht__))

    zs_optht_wo__[outlier_mask != 0] = outlier_mask[outlier_mask != 0]
    #noise_wo__ = sweep(zs_optht_wo__,1,matrix(std_optht_,nrow=J,ncol=1,byrow=TRUE),"*")
    #noise_wo__ = sweep(noise_wo__,1,matrix(mu_optht_, nrow = J, ncol = 1, byrow = TRUE),"+")
    noise_wo__ = zs_optht_wo__ * std_optht_ + mu_optht_
    zs_wo__ = zs_lr_optht__ + noise_wo__

    l_ji_wo__ = zs_wo__ * as.vector(l_j_std_) + as.vector(l_j_)
    data_sf_wo__ = 2^l_ji_wo__ * as.vector(mu__ +1) -1

    data_wo__ = data_sf_wo__ *as.vector(sf_)
    data_wo__[data_wo__ <0] = 0

    outlier_mask_df = as.data.frame(outlier_mask)
    rownames(outlier_mask_df) = rownames(df)
    colnames(outlier_mask_df) = colnames(df)
    fname_new = paste0(name,id_,ext)
    fname_outliers_new = sub(paste0(ext, "$"), paste0("-omask",ext), fname_new)
    write.table(outlier_mask_df,fname_outliers_new)

    data_with_outliers = data_wo__
    data_with_outliers[is.nan(data_with_outliers)] = data__[is.nan(data_with_outliers)]
    data_with_outliers_df = as.data.frame(round(data_with_outliers))
    rownames(data_with_outliers_df) = rownames(df)
    colnames(data_with_outliers_df) = colnames(df)
    write.table(data_with_outliers_df,fname_new)
}
