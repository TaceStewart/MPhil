# p_calculator takes in a dataframe df of reefs and calculates average
#   probability of the reef being in a recovered state, then returns the df with 
#   these new values. The input df must have a `d_single`, `d_comp`, `r_single`,
#   and `r_comp` column
p_calculator <- function (df, mgmt_benefit) {
  df <- mutate(df, 
               r_single_unmgd = r_single,
               r_single_mgd = 0, 
               r_comp_unmgd = r_comp,
               r_comp_mgd = 0,
               pr_recov_sing_unmgd = 0,
               pr_recov_sing_mgd = 0,
               pr_recov_comp_unmgd = 0,
               pr_recov_comp_mgd = 0)
  
  # For each reef,
  for (reef in 1:nrow(df)) {
    df$r_single_mgd[reef] <- ifelse(is.na(as.numeric(df$r_single[reef])),
                                    mean(as.numeric(df$r_single), na.rm = TRUE) * 
                                           (1 + mgmt_benefit), 
                                    min(as.numeric(df$r_single[reef]) * 
                                          (1 + mgmt_benefit), 1))
    df$r_comp_mgd[reef] <- ifelse(is.na(as.numeric(df$r_comp[reef])),
                                  mean(as.numeric(df$r_comp), na.rm = TRUE) * 
                                         (1 + mgmt_benefit), 
                                  min(as.numeric(df$r_comp[reef]) * 
                                        (1 + mgmt_benefit), 1))
    
    # Calculate probability of recovered state for area if managed in single model
    # Prob disturbance at area
    di <- df$d_single[reef] %>%
      as.numeric()
    # Prob recovering post-disturbance
    ri <- df$r_single_mgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(df$r_single_mgd), na.rm = TRUE),
             df$r_single_mgd[reef]) %>%
      as.numeric()
    nrows <- 2
    ncols <- 2
    A_stewart_num <- matrix(c((1-di), ri, 
                              di, (1-ri)),
                            nrows, ncols, byrow = TRUE)
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values  # eigenvalues
    ind <- which.max(D_stewart_num)
    df$pr_recov_sing_mgd[reef] <- V_stewart_num[1,ind]/
      (sum(V_stewart_num[,ind])) # scale the first (i.e. 'recov') entry
    
    # Calculate probability of recovered state for area if not managed in single model
    ri <- df$r_single_unmgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(df$r_single_unmgd), na.rm = TRUE),
             df$r_single_unmgd[reef]) %>%
      as.numeric()
    A_stewart_num <- matrix(c((1-di), ri, 
                              di, (1-ri)),
                            nrows, ncols, byrow = TRUE)
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values  # eigenvalues
    ind <- which.max(D_stewart_num)
    df$pr_recov_sing_unmgd[reef] <- V_stewart_num[1,ind]/
      (sum(V_stewart_num[,ind])) # scale the first (i.e. 'recov') entry
    
    # Calculate probability of recovered state for area if managed
    # Prob disturbance in area
    d1 <- df$d_single[reef] %>%
      as.numeric()
    d2 <- df$d_comp[reef] %>%
      as.numeric()
    # Prob recovering post disturbance
    r1 <- df$r_single_mgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(df$r_single_mgd), na.rm = TRUE),
             df$r_single_mgd[reef]) %>%
      as.numeric()
    r2 <- df$r_comp_mgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(df$r_comp_mgd), na.rm = TRUE),
             df$r_comp_mgd[reef]) %>%
      as.numeric()
    nrows <- 3
    ncols <- 3
    A_stewart_num <- matrix(c((1-d1)*(1-d2), r1, r2, 
                              d1*(1-d2), 1-r1, 0, 
                              d2, 0, 1-r2),
                            nrows, ncols, byrow = TRUE)
    eigs <- eigen(A_stewart_num) 
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values  # eigenvalues
    ind <- which.max(D_stewart_num)
    df$pr_recov_comp_mgd[reef] <- V_stewart_num[1,ind]/
      (sum(V_stewart_num[,ind])) # scale the first (i.e. 'recov') entry
    
    # Calculate probability of recovered state for area if not managed
    r1 <- df$r_single_unmgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(df$r_single_unmgd), na.rm = TRUE),
             df$r_single_unmgd[reef]) %>%
      as.numeric()
    r2 <- df$r_comp_unmgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(df$r_comp_unmgd), na.rm = TRUE),
             df$r_comp_unmgd[reef]) %>%
      as.numeric()
    A_stewart_num <- matrix(c((1-d1)*(1-d2), r1, r2, 
                              d1*(1-d2), 1-r1, 0, 
                              d2, 0, 1-r2),
                            nrows, ncols, byrow = TRUE)
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values  # eigenvalues
    ind <- which.max(D_stewart_num)
    df$pr_recov_comp_unmgd[reef] <- V_stewart_num[1,ind]/
      (sum(V_stewart_num[,ind])) # scale the first (i.e. 'recov') entry
  }
  
  return(df)
}
