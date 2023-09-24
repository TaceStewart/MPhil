# p_calculator takes in a dataframe df of reefs and calculates average
#   probability of the reef being in a recovered state, then returns the df with 
#   these new values. 
p_calculator <- function (reef_df, mgmt_benefit) {
  reef_df <- mutate(reef_df, 
                    r_single_unmgd = prob_s_recov,
                    r_single_mgd = 0, 
                    r_comp_unmgd = prob_c_recov,
                    r_comp_mgd = 0,
                    pr_recov_sing_unmgd = 0,
                    pr_recov_sing_mgd = 0,
                    pr_recov_comp_unmgd = 0,
                    pr_recov_comp_mgd = 0)
  
  # For each reef,
  for (reef in 1:nrow(reef_df)) {
    reef_df$r_single_mgd[reef] <- ifelse(is.na(as.numeric(reef_df$prob_s_recov[reef])),
                                         mean(as.numeric(reef_df$prob_s_recov), na.rm = TRUE) * 
                                           (1 + mgmt_benefit), 
                                         min(as.numeric(reef_df$prob_s_recov[reef]) * 
                                               (1 + mgmt_benefit), 1))
    reef_df$r_comp_mgd[reef] <- ifelse(is.na(as.numeric(reef_df$prob_c_recov[reef])),
                                       mean(as.numeric(reef_df$prob_c_recov), na.rm = TRUE) * 
                                         (1 + mgmt_benefit), 
                                       min(as.numeric(reef_df$prob_c_recov[reef]) * 
                                             (1 + mgmt_benefit), 1))
    
    # Calculate probability of recovered state for area if managed in single model
    # Prob disturbance at area
    di <- reef_df$prob_s_dist[reef] %>%
      as.numeric()
    
    # Prob Impacted given disturbance
    etai <- reef_df$prob_s_impact[reef] %>%
      as.numeric()
    if(is.na(etai)) {
      etai <- 0
    } 
    
    # Prob recovering post-disturbance
    ri <- reef_df$r_single_mgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(reef_df$r_single_mgd), na.rm = TRUE),
             reef_df$r_single_mgd[reef]) %>%
      as.numeric()
    nrows <- 2
    ncols <- 2
    A_stewart_num <- matrix(c((1-di*etai), ri, 
                              di*etai, (1-ri)),
                            nrows, ncols, byrow = TRUE)
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values  # eigenvalues
    ind <- which.max(D_stewart_num)
    reef_df$pr_recov_sing_mgd[reef] <- V_stewart_num[1,ind]/
      (sum(V_stewart_num[,ind])) # scale the first (i.e. 'recov') entry
    
    # Calculate probability of recovered state for area if not managed in single model
    ri <- reef_df$r_single_unmgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(reef_df$r_single_unmgd), na.rm = TRUE),
             reef_df$r_single_unmgd[reef]) %>%
      as.numeric()
    A_stewart_num <- matrix(c((1-di*etai), ri, 
                              di*etai, (1-ri)),
                            nrows, ncols, byrow = TRUE)
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values  # eigenvalues
    ind <- which.max(D_stewart_num)
    reef_df$pr_recov_sing_unmgd[reef] <- V_stewart_num[1,ind]/
      (sum(V_stewart_num[,ind])) # scale the first (i.e. 'recov') entry
    
    # Calculate probability of recovered state for area if managed
    # Prob disturbance in area
    d1 <- reef_df$prob_s_dist[reef] %>%
      as.numeric()
    d2 <- reef_df$prob_c_dist[reef] %>%
      as.numeric()
    # Prob Impacted given disturbance
    eta1 <- reef_df$prob_s_impact[reef] %>%
      as.numeric()
    if(is.na(eta1)) {
      eta1 <- 0
    } 
    eta2 <- reef_df$prob_c_impact[reef] %>%
      as.numeric()
    if(is.na(eta2)) {
      eta2 <- 0
    } 
    # Prob recovering post disturbance
    r1 <- reef_df$r_single_mgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(reef_df$r_single_mgd), na.rm = TRUE),
             reef_df$r_single_mgd[reef]) %>%
      as.numeric()
    r2 <- reef_df$r_comp_mgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(reef_df$r_comp_mgd), na.rm = TRUE),
             reef_df$r_comp_mgd[reef]) %>%
      as.numeric()
    nrows <- 3
    ncols <- 3
    A_stewart_num <- matrix(c((1-d1*eta1)*(1-d2*eta2),  r1,    r2, 
                              d1*eta1*(1-d2*eta2),      1-r1,  0, 
                              d2*eta2,                  0,     1-r2),
                            nrows, ncols, byrow = TRUE)
    eigs <- eigen(A_stewart_num) 
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values  # eigenvalues
    ind <- which.max(D_stewart_num)
    reef_df$pr_recov_comp_mgd[reef] <- V_stewart_num[1,ind]/
      (sum(V_stewart_num[,ind])) # scale the first (i.e. 'recov') entry
    
    # Calculate probability of recovered state for area if not managed
    r1 <- reef_df$r_single_unmgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(reef_df$r_single_unmgd), na.rm = TRUE),
             reef_df$r_single_unmgd[reef]) %>%
      as.numeric()
    r2 <- reef_df$r_comp_unmgd[reef] %>%
      {is.na(.) || is.nan(.) || . == "NaN"} %>%
      ifelse(mean(as.numeric(reef_df$r_comp_unmgd), na.rm = TRUE),
             reef_df$r_comp_unmgd[reef]) %>%
      as.numeric()
    A_stewart_num <- matrix(c((1-d1*eta1)*(1-d2*eta2),  r1,    r2, 
                              d1*eta1*(1-d2*eta2),      1-r1,  0, 
                              d2*eta2,                  0,     1-r2),
                            nrows, ncols, byrow = TRUE)
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values  # eigenvalues
    ind <- which.max(D_stewart_num)
    reef_df$pr_recov_comp_unmgd[reef] <- V_stewart_num[1,ind]/
      (sum(V_stewart_num[,ind])) # scale the first (i.e. 'recov') entry
  }
  
  return(reef_df)
}
