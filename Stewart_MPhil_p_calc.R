# p_calculator takes in a dataframe reef_data of reefs and calculates average
#   probability of the reef being in a recovered state, then returns the reef_data with
#   these new values.
p_calculator <- function(reef_data, mgmt_benefit) {
  reef_data <- mutate(reef_data,
                      r_single_unmgd = prob_s_recov,
                      r_single_mgd = 0,
                      r_comp_unmgd = prob_c_recov,
                      r_comp_mgd = 0,
                      pr_recov_sing_unmgd = 0,
                      pr_recov_sing_mgd = 0,
                      pr_recov_comp_unmgd = 0,
                      pr_recov_comp_mgd = 0
  )
  
  # Convert all columns except reef_name and sector to numeric
  cols_names <- c("prob_s_dist", "prob_c_dist", 
                  "prob_s_impact", "prob_c_impact", 
                  "prob_s_recov", "prob_c_recov", 
                  "r_single_unmgd", "r_single_mgd", 
                  "r_comp_unmgd", "r_comp_mgd", 
                  "pr_recov_sing_unmgd", "pr_recov_sing_mgd", 
                  "pr_recov_comp_unmgd", "pr_recov_comp_mgd")
  reef_data[cols_names] <- sapply(reef_data[cols_names], as.numeric)
  
  # For each reef,
  for (reef in 1:nrow(reef_data)) {
    sector_indx <- which(reef_data$sector == reef_data$sector[reef])
    if (is.na(reef_data$prob_s_recov[reef])) {
      # Get recovery rates in the management area/sector
      sector_s_recovs <- reef_data$prob_s_recov[sector_indx] %>%
        as.numeric()
      sector_s_recovs <- sector_s_recovs[!is.na(sector_s_recovs)]
      
      # Sample from the distribution of non-NA recovery rates in the management area/sector
      reef_data$r_single_unmgd[reef] <- ifelse(length(sector_s_recovs) == 0 || is.na(sd(sector_s_recovs)),
                                               rnorm(n = 1,
                                                     mean = mean(reef_data$prob_s_recov,
                                                                 na.rm = TRUE),
                                                     sd = sd(reef_data$prob_s_recov,
                                                             na.rm = TRUE)) %>%
                                                 as.numeric(),
                                               rnorm(n = 1, 
                                                     mean = mean(sector_s_recovs,
                                                                 na.rm = TRUE), 
                                                     sd = sd(sector_s_recovs,
                                                             na.rm = TRUE)) %>%
                                                 as.numeric())
      if (is.na(reef_data$r_single_unmgd[reef]) || 
          is.nan(reef_data$r_single_unmgd[reef])) {
        reef_data$pr_recov_sing_unmgd[reef] <- NA
        reef_data$pr_recov_comp_unmgd[reef] <- NA
        reef_data$pr_recov_sing_mgd[reef] <- NA
        reef_data$pr_recov_comp_mgd[reef] <- NA
        next
      }
    }
    reef_data$r_single_unmgd[reef_data$r_single_unmgd < 0] <- 0
    reef_data$r_single_unmgd[reef_data$r_single_unmgd > 1] <- 1
    reef_data$r_single_mgd[reef] <- 1 - mgmt_benefit * (1 - as.numeric(reef_data$r_single_unmgd[reef]))
    
    if (is.na(reef_data$prob_c_recov[reef])) {
      # Get recovery rates in the management area/sector
      sector_c_recovs <- reef_data$prob_c_recov[sector_indx] %>%
        as.numeric()
      sector_c_recovs <- sector_c_recovs[!is.na(sector_c_recovs)]
      
      # Sample from the non-NA recovery rates in the management area/sector
      # Sample from the distribution of non-NA recovery rates in the management area/sector
      reef_data$r_comp_unmgd[reef] <- ifelse(length(sector_c_recovs) == 0 || is.na(sd(sector_c_recovs)),
                                             rnorm(n = 1,
                                                   mean = mean(reef_data$prob_c_recov,
                                                               na.rm = TRUE),
                                                   sd = sd(reef_data$prob_c_recov,
                                                           na.rm = TRUE)) %>%
                                               as.numeric(),
                                             rnorm(1, 
                                                   mean = mean(sector_c_recovs,
                                                               na.rm = TRUE), 
                                                   sd = sd(sector_c_recovs,
                                                           na.rm = TRUE)) %>%
                                               as.numeric())
      if (is.na(reef_data$r_comp_unmgd[reef]) || 
          is.nan(reef_data$r_comp_unmgd[reef])) {
        reef_data$pr_recov_sing_unmgd[reef] <- NA
        reef_data$pr_recov_comp_unmgd[reef] <- NA
        reef_data$pr_recov_sing_mgd[reef] <- NA
        reef_data$pr_recov_comp_mgd[reef] <- NA
        next
      }
    }
    reef_data$r_comp_unmgd[reef_data$r_comp_unmgd < 0] <- 0
    reef_data$r_comp_unmgd[reef_data$r_comp_unmgd > 1] <- 1
    reef_data$r_comp_mgd[reef] <-  1 - mgmt_benefit * (1 - as.numeric(reef_data$r_comp_unmgd[reef]))
    
    # Calculate probability of recovered state for area if managed in single model
    # Prob disturbance at area
    di <- reef_data$d_tot[reef] %>%
      as.numeric()
    
    # Prob Impacted given disturbance
    etai <- reef_data$prob_s_impact[reef] %>%
      as.numeric()
    if (is.na(etai)) {
      sector_s_etas <- as.numeric(reef_data$prob_s_impact[sector_indx])
      etai <- sample(sector_s_etas[!is.na(sector_s_etas)], 1) %>%
        as.numeric()
    }
    
    # Prob recovering post-disturbance
    ri <- reef_data$r_single_mgd[reef] 
    nrows <- 2
    ncols <- 2
    A_stewart_num <- matrix(
      c(
        (1 - di * etai), ri,
        di * etai, (1 - ri)
      ),
      nrows, ncols,
      byrow = TRUE
    )
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values # eigenvalues
    ind <- which.max(D_stewart_num)
    reef_data$pr_recov_sing_mgd[reef] <- V_stewart_num[1, ind] /
      (sum(V_stewart_num[, ind],
           na.rm = TRUE)) # scale the first (i.e. 'recov') entry
    
    # Calculate probability of recovered state for area if not managed in single model
    ri <- reef_data$r_single_unmgd[reef]
    A_stewart_num <- matrix(
      c(
        (1 - di * etai), ri,
        di * etai, (1 - ri)
      ),
      nrows, ncols,
      byrow = TRUE
    )
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values # eigenvalues
    ind <- which.max(D_stewart_num)
    reef_data$pr_recov_sing_unmgd[reef] <- V_stewart_num[1, ind] /
      (sum(V_stewart_num[, ind],
           na.rm = TRUE)) # scale the first (i.e. 'recov') entry
    
    # Calculate probability of recovered state for area if managed
    # Prob disturbance in area
    d1 <- reef_data$prob_s_dist[reef] %>%
      as.numeric()
    d2 <- reef_data$prob_c_dist[reef] %>%
      as.numeric()
    # Prob Impacted given disturbance
    eta1 <- reef_data$prob_s_impact[reef] %>%
      as.numeric()
    if (is.na(eta1)) {
      sector_s_etas <- as.numeric(reef_data$prob_s_impact[sector_indx])
      eta1 <- sample(sector_s_etas[!is.na(sector_s_etas)], 1) %>%
        as.numeric()
    }
    eta2 <- reef_data$prob_c_impact[reef] %>%
      as.numeric()
    if (is.na(eta2)) {
      sector_c_etas <- as.numeric(reef_data$prob_c_impact[sector_indx])
      eta2 <- sample(sector_c_etas[!is.na(sector_c_etas)], 1) %>%
        as.numeric()
    }
    # Prob recovering post disturbance
    r1 <- reef_data$r_single_mgd[reef] 
    r2 <- reef_data$r_comp_mgd[reef] 
    nrows <- 3
    ncols <- 3
    A_stewart_num <- matrix(
      c(
        (1 - d1 * eta1) * (1 - d2 * eta2), r1, r2,
        d1 * eta1 * (1 - d2 * eta2), 1 - r1, 0,
        d2 * eta2, 0, 1 - r2
      ),
      nrows, ncols,
      byrow = TRUE
    )
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values # eigenvalues
    ind <- which.max(D_stewart_num)
    reef_data$pr_recov_comp_mgd[reef] <- V_stewart_num[1, ind] /
      (sum(V_stewart_num[, ind],
           na.rm = TRUE)) # scale the first (i.e. 'recov') entry
    
    # Calculate probability of recovered state for area if not managed
    r1 <- reef_data$r_single_unmgd[reef] #%>%
    # {
    #   is.na(.) || is.nan(.) || . == "NaN"
    # } %>%
    # ifelse(mean(as.numeric(reef_data$r_single_unmgd), na.rm = TRUE),
    #   reef_data$r_single_unmgd[reef]
    # ) %>%
    # as.numeric()
    r2 <- reef_data$r_comp_unmgd[reef] #%>%
    # {
    #   is.na(.) || is.nan(.) || . == "NaN"
    # } %>%
    # ifelse(mean(as.numeric(reef_data$r_comp_unmgd), na.rm = TRUE),
    #   reef_data$r_comp_unmgd[reef]
    # ) %>%
    # as.numeric()
    A_stewart_num <- matrix(
      c(
        (1 - d1 * eta1) * (1 - d2 * eta2), r1, r2,
        d1 * eta1 * (1 - d2 * eta2), 1 - r1, 0,
        d2 * eta2, 0, 1 - r2
      ),
      nrows, ncols,
      byrow = TRUE
    )
    eigs <- eigen(A_stewart_num)
    V_stewart_num <- eigs$vectors # eigenvectors
    D_stewart_num <- eigs$values # eigenvalues
    ind <- which.max(D_stewart_num)
    reef_data$pr_recov_comp_unmgd[reef] <- V_stewart_num[1, ind] /
      (sum(V_stewart_num[, ind],
           na.rm = TRUE)) # scale the first (i.e. 'recov') entry
  }
  
  return(reef_data)
}
