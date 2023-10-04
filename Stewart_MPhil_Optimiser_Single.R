optimiser_single <- function(df, mgmt_constraint) {
  # Set up single linear programming model
  num_reefs <- nrow(df)
  sing_result <- MIPModel() %>% 
    add_variable(y[i], i = 1:num_reefs, type = "binary") %>% 
    set_objective(sum_expr(y[i] * df$pr_recov_sing_mgd[i] + 
                             (1 - y[i]) * df$pr_recov_sing_unmgd[i],
                           i = 1:num_reefs),
                  sense = "max") %>% 
    add_constraint(sum_expr(y[i], i = 1:num_reefs) <= ceiling(mgmt_constraint * 
                                                                num_reefs)) %>% 
    solve_model(with_ROI("glpk"))
  
  # Return single results
  return(sing_result)
}