

my_output <- readRDS("draws_df_backup.RDS")

head(my_output)
names(my_output)

# beta11 <- draws_df$`mu[1]` + draws_df$`beta_star[1,1]`
cols <- grep("beta_star\\[[0-9]+,1\\]", names(my_output))
names(my_output)[cols]
beta_star <- my_output[, cols]
head(beta_star)

mu_cols <- grep("mu", names(my_output))
names(my_output)[mu_cols]
mu_df <- my_output[, mu_cols]

beta <- mu_df + beta_star
dim(beta)
summary(beta)

# THEIR
load("past_output/final_simsmatrix_model3_spatial_casecrossover.RData")

## so the indices are reversed
i_reg = 1
beta_reg <- winbugs_res[,grepl(paste0("^beta\\[", i_reg,","), 
                               colnames(winbugs_res))]

summary(beta_reg)


#hmm these are pretty different? 
# might be because they parameterized differently, with 
# the sigma and beta star differences ....
# i wonder how this plays out