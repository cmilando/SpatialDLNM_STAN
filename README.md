# SpatialDLNM_STAN

check out `03_test_SB_CondPois.R` for a 2 case test


# HAHA and then you have to do crosspred by hand lol


# THOUGHTS
-- TO EXTEND TO MULTIPLE REGIONS
-- PASS In a 3-D array for beta [N, K, J] where J is regions

*** the reason INLA takes so long is that there are SO MANY intercepts

****
***** hmm but you'd have to manually code the gnm
*** THATS WHERE THE SPREADSHEET COMES IN BABY

****
***** so you'd need to send in the strata variable somehow, maybe as K
***** might be easier to start with time-series




-- UM JUST LOOK HERE AT PAGE 20
https://cran.r-project.org/web/packages/CARBayes/vignettes/CARBayes.pdf


hmmm bc of dlnm you may want to just do this yourself
and you trust STAN more
and this Laurexz thing seems easy
https://academic.oup.com/ije/article/53/3/dyae061/7654027

inla_formula <- visits_n ~ -1 +
  cb1 + cb2 + cb3 + cb4 + cb5 + cb6 +
  cb7 + cb8 + cb9 +
  f(strata, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) +
  f(id_cb1, cb1, model = "bym", graph = list_neig,
    hyper = list(prec.unstruct = list(prior = sdunif),
                 prec.spatial = list(prior = sdunif)))


# REFERNECES:
* https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/fit_scotland_bym2.R
* https://mc-stan.org/learn-stan/case-studies/icar_stan.html




