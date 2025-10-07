#' Bayesian Immune Cell Abundance Model (BICAM)
#'
#' @param dat data frame with dataset (proper setup displayed in tutorial)
#' @param M number of cell types/parameters of interest
#' @param it number of sampling iterations (after burn-in)
#' @param adapt number of adaptation iterations (for compiling model)
#' @param burn number of burn-in iterations
#' @param thin number of thinning samples
#' @param ran_eff indicate whether to use random subject effect (repeated measurements)
#' @param ncov number of covariates input into the model
#' @param chains number of chains to run
#' @param cores number of cores
#' @param v0_mu_logit anticipated proportion of cell types/parameters
#' @param model covariance model selection
#' @param dis distance matrix for Exp. Decay model
#' @param tree tree-structured covariance matrix for Tree and Scaled Tree models
#' @param treelevels list of matrices for multilevel, tree-structured covariance matrix for TreeLevels model
#'
#' @returns A list of inputs and results
#'
#' @import rjags
#'
#' @examples
#' data(dat)
#' BICAM(dat,2,1500,250,250)
#'
#'
#'@export
BICAM <- function(dat,M,adapt,burn,it,
                  thin=1,ran_eff=1,chains=4,cores=4,v0_mu_logit=0.01,ncov=1,
                  model="Unstr",dis=NULL,tree=NULL,treelevels=NULL
){
  #---------- beginning of function ----------
  start_time = Sys.time()


  #---------- Model ----------
  if(!(model == "Unstr" || model == "ExpDecay" || model == "Tree" || model == "TreeLevels" || model == "TreeScaled")){
    stop("Invalid input for 'model' parameter.\nPlease enter model as 'Unstr', 'ExpDecay', 'Tree', or 'TreeScaled'")
  }

  # Check if required matrices are missing
  if (model == "ExpDecay" && is.null(dis)) {
    stop("Distance matrix 'dis' must be provided when model is 'ExpDecay'.")
  }
  if (model == "Tree" && is.null(tree)) {
    stop("Tree-structured covariance matrix 'tree' must be provided when model is 'Tree'.")
  }
  if (model == "TreeLevels" && is.null(treelevels)) {
    stop("List of bases for the tree-structured covariance matrix 'treelevels' must be provided when model is 'TreeLevels'.")
  }
  if (model == "TreeScaled" && is.null(tree)) {
    stop("Tree-structured covariance matrix 'tree' must be provided when model is 'TreeScaled'.")
  }

  ### Initializing potential parameters
  treemat <- 0
  nlevels <- 0
  tree_inv <- diag(2)

  # Check for random subject effect
  if(ran_eff == 1){
    ran_mod <- " + ranS[id[i]]"
    ran_for <- "for(sub in 1:id_unique){
        ranS[sub] ~ dnorm(0,ranS_var[sub])
        ranS_var[sub] ~ dgamma(0.1,0.1)
      }"
  } else {
    ran_mod <- ""
    ran_for <- ""
  }

  if(ncov == 1){
    log_mod <- paste0("bet0[k] + bet1[k] * X[i]",ran_mod)
    X <- dat[,3]
    betas_mod <- "bet0[1:M] ~ dmnorm(v0,bet0_cov_prior)
                  bet1[1:M] ~ dmnorm(v1,bet1_cov_prior)"
    v_mod <- "v0[k] ~ dnorm(v0_mu[k],eta0[k])
              v1[k] ~ dnorm(v1_mu[k],eta1[k])"
    v_mu_mod <- "v0_mu[k] ~ dnorm(v0_mu_init[k],1)
                 v1_mu[k] ~ dnorm(0,1)"
    eta_mod <- "eta0[k] ~ dgamma(2,0.2)
                eta1[k] ~ dgamma(2,0.2)"
  } else {
    bet_mod <- "bet0[k]"
    betas_mod <- "bet0[1:M] ~ dmnorm(v0,bet0_cov_prior)"
    v_mod <- "v0[k] ~ dnorm(v0_mu[k],eta0[k])"
    v_mu_mod <- "v0_mu[k] ~ dnorm(v0_mu_init[k],1)"
    eta_mod <- "eta0[k] ~ dgamma(2,0.2)"
    for(j in 1:ncov){
      bet_mod <- paste0(bet_mod," + bet",j,"[k] * X[i,",j,"]")
      betas_mod <- paste0(betas_mod,"\n",
                          "bet",j,"[1:M] ~ dmnorm(v",j,",bet",j,"_cov_prior)")
      v_mod <- paste0(v_mod,"\n",
                      "v",j,"[k] ~ dnorm(v",j,"_mu[k],eta",j,"[k])")
      v_mu_mod <- paste0(v_mu_mod,"\n",
                         "v",j,"_mu[k] ~ dnorm(0,1)")
      eta_mod <- paste0(eta_mod,"\n",
                        "eta",j,"[k] ~ dgamma(2,0.2)")
    }
    log_mod <- paste0(bet_mod,ran_mod)
    X <- as.matrix(dat[,3:(2+ncov)])
  }


  # Setting up Covariance structure based on model input
  if(model == "Unstr"){
    cov_mod <-
      "bet0_cov_prior ~ dwish(R,df)
      bet1_cov_prior ~ dwish(R,df)"
    if(ncov > 1){
      for(j in 2:ncov){
        cov_mod <- paste0(cov_mod,"\n",
                          "bet",j,"_cov_prior ~ dwish(R,df)")
      }
    }
  } else if(model == "ExpDecay"){
    tau_mod <- "tau0[k] ~ dgamma(0.5,1) T(1e-10,)
                tau1[k] ~ dgamma(0.5,1) T(1e-10,)"
    tau_mat_mod <- "tau0_mat[rowm,colm] <- ifelse(rowm == colm, (1/sqrt(tau0[rowm])),0)
                    tau1_mat[rowm,colm] <- ifelse(rowm == colm, (1/sqrt(tau1[rowm])),0)"
    bet_cor_prior_mod <- "bet0_cor_prior[rowm,colm] <- exp(-(dis[rowm,colm]/zeta0))
                          bet1_cor_prior[rowm,colm] <- exp(-(dis[rowm,colm]/zeta1))"
    zeta_mod <- "zeta0 ~ dgamma(1,0.01) T(1e-10,)
                 zeta1 ~ dgamma(1,0.01) T(1e-10,)"
    bet_cov_prior_mod <- "bet0_cov_prior_inv <- tau0_mat %*% bet0_cor_prior %*% tau0_mat
                          bet0_cov_prior <- inverse(bet0_cov_prior_inv)
                          bet1_cov_prior_inv <- tau1_mat %*% bet1_cor_prior %*% tau1_mat
                          bet1_cov_prior <- inverse(bet1_cov_prior_inv)"
    if(ncov > 1){
      for(j in 2:ncov){
        tau_mod <- paste0(tau_mod,"\n",
                          "tau",j,"[k] ~ dgamma(5,1) T(1e-10,)")
        tau_mat_mod <- paste0(tau_mat_mod,"\n",
                              "tau",j,"_mat[rowm,colm] <- ifelse(rowm == colm, (1/sqrt(tau",j,"[rowm])),0)")
        bet_cor_prior_mod <- paste0(bet_cor_prior_mod,"\n",
                                    "bet",j,"_cor_prior[rowm,colm] <- exp(-(dis[rowm,colm]/zeta",j,"))")
        zeta_mod <- paste0(zeta_mod,"\n",
                           "zeta",j," ~ dgamma(1,0.1) T(1e-10,)")
        bet_cov_prior_mod <- paste0(bet_cov_prior_mod,"\n",
                                    "bet",j,"_cov_prior_inv <- tau",j,"_mat %*% bet",j,"_cor_prior %*% tau",j,"_mat\n",
                                    "bet",j,"_cov_prior <- inverse(bet",j,"_cov_prior_inv)")
      }
    }

    cov_mod <- paste0("
    for (k in 1:M) {",
                      tau_mod,
                      "}
    for(rowm in 1:M){
        for(colm in 1:M){",
                      tau_mat_mod,
                      bet_cor_prior_mod,
                      "
        }
      }",

                      zeta_mod,
                      bet_cov_prior_mod)
  } else if(model == "Tree"){
    tree_inv <- matlib::inv(tree)
    cov_mod <-
      paste("bet0_cov_prior = tree_inv",
            "bet1_cov_prior = tree_inv",sep="\n")
    if(ncov > 1){
      for(j in 2:ncov){
        cov_mod <- paste0(cov_mod,"\n",
                          "bet",j,"_cov_prior = tree_inv")
      }
    }
  }  else if(model == "TreeScaled"){
    tree_inv <- matlib::inv(tree)
    cov_mod <-
      "
      lambda0 ~ dgamma(0.5,0.1) T(0.00001,)
      lambda1 ~ dgamma(0.5,0.1) T(0.00001,)

      bet0_cov_prior <- (1/lambda0) * tree_inv
      bet1_cov_prior <- (1/lambda1) * tree_inv"
    if(ncov > 1){
      for(j in 2:ncov){
        cov_mod <- paste0(cov_mod,"\n",
                          "lambda",j," ~ dgamma(0.5,0.1) T(0.00001,)","\n",
                          "bet",j,"_cov_prior <- (1/lambda",j,") * tree_inv")
      }
    }
  } else if(model == "TreeLevels"){
    nlevels <- length(treelevels)
    treemat <- array(unlist(treelevels), dim = c(M, M, nlevels))
    omega_mod <- "omega0[l] ~ dgamma(0.5,0.1) T(0.0001,)
                  omega1[l] ~ dgamma(0.5,0.1) T(0.0001,)"
    base_mod <- "base0[lrow,lcol,l] <- omega0[l]*treemat[lrow,lcol,l]
                 base1[lrow,lcol,l] <- omega1[l]*treemat[lrow,lcol,l]"
    tree_sum_temp_mod <- "tree_sum0_temp[lrow,lcol,1] <- base0[lrow,lcol,1]
                          tree_sum1_temp[lrow,lcol,1] <- base1[lrow,lcol,1]"
    tree_sum_temp_mod2 <- "tree_sum0_temp[lrow,lcol,l] <- tree_sum0_temp[lrow,lcol,l-1] + base0[lrow,lcol,l]
                           tree_sum1_temp[lrow,lcol,l] <- tree_sum1_temp[lrow,lcol,l-1] + base1[lrow,lcol,l]"
    tree_sum_mod <- "tree_sum0[lrow,lcol] <- tree_sum0_temp[lrow,lcol,nlevels]
                     tree_sum1[lrow,lcol] <- tree_sum1_temp[lrow,lcol,nlevels]"
    bet_cov_prior_mod <- "bet0_cov_prior <- inverse(tree_sum0)
                          bet1_cov_prior <- inverse(tree_sum1)"
    if(ncov > 1){
      for(j in 2:ncov){
        omega_mod <- paste0(omega_mod,"\n",
                            "omega",j,"[l] ~ dgamma(0.5,0.1) T(0.0001,)")
        base_mod <- paste0(base_mod,"\n",
                           "base",j,"[lrow,lcol,l] <- omega",j,"[l]*treemat[lrow,lcol,l]")
        tree_sum_temp_mod <- paste0(tree_sum_temp_mod,"\n",
                                    "tree_sum",j,"_temp[lrow,lcol,1] <- base",j,"[lrow,lcol,1]")
        tree_sum_temp_mod2 <- paste0(tree_sum_temp_mod2,"\n",
                                     "tree_sum",j,"_temp[lrow,lcol,l] <- tree_sum",j,"_temp[lrow,lcol,l-1] + base",j,"[lrow,lcol,l]")
        tree_sum_mod <- paste0(tree_sum_mod,"\n",
                               "tree_sum",j,"[lrow,lcol] <- tree_sum",j,"_temp[lrow,lcol,nlevels]")
        bet_cov_prior_mod <- paste0(bet_cov_prior_mod,"\n",
                                    "bet",j,"_cov_prior <- inverse(tree_sum",j,")")
      }
    }
    cov_mod <- paste0("
    for(l in 1:nlevels) {",
                      omega_mod,
                      "for(lrow in 1:M){
          for(lcol in 1:M){",
                      base_mod,
                      "}
        }
    }

     for(lrow in 1:M){
      for(lcol in 1:M){",
                      tree_sum_temp_mod,
                      "for(l in 2:nlevels) {",
                      tree_sum_temp_mod2,
                      "}",
                      tree_sum_mod,
                      "}
    }",
                      bet_cov_prior_mod)
  }


  # Constructing full model string for JAGS
  model_string <- paste("
    model {
      for (k in 1:M) {
        for(i in 1:N){
          y[i,k] ~ dbinom(p_temp[i,k], tot_cells[i])
          mu_temp[i,k] <- ilogit(",log_mod,")
          p_temp[i,k] ~ dbeta(alp[i,k],bet[i,k]) T(1e-10,(1-1e-10))
          alp[i,k] <- mu_temp[i,k]*gam[k]
          bet[i,k] <- (1-mu_temp[i,k])*gam[k]
        }
        mu[k] <- mean(mu_temp[,k])
        p[k] <- mean(p_temp[,k])
        gam[k] ~ dgamma(7,0.1)",
                        eta_mod,
                        v_mu_mod,
                        v_mod,
                        "}",
                        ran_for,
                        cov_mod,
                        betas_mod,
                        "
    }
  ",sep="\n")


  #---------- Data Prepping ----------
  dat[,1] <- as.numeric(as.factor(dat[,1]))
  I_M <- diag(M)
  Y <- as.matrix(cbind(dat[,(3+ncov):(2+ncov+M)]))
  R = diag(M)
  df=M+1
  v0_mu_init = VGAM::logitlink(v0_mu_logit)
  data1 <- list(tot_cells=dat[,2],
                N=dim(dat)[1],
                y=Y,
                X=X,
                M=M,
                id=dat[,1],
                id_unique=length(unique(dat[,1])),
                v0_mu_init=rep(v0_mu_init,M),
                dis=dis,
                tree_inv=tree_inv,
                treemat=treemat,
                nlevels=nlevels,
                R=R,
                df=df
  )

  # Removing items from data list
  if(ran_eff != 1){
    data1 <- data1[names(data1) != "id"]
    data1 <- data1[names(data1) != "id_unique"]
  }
  if(model != "Unstr"){
    data1 <- data1[names(data1) != "R"]
    data1 <- data1[names(data1) != "df"]
  }
  if(model != "ExpDecay"){
    data1 <- data1[names(data1) != "dis"]
  }
  if(!(model == "Tree" || model == "TreeScaled")){
    data1 <- data1[names(data1) != "tree_inv"]
  }
  if(model != "TreeLevels"){
    data1 <- data1[names(data1) != "treemat"]
    data1 <- data1[names(data1) != "nlevels"]
  }


  #---------- Variable Monitoring ----------
  if(model == "Unstr"){
    var_names <- c("p","gam","bet0","bet1","mu","v0","v1",
                   "eta0","eta1","v0_mu","v1_mu",
                   "bet0_cov_prior","bet1_cov_prior"
    )
    if(ncov > 1){
      for(j in 2:ncov){
        var_names <- c(var_names,paste0("bet",j),paste0("v",j),paste0("eta",j),
                       paste0("v",j,"_mu"),paste0("bet",j,"_cov_prior"))
      }
    }
  } else if(model == "ExpDecay"){
    var_names <- c("p","gam","bet0","bet1","mu","v0","v1","v0_mu","v1_mu",
                   "eta0","eta1","tau0","tau1","zeta0","zeta1",
                   "bet0_cor_prior","bet1_cor_prior",
                   "bet0_cov_prior","bet1_cov_prior"
    )
    if(ncov > 1){
      for(j in 2:ncov){
        var_names <- c(var_names,paste0("bet",j),paste0("v",j),paste0("eta",j),
                       paste0("v",j,"_mu"),paste0("bet",j,"_cov_prior"),
                       paste0("tau",j),paste0("zeta",j),paste0("bet",j,"_cor_prior"))
      }
    }
  } else if(model == "Tree"){
    var_names <- c("p","gam","bet0","bet1","mu","v0","v1",
                   "eta0","eta1","v0_mu","v1_mu"
    )
    if(ncov > 1){
      for(j in 2:ncov){
        var_names <- c(var_names,paste0("bet",j),paste0("v",j),paste0("eta",j),
                       paste0("v",j,"_mu"))
      }
    }
  } else if(model == "TreeScaled"){
    var_names <- c("p","gam","bet0","bet1","mu","v0","v1",
                   "eta0","eta1","v0_mu","v1_mu",
                   "lambda0","lambda1"
    )
    if(ncov > 1){
      for(j in 2:ncov){
        var_names <- c(var_names,paste0("bet",j),paste0("v",j),paste0("eta",j),
                       paste0("v",j,"_mu"),paste0("lambda",j))
      }
    }
  } else if(model == "TreeLevels"){
    var_names <- c("p","gam","bet0","bet1","mu","v0","v1",
                   "eta0","eta1","v0_mu","v1_mu",
                   "tree_sum0","tree_sum1","omega0","omega1"
    )
    if(ncov > 1){
      for(j in 2:ncov){
        var_names <- c(var_names,paste0("bet",j),paste0("v",j),paste0("eta",j),
                       paste0("v",j,"_mu"),paste0("tree_sum",j),paste0("omega",j))
      }
    }
  }


  #---------- Initial Values for Sampling ----------
  RNG_name_list <- rep(c("base::Wichmann-Hill","base::Marsaglia-Multicarry","base::Super-Duper",
                         "base::Mersenne-Twister","lecuyer::RngStream","user::baseRNG"),times=ceiling(chains/6))
  inits <- vector(mode='list', length=chains)
  for(j in 1:chains){
    inits_temp <-
      list(p_temp = matrix(rep(v0_mu_logit,M*dim(dat)[1]),nrow=dim(dat)[1]), gam = rep(50,M),
           bet0 = rep(v0_mu_init,M), bet1 = rep(0,M),
           v0 = rep(0,M), v1 = rep(0,M), .RNG.name = RNG_name_list[j], .RNG.seed = j,
           eta0 = rep(1,M), eta1 = rep(1,M),
           v0_mu = rep(v0_mu_init,M),v1_mu = rep(0,M))
    if(ncov > 1){
      for(l in 2:ncov){
        inits_temp[[paste0("bet",l)]] <- rep(0,M)
        inits_temp[[paste0("v",l)]] <- rep(0,M)
        inits_temp[[paste0("eta",l)]] <- rep(1,M)
        inits_temp[[paste0("v",l,"_mu")]] <- rep(0,M)
      }
    }
    if(ran_eff == 1){
      inits_temp[["ranS"]] = rep(0,length(unique(dat[,1])))
      inits_temp[["ranS_var"]] = rep(1,length(unique(dat[,1])))
    }
    if(model == "Unstr"){
      inits_temp[["bet0_cov_prior"]] <- diag(M)
      inits_temp[["bet1_cov_prior"]] <- diag(M)
      if(ncov > 1){
        for(jcov in 2:ncov){
          inits_temp[[paste0("bet",jcov,"_cov_prior")]] <- diag(M)
        }
      }
    }
    if(model == "ExpDecay"){
      inits_temp[["zeta0"]] <- 1
      inits_temp[["zeta1"]] <- 1
      if(ncov > 1){
        for(jcov in 2:ncov){
          inits_temp[[paste0("zeta",jcov)]] <- 1
        }
      }
    }
    if(model == "TreeScaled"){
      inits_temp[["lambda0"]] <- 1
      inits_temp[["lambda1"]] <- 1
      if(ncov > 1){
        for(jcov in 2:ncov){
          inits_temp[[paste0("lambda",jcov)]] <- 1
        }
      }
    }
    if(model == "TreeLevels"){
      inits_temp[["omega0"]] <- rep(1,nlevels)
      inits_temp[["omega1"]] <- rep(1,nlevels)
      if(ncov > 1){
        for(jcov in 2:ncov){
          inits_temp[[paste0("omega",jcov)]] <- rep(1,nlevels)
        }
      }
    }
    inits[[j]] <- inits_temp
  }


  #---------- Compiling/Running Model ----------

  message("\nSampling From Posterior\n\n")

  samples_temp <- runjags::run.jags(model = model_string,
                                    data = data1,
                                    n.chains = chains,
                                    inits = inits,
                                    adapt = adapt,
                                    burnin = burn,
                                    thin = thin,
                                    sample = it,
                                    monitor = var_names,
                                    method = "parallel")

  message("\nSampling Finished\n\n")

  samples <- coda::as.mcmc.list(samples_temp)
  end_time = Sys.time()
  my_time = end_time-start_time

  BHBBM_results <- list(samples,my_time,model_string,inits,data1,var_names)
}
