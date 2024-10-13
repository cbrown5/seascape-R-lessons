#Code from https://github.com/emitanaka/boral/blob/master/R/makejagsboralmodel.R to make jags code for a boral model

make.jagsboralmodel <- function(family, num.X = 0, X.ind = NULL, num.traits = 0, which.traits = NULL,
                                lv.control = list(num.lv = 2, type = "independent"), row.eff = "none", row.ids = NULL,
                                offset = NULL, trial.size = 1, n, p, model.name = NULL,
                                prior.control = list(type = c("normal", "normal", "normal", "uniform"), hypparams = c(10, 10, 10, 30), ssvs.index = -1, ssvs.g = 1e-6, ssvs.traitsindex = -1), num.lv = NULL) {
  check_which_traits(num.traits = num.traits, which.traits = which.traits, num.X = num.X, makejagsboralfile_messages = TRUE)
  if (is.null(which.traits)) {
    which.traits <- vector("list", num.X + 1)
    for (k in 1:(num.X + 1)) {
      which.traits[[k]] <- 0
    }
  }
  
  lv.control <- check_lv_control(num.lv = num.lv, lv.control = lv.control, need.distmat = FALSE)
  num.lv <- lv.control$num.lv
  
  complete_trial_size <- check_trial_size(family = family, trial.size = trial.size, create.complete.trial.size = TRUE, y = matrix(NA, nrow = 1, ncol = p))
  
  complete_family <- check_family(family = family, y = matrix(NA, nrow = 1, ncol = p), traits = NULL)
  
  if (row.eff != "none" && is.null(row.ids)) {
    row.ids <- matrix(1:n, ncol = 1)
    message("row.ids assumed to be matrix with one column and elements 1,2,...n i.e., a row-specific intercept")
  }
  if (!is.null(row.ids)) {
    row.ids <- as.matrix(row.ids)
    if (nrow(row.ids) != n) {
      stop("Number of rows in the matrix row.ids should be equal to n")
    }
    if (is.null(colnames(row.ids))) {
      colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
    }
  }
  
  
  if (!is.null(offset)) {
    if (!is.matrix(offset)) {
      stop("offset could be a matrix with the same dimensions as y")
    }
  }
  
  
  prior.control <- fillin_prior_control(x = prior.control)
  check_prior_control(prior.control = prior.control)
  
  if (length(prior.control$ssvs.index) == 1 & num.X > 0) {
    prior.control$ssvs.index <- rep(prior.control$ssvs.index, num.X)
  }
  if (num.traits > 0) {
    if (!is.list(prior.control$ssvs.traitsindex)) {
      prior.control$ssvs.traitsindex <- vector("list", num.X + 1)
      for (k in 1:(num.X + 1)) {
        prior.control$ssvs.traitsindex[[k]] <- rep(-1, length(which.traits[[k]]))
      }
    }
    if (is.list(prior.control$ssvs.traitsindex)) {
      check_ssvstraits(prior.control$ssvs.traitsindex, which.traits)
    }
  }
  
  X.ind <- check_X_ind(X.ind = X.ind, p = p, num.X = num.X, prior.control = prior.control)
  
  index.ord.cols <- which(complete_family == "ordinal")
  index.tweed.cols <- which(complete_family == "tweedie")
  
  
  ## Checks done; starting writing JAGS script!
  
  
  model_script <- paste0("## JAGS model written for boral version ", packageDescription("boral")$Version, " on ", as.character(Sys.time()), " ##\n\n model {")
  model_script <- c(model_script, "\t ## Data Level ## \n\t for(i in 1:n) {")
  
  write.resp.script <- setup.resp.families.lv(p = p, complete.family = complete_family, num.lv = num.lv, row.eff = row.eff, row.ids = row.ids, offset = offset, num.X = num.X, complete.trial.size = complete_trial_size, index.tweed.cols = index.tweed.cols, index.ord.cols = index.ord.cols)
  model_script <- c(model_script, write.resp.script)
  model_script <- c(model_script, paste0("\t\t }"))
  
  model_script <- c(model_script, paste0("\t ## Latent variables ##"))
  if (lv.control$type == "independent") {
    model_script <- c(model_script, paste0("\t for(i in 1:n) { for(k in 1:num.lv) { lvs[i,k] ~ dnorm(0,1) } } \n\n\t ## Process level and priors ##"))
  }
  if (lv.control$type == "exponential") {
    model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- exp(-distmat[k1,k2]/lv.covparams[1]) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))
  }
  if (lv.control$type == "squared.exponential") {
    model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- exp(-pow(distmat[k1,k2]/lv.covparams[1],2)) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))
  }
  if (lv.control$type == "powered.exponential") {
    model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- exp(-pow(distmat[k1,k2]/lv.covparams[1],lv.covparams[2])) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))
  }
  if (lv.control$type == "spherical") {
    model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- step(lv.covparams[1] - distmat[k1,k2])*(1 - 1.5*distmat[k1,k2]/lv.covparams[1] + 0.5*pow(distmat[k1,k2]/lv.covparams[1],3)) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))
  }
  # if(lv.control$type == "cauchy") ## DOES NOT WORK VERY WELL!
  # model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- pow(1+ pow(distmat[k1,k2]/lv.covparams[1],2), -lv.covparams[2]) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))
  ## Matern not implemented to due complications/lack of direct availability of a BesselK function
  rm(write.resp.script)
  
  ## Build prior strings for all priors distributions
  prior.strings <- construct.prior.strings(x = prior.control)
  
  
  ## Code for column-specific intercept. Note this is set up different to how X variables are set up to save some coding space!
  ## No traits or traits included but not regressed against intercept
  if (num.traits == 0 || (num.traits > 0 & which.traits[[1]][1] == 0)) {
    ## Not ordinal columns, then as per usual
    if (length(index.ord.cols) == 0) {
      model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ ", prior.strings$p1, " } ## Separate species intercepts"))
    }
    ## If 1 ordinal column, then intercept for this column equal 0
    if (length(index.ord.cols) == 1) {
      model_script <- c(model_script, paste0("\t lv.coefs[", index.ord.cols, ",1] <- 0 ## Single ordinal species intercept"))
      for (j in (1:p)[-index.ord.cols]) {
        model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ ", prior.strings$p1, "All other species intercepts"))
      }
    }
    ## More than 1 ordinal column, then set up random intercept for this species
    if (length(index.ord.cols) > 1) {
      if (length(index.ord.cols) == p) {
        model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ dnorm(0,pow(ordinal.sigma,-2)) } ## Random intercept for all ordinal species"))
      } else {
        for (j in index.ord.cols) {
          model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ dnorm(0,pow(ordinal.sigma,-2)) ## Random intercept for all ordinal species"))
        }
        for (j in (1:p)[-index.ord.cols]) {
          model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ ", prior.strings$p1, "All other species intercepts"))
        }
      }
      model_script <- c(model_script, paste0("\t ordinal.sigma ~ ", prior.strings$p4))
    }
    if ((num.traits > 0 & which.traits[[1]][1] == 0)) {
      model_script <- c(model_script, paste0("\t traits.int[1] <- 0; for(l in 1:num.traits) { traits.coefs[1,l] <- 0 } \n\t trait.sigma[1] <- 0 ## Traits not used for intercept"))
    }
  }
  ## Traits included in model and regressed against intercept
  if (num.traits > 0 & all(which.traits[[1]] > 0)) {
    ## If there are 0 or > 1 ordinal columns, then regress all intercepts against traits
    if (length(index.ord.cols) != 1) {
      model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ dnorm(traits.int[1] + inprod(traits[j,],traits.coefs[1,1:num.traits]),pow(trait.sigma[1],-2)) } ## Species intercepts regressed against traits"))
    }
    ## If there is 1 ordinal column, do not regress this intercept against trait
    if (length(index.ord.cols) == 1) {
      model_script <- c(model_script, paste0("\t lv.coefs[", index.ord.cols, ",1] <- 0 ## Ordinal species intercept"))
      for (j in (1:p)[-index.ord.cols]) {
        model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ dnorm(traits.int[1] + inprod(traits[", j, ",],traits.coefs[1,1:num.traits]),pow(trait.sigma[1],-2)) ## All other intercepts"))
      }
    }
    model_script <- c(model_script, paste0("\t traits.int[1] ~ ", prior.strings$p1))
    for (l in which.traits[[1]]) {
      if (prior.control$ssvs.traitsindex[[1]][which(which.traits[[1]] == l)] == -1) {
        model_script <- c(model_script, paste0("\t traits.coefs[", 1, ",", l, "] ~ ", prior.strings$p1, " ## Traits used for intercept"))
      }
      if (prior.control$ssvs.traitsindex[[1]][which(which.traits[[1]] == l)] == 0) {
        ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[1], "*((1-ssvs.traitscoefs1", l, ")*", prior.control$ssvs.g, " + ssvs.traitscoefs1", l, "),-1)); ssvs.traitscoefs1", l, " ~ dbern(0.5)")
        model_script <- c(model_script, paste0("\t traits.coefs[", 1, ",", l, "] ~ ", ssvs.prior.string, " ## Traits used for intercept"))
      }
    }
    if (length((1:num.traits)[-which.traits[[1]]]) > 0) {
      for (l in (1:num.traits)[-which.traits[[1]]]) {
        model_script <- c(model_script, paste0("\t traits.coefs[", 1, ",", l, "] <- 0 ## Traits not used for intercept"))
      }
    }
    model_script <- c(model_script, paste0("\t trait.sigma[1] ~ ", prior.strings$p4))
  }
  
  
  if (any(complete_family == "tweedie")) {
    model_script <- c(model_script, paste0("\t powerparam ~ dunif(1,2) ## Tweedie power parameter"))
  }
  if (any(complete_family == "ordinal")) {
    model_script <- c(model_script, paste0("\t for(k in 1:(num.ord.levels-1)) { cutoffs0[k] ~ ", prior.strings$p1, " }"))
    model_script <- c(model_script, paste0("\t cutoffs[1:(num.ord.levels-1)] <- sort(cutoffs0) ## Ordinal cutoffs"))
  }
  
  
  ## Priors on row effects
  if (row.eff == "fixed") {
    for (k in 1:ncol(row.ids))
    {
      model_script <- c(model_script, paste0("\n\t row.coefs.ID", k, "[1] <- 0"))
      model_script <- c(model_script, paste0("\n\t for(i in 2:n.ID[", k, "]) { row.coefs.ID", k, "[i] ~ ", prior.strings$p1, " } "))
    }
  }
  if (row.eff == "random") {
    for (k in 1:ncol(row.ids)) {
      model_script <- c(model_script, paste0("\n\t for(i in 1:n.ID[", k, "]) { row.coefs.ID", k, "[i] ~ dnorm(0, pow(row.sigma.ID", k, ",-2)) } "))
      model_script <- c(model_script, paste0("\t row.sigma.ID", k, " ~ ", prior.strings$p4))
    }
  }
  
  
  ## Priors on latent variables if required, controlled by prior.control$hypparams[2]
  if (lv.control$type %in% c("exponential", "squared.exponential", "spherical")) {
    model_script <- c(model_script, paste0("\t lv.covparams[1] ~ ", prior.strings$p22))
  }
  if (lv.control$type %in% c("powered.exponential")) {
    model_script <- c(model_script, paste0("\t lv.covparams[1] ~ ", prior.strings$p22))
    model_script <- c(model_script, paste0("\t lv.covparams[2] ~ dunif(0,2)"))
  }
  
  ## Priors on Latent variable coefficients, controlled by prior.control$hypparams[2]
  model_script <- c(model_script, paste0("\n\t for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { lv.coefs[i,j] <- 0 } } ## Constraints to 0 on upper diagonal"))
  model_script <- c(model_script, paste0("\t for(i in 1:num.lv) { lv.coefs[i,i+1] ~ ", prior.strings$p22, " } ## Sign constraints on diagonal elements"))
  model_script <- c(model_script, paste0("\t for(i in 2:num.lv) { for(j in 2:i) { lv.coefs[i,j] ~ ", prior.strings$p2, " } } ## Free lower diagonals"))
  model_script <- c(model_script, paste0("\t for(i in (num.lv+1):p) { for(j in 2:(num.lv+1)) { lv.coefs[i,j] ~ ", prior.strings$p2, " } } ## All other elements"))
  
  
  ## Prior for X-coefficients, controlled by prior.control$hypparams[3]
  if (num.X > 0) {
    model_script <- c(model_script, paste0("\n"))
    ## Traits not included in model
    if (num.traits == 0) {
      for (i in 1:length(prior.control$ssvs.index)) {
        if (prior.control$ssvs.index[i] == -1) {
          if (!is.null(X.ind)) {
            model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", prior.strings$p3, "I(-X.ind[j,", i, "],X.ind[j,", i, "]) } "))
          }
          if (is.null(X.ind)) {
            model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", prior.strings$p3, " } "))
          }
        }
        if (prior.control$ssvs.index[i] == 0) {
          ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[3], "*((1-ssvs.indX", i, "[j])*", prior.control$ssvs.g, " + ssvs.indX", i, "[j]),-1)); ssvs.indX", i, "[j] ~ dbern(0.5)")
          model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", ssvs.prior.string, " }"))
        }
        if (prior.control$ssvs.index[i] > 0) {
          ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[3], "*((1-ssvs.gp", prior.control$ssvs.index[i], ")*", prior.control$ssvs.g, " + ssvs.gp", prior.control$ssvs.index[i], "),-1))")
          model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", ssvs.prior.string, " } "))
        }
      }
    }
    
    if (num.traits > 0) {
      for (i in 1:num.X) {
        ## Traits included but X coefs not regressed against them
        if (which.traits[[i + 1]][1] == 0) {
          model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", prior.strings$p3, " } ## Coefficient not regressed against any traits"))
          model_script <- c(model_script, paste0("\t traits.int[", i + 1, "] <- 0; trait.sigma[", i + 1, "] <- 0; for(l in 1:num.traits) { traits.coefs[", i + 1, ",l] <- 0 } \n"))
        }
        ## Traits included and X coefs regressed against some of them
        if (all(which.traits[[i + 1]] > 0)) {
          model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ dnorm(traits.int[", i + 1, "] + inprod(traits[j,],traits.coefs[", i + 1, ",1:num.traits]),pow(trait.sigma[", i + 1, "],-2)) } "))
          model_script <- c(model_script, paste0("\t traits.int[", i + 1, "] ~ ", prior.strings$p3))
          for (l in which.traits[[i + 1]]) {
            if (prior.control$ssvs.traitsindex[[i + 1]][which(which.traits[[i + 1]] == l)] == -1) {
              model_script <- c(model_script, paste0("\t traits.coefs[", i + 1, ",", l, "] ~ ", prior.strings$p3, " ## Traits used for this X.coefs"))
            }
            if (prior.control$ssvs.traitsindex[[i + 1]][which(which.traits[[i + 1]] == l)] == 0) {
              ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[3], "*((1-ssvs.traitscoefs", i + 1, l, ")*", prior.control$ssvs.g, " + ssvs.traitscoefs", i + 1, l, "),-1)); ssvs.traitscoefs", i + 1, l, " ~ dbern(0.5)")
              model_script <- c(model_script, paste0("\t traits.coefs[", i + 1, ",", l, "] ~ ", ssvs.prior.string, " ## Traits used for this X.coefs"))
            }
          }
          if (length((1:num.traits)[-which.traits[[i + 1]]]) > 0) {
            for (l in (1:num.traits)[-which.traits[[i + 1]]]) {
              model_script <- c(model_script, paste0("\t traits.coefs[", i + 1, ",", l, "] <- 0 ## traits not used for this X.coefs"))
            }
          }
          model_script <- c(model_script, paste0("\t trait.sigma[", i + 1, "] ~ ", prior.strings$p4, "\n"))
        }
      }
    }
    
    model_script <- c(model_script, paste0(""))
    if (any(prior.control$ssvs.index > 0)) {
      for (i in unique(prior.control$ssvs.index[prior.control$ssvs.index > 0])) {
        model_script <- c(model_script, paste0("\t ssvs.gp", i, " ~ dbern(0.5)"))
      }
    }
  }
  # 	if(num.X > 0 & any(family == "multinom")) {
  # 		model_script <- c(model_script, paste0("\t for(j in 1:",length(index_multinom_cols),") { for(i in 1:num.X) { X.multinom.params[j,i,1] <- 0 } } "))
  # 		model_script <- c(model_script, paste0("\t for(j in 1:",length(index_multinom_cols),") { for(i in 1:num.X) { for(k in 2:num.multinom.levels) { X.multinom.params[j,i,k] ~ dnorm(0,",1/prior.control$hypparams[1],") } } } "))
  # 		}
  
  
  ## Prior on dispersion parameters, controlled by prior.control$hypparams[4]
  if (!all(complete_family %in% c("poisson", "binomial", "ordinal", "multinom", "exponential"))) { # "bernoulli"
    model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,num.lv+2] ~ ", prior.strings$p4, " } ## Dispersion parameters"))
  }
  
  model_script <- c(model_script, "\n\t }")
  
  
  if (!is.null(model.name)) {
    write(model_script, file = model.name)
  }
  if (is.null(model.name)) {
    write(model_script, file = "jagsboralmodel.txt")
  }
}

###############################
## Unseen check functions
###############################

check_domarglik_ics <- function(fit.mcmc.names, index.ordinal.cols) {
  out <- TRUE
  
  if (length(grep("traits.coefs", fit.mcmc.names)) > 1) {
    out <- FALSE
  }
  if (length(index.ordinal.cols) > 1) {
    out <- FALSE
  }
  if (length(grep("lv.covparams", fit.mcmc.names)) > 1) {
    out <- FALSE
  }
  
  return(out)
}


check_family <- function(family, y, traits = NULL) {
  if (length(family) != ncol(y) & length(family) != 1) {
    stop("Number of elements in family must either one or the # of columns in y.")
  }
  if (length(family) == 1) {
    complete_family <- rep(family, ncol(y))
  }
  if (length(family) == ncol(y)) {
    complete_family <- family
  }
  complete_family <- match.arg(complete_family, choices = c(
    "negative.binomial", "poisson", "binomial",
    "normal", "lnormal", "tweedie", "ordinal", "exponential", "gamma", "beta"
  ), several.ok = TRUE)
  if (length(complete_family) != ncol(y)) {
    stop("At least one of the elements in family is not supported in current version of boral...sorry!")
  }
  
  if (any(complete_family == "ordinal")) {
    if (sum(y[, complete_family == "ordinal"] == 0) > 0) {
      stop("For ordinal data, please shift minimum level to 1")
    }
    if (!is.null(traits) & (sum(complete_family == "ordinal") == 1)) {
      message("The intercept for the single ordinal response is set to zero and not regressed traits for parameter identifiability reasons")
    }
  }
  
  return(complete_family)
  #     if(all(complete_family == "binomial") & all(complete_trial_size == 1))
  #         family <- rep("bernoulli",p)
}


check_lv_control <- function(num.lv, lv.control, need.distmat = TRUE) {
  if (is.null(lv.control$type)) {
    lv.control$type <- "independent"
  }
  lv.control$type <- match.arg(lv.control$type, choices = c("independent", "exponential", "squared.exponential", "powered.exponential", "spherical"))
  
  if (!is.null(num.lv)) {
    warning("num.lv is now a redundant argument and replaced by lv.control$num.lv. Please set num.lv = NULL.", immediate. = TRUE)
    # lv.control$num.lv <- num.lv
  }
  if (is.null(lv.control$num.lv)) {
    lv.control$num.lv <- 0
  }
  
  if (need.distmat) {
    if (lv.control$type != "independent" & is.null(lv.control$distmat)) {
      stop("If structured latent variables are used, then please supply a distance matrix to lv.control$distmat.")
    }
  }
  
  if (lv.control$num.lv > 5) {
    warning("We won't stop you, but please consider if you really want more than five latent variables in the model!", immediate. = TRUE)
  }
  
  return(lv.control)
}


check_offset <- function(offset = NULL, y) {
  if (!is.null(offset)) {
    if (!is.matrix(offset)) {
      stop("offset should be a matrix with the same dimensions as y")
    }
    if (nrow(offset) != nrow(y)) {
      stop("offset should be a matrix with the same dimensions as y")
    }
    if (ncol(offset) != ncol(y)) {
      stop("offset should be a matrix with the same dimensions as y")
    }
  }
}


check_prior_control <- function(prior.control) {
  if (length(prior.control$hypparams) != 4 || length(prior.control$type) != 4) {
    stop("prior.control$type and prior.control$hypparams must be a vector of four elements. Please see boral help file as to what the elements correspond to.\n")
  }
  if (!all(prior.control$type[-4] %in% c("normal", "uniform", "cauchy"))) {
    stop("At least one of the first three elements of prior.control$type is not supported in the current version of boral...sorry!")
  }
  if (!(prior.control$type[4] %in% c("uniform", "halfcauchy", "halfnormal"))) {
    stop("The fourth element of prior.control$type is not supported in the current version of boral...sorry!")
  }
}


check_row_ids <- function(row.ids = NULL, y) {
  if (!is.null(row.ids)) {
    row.ids <- as.matrix(row.ids)
    if (nrow(row.ids) != nrow(y)) {
      stop("Number of rows in the matrix row.ids should be equal to number of rows in y")
    }
    if (is.null(colnames(row.ids))) {
      colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
    }
  }
  
  return(row.ids)
}


check_row_params <- function(row.params = NULL, y, row.ids = NULL) {
  if (!is.null(row.params)) {
    if (is.null(row.ids)) {
      row.ids <- matrix(1:nrow(y), ncol = 1)
      colnames(row.ids) <- "ID1"
    }
    if (!is.list(row.params)) {
      stop("row.params should be a list with length equal to the number of columns in row.ids")
    }
    if (length(row.params) != ncol(row.ids)) {
      stop("row.params should be a list with length equal to the number of columns in row.ids")
    }
  }
}

check_ssvstraits <- function(ssvs.traitsindex, which.traits) {
  if (!is.null(which.traits)) {
    if (length(ssvs.traitsindex) != length(which.traits)) {
      stop("Both prior.control$ssvs.traitsindex and which.traits should be lists of equal length")
    }
    if (!all(unlist(which.traits) >= 0)) {
      stop("All elements of which.traits must be non-negative")
    }
    if (!all(unlist(ssvs.traitsindex) %in% c(-1, 0))) {
      stop("All elements in the list prior.control$ssvs.traitsindex should be either equal to -1 (no SSVS) or 0 (SSVS applied)")
    }
    for (k in 1:length(which.traits))
    {
      if (which.traits[[k]][1] == 0 & !all(ssvs.traitsindex[[k]] == -1)) {
        stop(paste0("If which.traits[[", k, "]][1] == 0, then all the elements of ssvs.traitsindex[[", k, "]] must equal to -1. That is, if traits are not used then for covariate ", k, " and no SSVS can be done on this"))
      }
      if (any(which.traits[[k]] > 0) & length(ssvs.traitsindex[[k]]) != length(which.traits[[k]])) {
        stop(paste0("If the elements of which.traits[[", k, "]] are positive, then the length of ssvs.traitsindex[[", k, "]] must match the length of which.traits[[", k, "]]. That is, if traits are used then for covariate ", k, " then a corresponding index needs to be supplied to determine if SSVS needs to be done for this"))
      }
    }
  }
}


check_traits <- function(traits, y) {
  if (!is.null(traits)) {
    if (!is.matrix(traits)) {
      traits <- as.matrix(traits)
    }
    if (nrow(traits) != ncol(y)) {
      stop("If traits are supplied, then please ensure the number of rows in traits i.e., number of species, is equal to the number of columns in the response matrix")
    }
    if (any(apply(traits, 2, function(x) all(x == 1)))) {
      stop("No intercept column should be included in traits. It will be included automatically.")
    }
  }
}


check_trial_size <- function(family, trial.size, create.complete.trial.size = FALSE, y = NULL) {
  if (any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) {
    stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y..")
  }
  
  if (create.complete.trial.size) {
    if (any(family == "binomial") & length(trial.size) == 1) {
      complete_trial_size <- rep(0, ncol(y))
      complete_trial_size[which(family == "binomial")] <- trial.size
    }
    if (any(family == "binomial") & length(trial.size) == ncol(y)) {
      complete_trial_size <- trial.size
    }
    if (all(family != "binomial")) {
      complete_trial_size <- rep(0, ncol(y))
    }
    return(complete_trial_size)
  }
}


check_X_ind <- function(X.ind = NULL, p, num.X, prior.control) {
  if (!is.null(X.ind)) {
    X.ind <- as.matrix(X.ind)
    if (nrow(X.ind) != p || ncol(X.ind) != num.X) {
      stop("X.ind must be a matrix with the number of rows equal to the # of columns in y and the number of columns in X")
    }
    if (!all(X.ind %in% c(0, 1))) {
      stop("All elements of X.ind must either equal to 1 or 0, corresponding to whether a particular covariate is included or excluded for a particular column response, respectively")
    }
  }
  if (any(prior.control$ssvs.index != -1)) {
    message("X.ind is ignored for any columns on which SSVS is used")
  }
  
  return(X.ind)
}


check_which_traits <- function(num.traits, which.traits, traits = NULL, y = NULL, num.X, makejagsboralfile_messages = FALSE) {
  if (num.traits > 0 & makejagsboralfile_messages == FALSE) {
    if (num.X == 0 & num.traits > 0) {
      stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please supply X.")
    }
    if (is.null(which.traits)) {
      stop("If traits are supplied, then please also supply which.traits to inform what traits are regressed against which covariates.")
    }
    if (nrow(traits) != ncol(y)) {
      stop("If traits are supplied, then please ensure the number of rows in traits i.e., number of species, is equal to the number of columns in y.")
    }
    if ((num.X + 1) != length(which.traits)) {
      stop("which.traits should have equal to 1+ncol(X).")
    }
    if (any(sapply(which.traits, length) > num.traits)) {
      stop("Each element in the list which.traits should have at most ncol(traits) elements.")
    }
    if (any(sapply(which.traits, function(x) any(x > ncol(traits))))) {
      stop("The values contained in the list which.traits can be takes from 1 to ncol(traits).")
    }
  }
  
  
  if (num.traits > 0 & makejagsboralfile_messages == TRUE) {
    if (num.X == 0) {
      stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please set num.X > 0.")
    }
    if (is.null(which.traits)) {
      stop("If num.traits > 0, then please supply which.traits to inform what traits are regressed against which covariates.")
    }
    if (!is.null(which.traits) & ((num.X + 1) != length(which.traits))) {
      stop("which.traits should be a list with length 1+num.X.")
    }
    if (!is.null(which.traits) & any(sapply(which.traits, length) > num.traits)) {
      stop("Each element in the list which.traits should have at most num.traits elements.")
    }
  }
}

###############################
## Unseen functions
###############################

## Construct strings for the prior distributions used in makejagsboralmodel and makejagsboralnullmodel
construct.prior.strings <- function(x) {
  x$type[1] <- match.arg(x$type[1], choices = c("normal", "cauchy", "uniform"))
  x$type[2] <- match.arg(x$type[2], choices = c("normal", "cauchy", "uniform"))
  x$type[3] <- match.arg(x$type[3], choices = c("normal", "cauchy", "uniform"))
  x$type[4] <- match.arg(x$type[4], choices = c("halfnormal", "halfcauchy", "uniform"))
  
  
  if (x$type[1] == "normal") {
    prior.string1 <- paste0("dnorm(0,", 1 / x$hypparams[1], ")")
  }
  if (x$type[1] == "cauchy") {
    prior.string1 <- paste0("dt(0,", 1 / x$hypparams[1], ",1)")
  }
  if (x$type[1] == "uniform") {
    prior.string1 <- paste0("dunif(-", x$hypparams[1], ",", x$hypparams[1], ")")
  }
  
  if (x$type[2] == "normal") {
    prior.string2 <- paste0("dnorm(0,", 1 / x$hypparams[2], ")")
    prior.string22 <- paste0("dnorm(0,", 1 / x$hypparams[2], ")I(0,)")
  }
  if (x$type[2] == "cauchy") {
    prior.string2 <- paste0("dt(0,", 1 / x$hypparams[2], ",1)")
    prior.string22 <- paste0("dt(0,", 1 / x$hypparams[4], ",1)I(0,)")
  }
  if (x$type[2] == "uniform") {
    prior.string2 <- paste0("dunif(-", x$hypparams[2], ",", x$hypparams[2], ")")
    prior.string22 <- paste0("dunif(0,", x$hypparams[4], ")")
  }
  
  if (x$type[3] == "normal") {
    prior.string3 <- paste0("dnorm(0,", 1 / x$hypparams[3], ")")
  }
  if (x$type[3] == "cauchy") {
    prior.string3 <- paste0("dt(0,", 1 / x$hypparams[3], ",1)")
  }
  if (x$type[3] == "uniform") {
    prior.string3 <- paste0("dunif(-", x$hypparams[3], ",", x$hypparams[3], ")")
  }
  
  if (x$type[4] == "uniform") {
    prior.string4 <- paste0("dunif(0,", x$hypparams[4], ")")
  }
  if (x$type[4] == "halfcauchy") {
    prior.string4 <- paste0("dt(0,", 1 / x$hypparams[4], ",1)I(0,)")
  }
  if (x$type[4] == "halfnormal") {
    prior.string4 <- paste0("dnorm(0,", 1 / x$hypparams[4], ")I(0,)")
  }
  # if(x$type[4] == "gamma") prior.string4 <- paste0("dgamma(",1/x$hypparams[4],",",1/x$hypparams[4],")")
  
  return(list(p1 = prior.string1, p2 = prior.string2, p22 = prior.string22, p3 = prior.string3, p4 = prior.string4))
}


## Fill in empty components of prior.control and mcmc.control
fillin_prior_control <- function(x) {
  if (!("type" %in% names(x))) {
    x$type <- c("normal", "normal", "normal", "uniform")
  }
  if (!("hypparams" %in% names(x))) {
    x$hypparams <- c(10, 10, 10, 30)
  }
  if (!("ssvs.index" %in% names(x))) {
    x$ssvs.index <- -1
  }
  if (!("ssvs.g" %in% names(x))) {
    x$ssvs.g <- 1e-6
  }
  if (!("ssvs.traitsindex" %in% names(x))) {
    x$ssvs.traitsindex <- -1
  }
  
  return(x)
}

fillin.mcmc.control <- function(x) {
  if (!("n.burnin" %in% names(x))) {
    x$n.burnin <- 10000
  }
  if (!("n.iteration" %in% names(x))) {
    x$n.iteration <- 40000
  }
  if (!("n.thin" %in% names(x))) {
    x$n.thin <- 30
  }
  if (!("seed" %in% names(x))) {
    x$seed <- 123
  }
  
  return(x)
}


## Given the lv, coefficients and cutoffs, return the multinomial probabilities for proportional odds regression for specific column of y
ordinal.conversion.spp <- function(n, lv = NULL, lv.coefs.j = NULL, num.lv = NULL,
                                   row.coefs = NULL, row.ids = NULL, X = NULL, X.coefs.j = NULL, offset.j = NULL, cutoffs, est = "median") {
  if (is.null(num.lv)) {
    num.lv <- 0
  }
  
  etas <- matrix(0, nrow = n, ncol = length(cutoffs)) ## num.ord.levels - 1
  for (k in 1:length(cutoffs)) {
    etas[, k] <- rep(cutoffs[k], n) - lv.coefs.j[1]
    if (!is.null(lv)) {
      etas[, k] <- etas[, k] - as.matrix(lv) %*% lv.coefs.j[2:(num.lv + 1)]
    }
    if (!is.null(row.coefs)) {
      if (est == "median") {
        for (k2 in 1:ncol(row.ids)) {
          etas[, k] <- etas[, k] - row.coefs[[k2]]$median[row.ids[, k2]]
        }
      }
      if (est == "mean") {
        for (k2 in 1:ncol(row.ids)) {
          etas[, k] <- etas[, k] - row.coefs[[k2]]$mean[row.ids[, k2]]
        }
      }
      if (est == "ignore") {
        for (k2 in 1:ncol(row.ids)) {
          etas[, k] <- etas[, k] - row.coefs[[k2]][row.ids[, k2]]
        }
      }
    }
    if (!is.null(offset.j)) {
      etas[, k] <- etas[, k] - matrix(offset.j, ncol = 1)
    }
    if (!is.null(X)) {
      etas[, k] <- etas[, k] - as.matrix(X) %*% X.coefs.j
    } ## Don't forget the negative sign!
  }
  probs <- matrix(0, n, length(cutoffs) + 1) ## num.ord.levels
  probs[, 1] <- pnorm(etas[, 1])
  for (k in 2:ncol(etas)) {
    probs[, k] <- pnorm(etas[, k]) - pnorm(etas[, k - 1])
  }
  probs[, length(cutoffs) + 1] <- 1 - pnorm(etas[, length(cutoffs)])
  
  rm(etas)
  return(probs)
}


## Process Geweke's convergence diagnotics from a single chain MCMC fit
process.geweke <- function(fit.mcmc, y, X = NULL, traits = NULL, family, num.lv,
                           row.eff, row.ids, num.ord.levels = NULL, prior.control) { # type = "independent"
  
  p <- ncol(y)
  num.X <- 0
  if (!is.null(X)) {
    num.X <- ncol(X)
  }
  num.traits <- 0
  if (!is.null(traits)) {
    num.traits <- ncol(traits)
  }
  
  fit_geweke <- geweke.diag(fit.mcmc)[[1]] ## Takes first chain only
  out_gewekelist <- list(lv.coefs = matrix(fit_geweke[grep("lv.coefs", names(fit_geweke))], nrow = p))
  if (num.lv > 0) {
    out_gewekelist$lv.coefs <- matrix(out_gewekelist$lv.coefs[, -c(2:(num.lv + 1))], nrow = p) ## Drop check on LV coefs
    fit_geweke <- fit_geweke[-grep("lvs", names(fit_geweke))] ## Drop check on lv
    # 			if(type != "independent")
    #                     out_gewekelist$lv.covparams <- fit_geweke[grep("lv.covparams",names(fit_geweke))]
  }
  rownames(out_gewekelist$lv.coefs) <- colnames(y)
  if (ncol(out_gewekelist$lv.coefs) > 1) {
    colnames(out_gewekelist$lv.coefs) <- c("beta0", "Disperson")
  }
  else {
    colnames(out_gewekelist$lv.coefs) <- c("beta0")
  }
  
  if (row.eff != "none") {
    out_gewekelist$row.coefs <- vector("list", ncol(row.ids))
    names(out_gewekelist$row.coefs) <- colnames(row.ids)
    for (k in 1:ncol(row.ids)) {
      out_gewekelist$row.coefs[[k]] <- fit_geweke[grep(paste0("row.coefs.ID", k, "\\["), names(fit_geweke))]
    }
  }
  if (row.eff == "random") {
    out_gewekelist$row.sigma <- vector("list", ncol(row.ids))
    names(out_gewekelist$row.sigma) <- colnames(row.ids)
    for (k in 1:ncol(row.ids)) {
      out_gewekelist$row.sigma[[k]] <- fit_geweke[grep(paste0("row.sigma.ID", k, "$"), names(fit_geweke))]
    }
  }
  
  if (num.X > 0) {
    out_gewekelist$X.coefs <- matrix(fit_geweke[grep("X.coefs", names(fit_geweke))], nrow = p)
    rownames(out_gewekelist$X.coefs) <- colnames(y)
    colnames(out_gewekelist$X.coefs) <- colnames(X)
    if (any(prior.control$ssvs.index > -1)) {
      out_gewekelist$X.coefs <- NULL
    } ## Drop check on X coefs if SSVS is used
  }
  
  if (num.traits > 0) {
    out_gewekelist$traits.coefs <- cbind(
      fit_geweke[grep("traits.int", names(fit_geweke))],
      matrix(fit_geweke[grep("traits.coefs", names(fit_geweke))], nrow = ncol(X) + 1, ncol = ncol(traits)),
      fit_geweke[grep("trait.sigma", names(fit_geweke))]
    )
    rownames(out_gewekelist$traits.coefs) <- c("beta0", colnames(X))
    colnames(out_gewekelist$traits.coefs) <- c("kappa0", colnames(traits), "sigma")
    if (any(unlist(prior.control$ssvs.index) > -1)) {
      out_gewekelist$traits.coefs <- NULL
    } ## Drop check on trait parameters if SSVS is used
  }
  
  if (any(family == "ordinal")) {
    out_gewekelist$cutoffs <- fit_geweke[grep("cutoffs", names(fit_geweke))]
    names(out_gewekelist$cutoffs) <- paste(1:(num.ord.levels - 1), "|", 2:num.ord.levels, sep = "")
    
    if (sum(family == "ordinal") > 2) {
      out_gewekelist$ordinal.sigma <- fit_geweke[grep("ordinal.sigma", names(fit_geweke))]
    }
  }
  
  if (any(family == "tweedie")) {
    out_gewekelist$powerparam <- fit_geweke[grep("powerparam", names(fit_geweke))]
  }
  
  
  prop.outside <- table(2 * pnorm(abs(unlist(out_gewekelist)), lower.tail = FALSE) < 0.05) / length(unlist(out_gewekelist))
  
  return(list(geweke.diag = out_gewekelist, prop.exceed = prop.outside))
}


## Sets up part of the JAGS script corresponding to family for responses; used in make.jagsboralmodel and make.jagsboralnullmodel.
setup.resp.families.lv <- function(p, complete.family, num.lv, row.eff, row.ids,
                                   offset, num.X, complete.trial.size, index.tweed.cols, index.ord.cols) {
  respfamily_script <- NULL
  
  for (j in 1:p) {
    if (complete.family[j] != "multinom") {
      if (length(unique(complete.family)) == 1) {
        if (j == 1) {
          if (num.lv == 0) {
            linpred_string <- paste("eta[i,j] <- 0", sep = "")
          }
          if (num.lv > 0) {
            linpred_string <- paste("eta[i,j] <- inprod(lv.coefs[j,2:(num.lv+1)],lvs[i,])", sep = "")
          }
          if (row.eff != "none") {
            for (k in 1:ncol(row.ids)) {
              linpred_string <- paste(linpred_string, " + row.coefs.ID", k, "[row.ids[i,", k, "]]", sep = "")
            }
          }
          if (num.X > 0) {
            linpred_string <- paste(linpred_string, " + inprod(X.coefs[j,],X[i,])", sep = "")
          }
          if (!is.null(offset)) {
            linpred_string <- paste(linpred_string, " + offset[i,j]", sep = "")
          }
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { ", linpred_string, " }", sep = ""))
        }
        if (j > 1) { }
      }
      if (length(unique(complete.family)) > 1) {
        if (num.lv == 0) {
          linpred_string <- paste("eta[i,", j, "] <- 0", sep = "")
        }
        if (num.lv > 0) {
          linpred_string <- paste("eta[i,", j, "] <- inprod(lv.coefs[", j, ",2:(num.lv+1)],lvs[i,])", sep = "")
        }
        if (row.eff != "none") {
          for (k in 1:ncol(row.ids)) {
            linpred_string <- paste(linpred_string, " + row.coefs.ID", k, "[row.ids[i,", k, "]]", sep = "")
          }
        }
        if (num.X > 0) {
          linpred_string <- paste(linpred_string, " + inprod(X.coefs[", j, ",],X[i,])", sep = "")
        }
        if (!is.null(offset)) {
          linpred_string <- paste(linpred_string, " + offset[i,", j, "]", sep = "")
        }
        respfamily_script <- c(respfamily_script, paste("\t\t ", linpred_string, sep = ""))
      }
    }
    
    if (complete.family[j] == "negative.binomial") {
      if (length(unique(complete.family)) == 1) {
        if (j == 1) {
          # respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { u[i,j] ~ dnorm(0, 1/lv.coefs[",j, ",2]) }", sep = ""))
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { u[i,j] ~ dgamma(1/lv.coefs[j,num.lv+2], 1/lv.coefs[j,num.lv+2]) }", sep = ""))
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { y[i,j] ~ dpois(exp(lv.coefs[j,1] + eta[i,j ])*(u[i,j])) } ## Parameterizing the NB as a multiplicative random effect models\n", sep = ""))
        }
        if (j > 1) { }
      }
      if (length(unique(complete.family)) > 1) {
        # respfamily_script <- c(respfamily_script, paste("\t\t u[i,",j, "] ~ dnorm(0, 1/lv.coefs[",j, ",2])", sep = ""))
        respfamily_script <- c(respfamily_script, paste("\t\t u[i,", j, "] ~ dgamma(1/lv.coefs[", j, ",num.lv+2], 1/lv.coefs[", j, ",num.lv+2])", sep = ""))
        respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dpois(exp(lv.coefs[", j, ",1] + eta[i,", j, "])*(u[i,", j, "])) ## Parameterizing the NB as a multiplicative random effect models, with size\n", sep = ""))
      }
    }
    
    # 		if(complete.family[j] == "negative.binomial2") {
    # 			if(length(unique(complete.family)) == 1) {
    # 				if(j == 1) {
    # 					respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { u[i,j] <- 1/(1 + lv.coefs[j,num.lv+2]*exp(lv.coefs[j,1] + eta[i,j])) }", sep = ""))
    # 					respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { y[i,j] ~ dnegbin(u[i,j],1/lv.coefs[j,num.lv+2]) } ## Parameterizing the NB as a function of prob and size \n", sep = ""))
    # 					}
    # 				if(j > 1) { }
    # 				}
    # 			if(length(unique(complete.family)) > 1) {
    # 				respfamily_script <- c(respfamily_script, paste("\t\t u[i,",j, "] <- 1/(1 + lv.coefs[",j, ",num.lv+2]*exp(lv.coefs[", j, ",1] + eta[i,",j, "]))", sep = ""))
    # 				respfamily_script <- c(respfamily_script, paste("\t\t y[i,",j, "] ~ dnegbin(u[i,", j, "],1/lv.coefs[",j, ",num.lv+2]) ## Parameterizing the NB as a function of prob and size \n", sep = ""))
    # 				}
    # 			}
    
    if (complete.family[j] == "normal") {
      if (length(unique(complete.family)) == 1) {
        if (j == 1) {
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { y[i,j] ~ dnorm(lv.coefs[j,1] + eta[i,j],pow(lv.coefs[j,num.lv+2],-2)) } \n", sep = ""))
        }
        if (j > 1) { }
      }
      if (length(unique(complete.family)) > 1) {
        respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dnorm(lv.coefs[", j, ",1] + eta[i,", j, "],pow(lv.coefs[", j, ",num.lv+2],-2)) \n", sep = ""))
      }
    }
    
    #         if(all(complete.family == "bernoulli")) { ## If all data are Bernoulli, then use step parameterization to speed up sampling!
    #             if(length(unique(complete.family)) == 1) {
    #                 if(j == 1) {
    #                     respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { Z[i,j] ~ dnorm(lv.coefs[j,1] + eta[i,j], 1) }", sep = ""))
    #                     respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { y[i,j] ~ dbern(step(Z[i,j])) }\n", sep = ""))
    #                     }
    #                 if(j > 1) { }
    #                 }
    #             if(length(unique(complete.family)) > 1) {
    #                 respfamily_script <- c(respfamily_script, paste("\t\t Z[i,",j, "] ~ dnorm(lv.coefs[", j, ",1] + eta[i,",j, "], 1)", sep = ""))
    #                 respfamily_script <- c(respfamily_script, paste("\t\t y[i,",j, "] ~ dbern(step(Z[i,",j, "]))\n", sep = ""))
    #                 }
    #             }
    
    if (complete.family[j] == "binomial") {
      respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dbin(phi(lv.coefs[", j, ",1] + eta[i,", j, "]),", complete.trial.size[j], ")\n", sep = ""))
    }
    
    if (complete.family[j] == "exponential") {
      if (length(unique(complete.family)) == 1) {
        if (j == 1) {
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { y[i,j] ~ dexp(pow(exp(lv.coefs[j,1] + eta[i,j]),-1)) }\n", sep = ""))
        }
        if (j > 1) { }
      }
      if (length(unique(complete.family)) > 1) {
        respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dexp(pow(exp(lv.coefs[", j, ",1] + eta[i,", j, "]),-1))\n", sep = ""))
      }
    }
    
    if (complete.family[j] == "gamma") {
      if (length(unique(complete.family)) == 1) {
        if (j == 1) {
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { y[i,j] ~ dgamma(exp(lv.coefs[j,1] + eta[i,j])*lv.coefs[j,num.lv+2], lv.coefs[j,num.lv+2]) } \n", sep = ""))
        }
        if (j > 1) { }
      }
      if (length(unique(complete.family)) > 1) {
        respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dgamma(exp(lv.coefs[", j, ",1] + eta[i,", j, "])*lv.coefs[", j, ",num.lv+2], lv.coefs[", j, ",num.lv+2])\n", sep = ""))
      }
    }
    
    if (complete.family[j] == "beta") {
      if (length(unique(complete.family)) == 1) {
        if (j == 1) {
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { y[i,j] ~ dbeta(ilogit(lv.coefs[j,1] + eta[i,j])*lv.coefs[j,num.lv+2],(1-ilogit(lv.coefs[j,1] + eta[i,j]))*lv.coefs[j,num.lv+2]) }\n", sep = ""))
        }
        if (j > 1) { }
      }
      if (length(unique(complete.family)) > 1) {
        respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dbeta(ilogit(lv.coefs[", j, ",1] + eta[i,", j, "])*lv.coefs[", j, ",num.lv+2],(1-ilogit(lv.coefs[", j, ",1] + eta[i,", j, "]))*lv.coefs[", j, ",num.lv+2])\n", sep = ""))
      }
    }
    
    if (complete.family[j] == "poisson") {
      if (length(unique(complete.family)) == 1) {
        if (j == 1) {
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { y[i,j] ~ dpois(exp(lv.coefs[j,1] + eta[i,j])) }\n", sep = ""))
        }
        if (j > 1) { }
      }
      if (length(unique(complete.family)) > 1) {
        respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dpois(exp(lv.coefs[", j, ",1] + eta[i,", j, "]))\n", sep = ""))
      }
    }
    
    if (complete.family[j] == "lnormal") {
      if (length(unique(complete.family)) == 1) {
        if (j == 1) {
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { y[i,j] ~ dlnorm(lv.coefs[j,1] + eta[i,j],pow(lv.coefs[j,num.lv+2],-2)) } \n", sep = ""))
        }
        if (j > 1) { }
      }
      if (length(unique(complete.family)) > 1) {
        respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dlnorm(lv.coefs[", j, ",1] + eta[i,", j, "],pow(lv.coefs[", j, ",num.lv+2],-2)) \n", sep = ""))
      }
    }
    
    if (complete.family[j] == "tweedie") {
      respfamily_script <- c(respfamily_script, paste("\t\t lambdanum[i,", which(index.tweed.cols == j), "] <- pow(exp(lv.coefs[", j, ",1] + eta[i,", j, "]),2-powerparam)/(lv.coefs[", j, ",num.lv+2]*(2-powerparam))", sep = ""))
      respfamily_script <- c(respfamily_script, paste("\t\t numfish[i,", which(index.tweed.cols == j), "] ~ dpois(lambdanum[i,", which(index.tweed.cols == j), "])", sep = ""))
      respfamily_script <- c(respfamily_script, paste("\t\t choose.shape[i,", which(index.tweed.cols == j), ",1] <- numfish[i,", which(index.tweed.cols == j), "]*(2-powerparam)/(powerparam-1)", sep = "")) ## If y > 0, then conditional on numfish, y is sum of independent gammas
      respfamily_script <- c(respfamily_script, paste("\t\t choose.rate[i,", which(index.tweed.cols == j), ",1] <- 1/(lv.coefs[", j, ",num.lv+2]*(powerparam-1)*pow(exp(lv.coefs[", j, ",1] + eta[i,", j, "]),powerparam-1))", sep = ""))
      respfamily_script <- c(respfamily_script, paste("\t\t choose.shape[i,", which(index.tweed.cols == j), ",2] <- 1", sep = "")) ## If y = 0, then Tweedie dist equals probability of Poisson equal 0
      respfamily_script <- c(respfamily_script, paste("\t\t choose.rate[i,", which(index.tweed.cols == j), ",2] <- exp(-lambdanum[i,", which(index.tweed.cols == j), "])", sep = ""))
      respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dgamma(choose.shape[i,", which(index.tweed.cols == j), ",1+equals(y[i,", which(index.tweed.cols == j), "],0)],choose.rate[i,", j, ",1+equals(y[i,", j, "],0)]) \n", sep = ""))
    }
    
    if (complete.family[j] == "ordinal") {
      if (length(index.ord.cols) == p) {
        if (j == 1) {
          respfamily_script <- c(respfamily_script, paste("\t\t for(j in 1:p) { \n\t\t\t prob[i,j,1] <- phi(cutoffs[1]-eta[i,j]-lv.coefs[j,1])", sep = ""))
          respfamily_script <- c(respfamily_script, paste("\t\t\t for(k in 2:(num.ord.levels-1)) { \n\t\t\t\t prob[i,j,k] <- phi(cutoffs[k]-eta[i,j]-lv.coefs[j,1]) - phi(cutoffs[k-1]-eta[i,j]-lv.coefs[j,1]) \n\t\t\t\t }", sep = ""))
          respfamily_script <- c(respfamily_script, paste("\t\t\t prob[i,j,num.ord.levels] <- 1-phi(cutoffs[num.ord.levels-1]-eta[i,j]-lv.coefs[j,1])", sep = ""))
          respfamily_script <- c(respfamily_script, paste("\t\t\t y[i,j] ~ dcat(prob[i,j,]) \n\t\t\t } \n", sep = ""))
        }
        if (j > 1) { }
      }
      if (length(index.ord.cols) < p) {
        respfamily_script <- c(respfamily_script, paste("\t\t prob[i,", which(index.ord.cols == j), ",1] <- phi(cutoffs[1]-eta[i,", j, "]-lv.coefs[", j, ",1])", sep = ""))
        respfamily_script <- c(respfamily_script, paste("\t\t for(k in 2:(num.ord.levels-1)) { prob[i,", which(index.ord.cols == j), ",k] <- phi(cutoffs[k]-eta[i,", j, "]-lv.coefs[", j, ",1]) - phi(cutoffs[k-1]-eta[i,", j, "]-lv.coefs[", j, ",1]) }", sep = ""))
        respfamily_script <- c(respfamily_script, paste("\t\t prob[i,", which(index.ord.cols == j), ",num.ord.levels] <- 1-phi(cutoffs[num.ord.levels-1]-eta[i,", j, "]-lv.coefs[", j, ",1])", sep = ""))
        respfamily_script <- c(respfamily_script, paste("\t\t y[i,", j, "] ~ dcat(prob[i,", which(index.ord.cols == j), ",])\n", sep = ""))
      }
    }
    
    if (complete.family[j] == "multinom") {
      stop("You shouldn't have gotten here!") ## Coefficients for lv are constrained to be same for all levels! Otherwise identifiability constraints are hard!
      #    		model_script <- c(model_script, paste("\t\t for(k in 1:num.multinom.levels[",j,"]) {",sep=""))
      # 			if(num.X == 0 & row.eff) model_script <- c(model_script, paste("\t\t\t mu[i,",which(index_multinom_cols == j),",k] <- exp(row.coefs[i] + inprod(lv.coefs[",j,",2:(num.lv+1)],lvs[i,]))",sep=""))
      # 			if(num.X > 0 & row.eff) model_script <- c(model_script, paste("\t\t\t mu[i,",which(index_multinom_cols == j),",k] <- exp(row.coefs[i] + inprod(lv.coefs[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.multinom.params[",which(index_multinom_cols == j),",,k],X[i,]))",sep=""))
      # 			if(num.X == 0 & !row.eff) model_script <- c(model_script, paste("\t\t\t mu[i,",which(index_multinom_cols == j),",k] <- exp(inprod(lv.coefs[",j,",2:(num.lv+1)],lvs[i,]))",sep=""))
      # 			if(num.X > 0 & !row.eff) model_script <- c(model_script, paste("\t\t\t mu[i,",which(index_multinom_cols == j),",k] <- exp(inprod(lv.coefs[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.multinom.params[",which(index_multinom_cols == j),",,k],X[i,]))",sep=""))
      # 			model_script <- c(model_script, paste("\t\t\t prob[i,",which(index_multinom_cols == j),",k] <- mu[i,",which(index_multinom_cols == j),",k]/sum(mu[i,",which(index_multinom_cols == j),",]) }",sep=""))
      # 			model_script <- c(model_script, paste("\t\t y[i,",j,"] ~ dcat(prob[i,",which(index_multinom_cols == j),",]+0.001)\n",sep=""))
    }
  }
  
  return(respfamily_script)
}


function() {
  ## Extract rhats from multiple chained MCMC fit
  rhats <- function(x, asc = FALSE) {
    if (asc) {
      x$BUGSoutput$summary[order(x$BUGSoutput$summary[, "Rhat"]), "Rhat", drop = FALSE]
    }
    else {
      x$BUGSoutput$summary[, "Rhat", drop = FALSE]
    }
  }
  
  
  # ## Process the rhats from multiple chained MCMC fit
  # process.rhats <- function(sims.matrix) {
  #  		combined_fit_mcmc <- as.mcmc(sims.matrix)
  # 		fit.rhats <- rhats(jagsfit, asc = FALSE)
  #     		make.rhatslist <- list(lv.coefs = matrix(fit.rhats[grep("lv.coefs", rownames(fit.rhats))], nrow = p))
  #     		#if(num.lv > 0) { fit.rhats <- fit.rhats[-grep("lvs",rownames(fit.rhats)),] } ## Drop check on lv
  #     		if(num.lv > 0) { make.rhatslist$lv.coefs <- as.matrix(make.rhatslist$lv.coefs[,-c(2:(num.lv+1))]) } ## Drop check on LV coefs
  #     		rownames(make.rhatslist$lv.coefs) <- colnames(y);
  #     		if(ncol(make.rhatslist$lv.coefs) > 1) { colnames(make.rhatslist$lv.coefs) <- c("beta0","Disperson") }
  # 		else { colnames(make.rhatslist$lv.coefs) <- c("beta0") }
  #
  #     		if(row.eff != "none") {
  #     			make.rhatslist$row.coefs <- fit.rhats[grep("row.coefs", rownames(fit.rhats))]
  #     			names(make.rhatslist$row.coefs) <- rownames(y)
  #     			}
  #     		if(row.eff == "random") {
  #     			make.rhatslist$row.sigma <- fit.rhats[grep("row.sigma.ID", rownames(fit.rhats))]
  #    				names(make.rhatslist$row.sigma) <- c("Row random effects sigma")
  #     			}
  #
  #    		if(num.X > 0) {
  #     			make.rhatslist$X.coefs <- matrix(fit.rhats[grep("X.coefs", rownames(fit.rhats))], nrow = p)
  #     			rownames(make.rhatslist$X.coefs) <- colnames(y); colnames(make.rhatslist$X.coefs) <- colnames(X)
  #     			}
  #
  # 		if(num.traits > 0) {
  #     			make.rhatslist$traits.coefs <- cbind(fit.rhats[grep("trait.int", rownames(fit.rhats))], matrix(fit.rhats[grep("traits.coefs", rownames(fit.rhats))], nrow = ncol(X)+1), fit.rhats[grep("trait.sigma", rownames(fit.rhats))])
  #     			rownames(make.rhatslist$traits.coefs) <- c("beta0",colnames(X)); colnames(make.rhatslist$traits.coefs) <- c("kappa0",colnames(traits),"sigma")
  #     			}
  #
  #    		if(any(family == "ordinal")) {
  #    			make.rhatslist$cutoffs <- fit.rhats[grep("cutoffs", rownames(fit.rhats))]
  #    			names(make.rhatslist$cutoffs) <- paste(1:(num.ord.levels - 1), "|", 2:num.ord.levels, sep = "")
  #  			if(sum(family == "ordinal") > 2) {
  #  				make.rhatslist$ordinal.sigma <- fit.rhats[grep("ordinal.sigma", rownames(fit.rhats))]
  #  				names(make.rhatslist$ordinal.sigma) <- "Species-specific random intercept sigma for ordinal responses"
  #  				}
  #    			}
  #
  #    		if(any(family == "tweedie")) {
  #    			make.rhatslist$powerparam <- fit.rhats[grep("powerparam", rownames(fit.rhats))]
  #    			names(make.rhatslist$powerparam) <- "Common power parameter"
  #    			}
  #
  # 	return(make.rhatslist)
  # 	}
}