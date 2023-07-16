require(tictoc);


## PATS function
# Arguments:
# ----------
# outcome = outcome variable
# partial.blip.mod = a list where each element represents the right hand side of the
#                    partial blip function, 
#                    for example list(~X11, ~X21); if treatment is to tailored
#                    according to X11 at first time point and X21 at second (H*).
# full.blip.mod = a list where each element represents the right hand size of the
#                 complete blip function (H), for example list(~X11+X12,~X21+X22+A1+X11+X12);
# partial.treat.mod = a list where each element represents the regression of the
#                     treatment according to H*. For example list(A1~X11, A2~X21);
# full.treat.mod = a list where each element represent the regression of the treatment
#                  according to H. For example list(A1~X11+X12, A2~X21+X22+A1+X12+X12);
# tf.mod = a list where each element is the right hand side of the treatment free model.
#          list(~X11+X12, ~X21+X22+A1+X11+X12);
# data = either NULL if all variables can be found in the global environment,
#        or the name of a data frame where variables can be found.
# method = "iptw dwols", "iptw gest", "ce dwols", "ce gest", "all"
# weight = "default" (overlap) or "iptw"
# var.estim = "bootstrap" (fixed m-out-of-n bootstrap or n-out-of-n bootstrap),
#             "none" or "adapt" (adaptive m-out-of-n bootstrap)
# B = The number of first stage bootstrap samples
# B2 = The number of second stage bootstrap samples
# M = numeric value representing a fixed value for the m-out-of-n bootstrap.
#     M = 0 is used to request the usual nonparametric bootstrap (n-out-of-n).
# verbose = TRUE or FALSE. If TRUE, the function informs the user of the current
#           running time and expected total running time.
# interrupt = TRUE or FALSE. If TRUE, the function offers the user to interrupt
#             the function if the expected running time is ­> 10 minutes.
#
# Value:
# ------
# opt.Y = The estimated expected outcome value under the estimated optimal PATS
# regret = The difference between the observed outcome and opt.Y
# covmat = The estimated variance-covariance matrix of the psi* parameters
# psi = The estimated psi* parameters
# nonreg = The estimated non-regularity, the p parameter in the adaptive m-out-of-n bootstrap 
# M = The m parameter in the adaptive m-out-of-n bootstrap
# alpha = The alpha parameter in the adaptive m-out-of-n bootstrap
# psi.boot = The bootstrap replicates of psi
# 
#
# Details:
# --------
# A must be binary 0/1
# An example is available at the end of the file
# 




PATS = function(outcome, partial.blip.mod, full.blip.mod, partial.treat.mod = NULL, full.treat.mod,
                tf.mod, data = NULL, method = "ce dwols", weight = "default",
                var.estim = "none", B = 200, B2 = 200, M = 0, verbose = FALSE, interrupt = FALSE){

  t0 = proc.time();
  expit = plogis;

  ## Error checks
  if(!is.character(method)) stop("method must be either 'iptw dwols', 'iptw gest', 'ce dwols' or 'ce gest'");
  if(length(method) != 1) stop("Only one method can be supplied");
  if(!(method %in% c("iptw dwols", "iptw gest", "ce dwols", "ce gest"))) stop("method must be either 'iptw dwols', 'iptw gest', 'ce dwols' or 'ce gest'");
  if(method %in% c("ce dwols", "ce gest")) partial.treat.mod = full.treat.mod;
  if(!is.list(partial.blip.mod)) stop("partial.blip.mod must be a list");
  if(!is.list(full.blip.mod)) stop("full.blip.mod must be a list");
  if(!is.list(partial.treat.mod)) stop("partial.treat.mod must be a list");
  if(!is.list(full.treat.mod)) stop("full.treat.mod must be a list");
  if(!is.list(tf.mod)) stop("tf.mod must be a list");
  if(!is.null(data) & !is.data.frame(data)) stop("data must either be NULL or a data frame");
  if (is.null(data)) {
        data <- cbind(outcome, do.call("cbind", lapply(partial.blip.mod, 
            get_all_vars)), do.call("cbind", lapply(full.blip.mod, 
            get_all_vars)), do.call("cbind", lapply(partial.treat.mod, 
            get_all_vars)), do.call("cbind", lapply(full.treat.mod, 
            get_all_vars)), do.call("cbind", lapply(tf.mod, get_all_vars)))
        data <- data[!duplicated(names(data))]
  } 
  if(!is.numeric(data$outcome)) stop("outcome must be a numeric variable");
  if(any(!sapply(partial.blip.mod, inherits, "formula"))) stop("partial.blip must be a list of formula objects");
  if(any(!sapply(full.blip.mod, inherits, "formula"))) stop("full.blip.mod must be a list of formula objects");
  if(any(!sapply(partial.treat.mod, inherits, "formula"))) stop("partial.treat.mod must be a list of formula objects");
  if(any(!sapply(full.treat.mod, inherits, "formula"))) stop("full.treat.mod must be a list of formula objects");
  if(any(!sapply(tf.mod, inherits, "formula"))) stop("tf.mod must be a list of formula objects");
  if(anyNA(data)) stop("Missing values are not currently allowed in this function");
  if(weight != "default" & weight != "iptw") stop("weight must either be default or iptw");
  if(var.estim != "bootstrap" & var.estim != "none" & var.estim != "adapt") stop("var.estim must either be 'bootstrap', 'none' or 'adapt'");
  if(var.estim == "boostrap" | var.estim == "adapt"){
    if(!is.numeric(B)) stop("B must be numeric");
    if(!is.numeric(M) & var.estim == "bootstrap") stop("M must be numeric");
    if(B < 2) stop ("B must be greater than 2");
    if(B%%1 != 0) {B = B - B%%1;
                      warning("B was not an integer; it has been rounded down.");};
    if(B2 < 2 & var.estim == "adapt") stop ("B2 must be greater than 2");
    if(B2%%1 != 0 & var.estim == "adapt") {B2 = B2 - B2%%1;
                      warning("B2 was not an integer; it has been rounded down.");};
    if(is.numeric(M)){
      if(M%%1 != 0){M = M - M%%1;
                      warning("M was not an integer; it has been rounded down.");};
    }
  }
  if(!(verbose == TRUE | verbose == FALSE)) stop("verbose must either be TRUE or FALSE");
  if(!(interrupt == TRUE | interrupt == FALSE)) stop("interrupt must either be TRUE or FALSE");
  if(identical(partial.blip.mod, full.blip.mod)) stop("The full blip model is the same as the partial blip; there is no point in using the PATS function");
  if(length(partial.blip.mod) != length(full.blip.mod) | 
     length(partial.blip.mod) != length(full.treat.mod) |
     length(partial.blip.mod) != length(partial.treat.mod) |
     length(partial.blip.mod) != length(tf.mod)) stop("At least two of the lists of models that were supplied are not of the same length");


  ## Number of follow-up times
  K = length(partial.blip.mod);


  if(K == 1 & var.estim == "adapt") M = n;


  if(K > 2 & var.estim == "adapt") stop("Adaptive m-out-of-n bootstrap has not been developped for more than two time-points.");

  analyze.f = function(data){
    n = nrow(data);

    ## Find the exposure variable
    A = lapply(full.treat.mod, function(x){model.response(model.frame(x, data = data))});


    ## Verify if A is coded 0/1
    for(j in 1:K){
      if(length(unique(A[[j]])) != 2) stop("The treatment is not binary");
      if(max(A[[j]]) != 1 | min(A[[j]]) != 0) stop("The treatment is not coded 0/1");
    }


    obj = list();
    obj$psi = obj$opt.treat = obj$covmat = obj$regret = obj$opt.Y = list();


    ## Propensity score models
    predA = NULL;
    predA.star = NULL;
    for(j in 1:K){
      modA = glm(full.treat.mod[[j]], family = "binomial", data = data);
      predA[[j]] = modA$fitted.values;
      if(!(method %in% c("ce dwols", "ce gest"))){
        modA.star = glm(partial.treat.mod[[j]], family = "binomial", data = data);
        predA.star[[j]] = modA.star$fitted;
      }
    }


    ## Weights
    if(any(method %in% c("iptw dwols", "ce dwols"))){
      w = NULL;
      w.dWOLS.star = NULL;
      w.dWOLS = NULL;
      for(j in 1:K){
        w[[j]] = A[[j]]*predA.star[[j]]/predA[[j]] + (1 - A[[j]])*(1 - predA.star[[j]])/(1 - predA[[j]]);
        if(weight == "default"){
          w.dWOLS.star[[j]] = abs(A[[j]] - predA.star[[j]]);
          w.dWOLS[[j]] = abs(A[[j]] - predA[[j]]);
        } else{
          w.dWOLS.star[[j]] = A[[j]]/predA.star[[j]] + (1 - A[[j]])/(1 - predA.star[[j]]);
          w.dWOLS[[j]] = A[[j]]/predA[[j]] + (1 - A[[j]])/(1 - predA[[j]]);
        }
      }
      psi.dwols = NULL;
      dopt.dwols = NULL;
      regret.dwols = matrix(, nrow = n, ncol = K);
      for(j in K:1){
        if(j == K){
          Y = data$outcome;
        } else{
          Y = data$outcome + rowSums(as.matrix(regret.dwols[,(j+1):K]));
        }
        Hpsi = model.matrix(full.blip.mod[[j]], data = data);
        X = cbind(model.matrix(tf.mod[[j]], data = data),
                  A[[j]]*Hpsi);
        W = w.dWOLS[[j]];
        est = solve(t(X)%*%(W*X))%*%t(X)%*%(W*Y);
        p2 = ncol(Hpsi);
        psi.dwols[[j]] = tail(est, p2);
        dopt.dwols[[j]] = 1*(Hpsi%*%psi.dwols[[j]] > 0);
        regret.dwols[,j] = dopt.dwols[[j]]*Hpsi%*%psi.dwols[[j]] - 
                           A[[j]]*Hpsi%*%psi.dwols[[j]];
      }
    }


    ## X(X'X)^{-1}X'
    if(any(method %in% c("iptw gest", "ce gest"))){
      XXXX = NULL;
      for(j in 1:K){
        X = model.matrix(tf.mod[[j]], data = data);
        XXXX[[j]] = X%*%solve((t(X)%*%X))%*%t(X);
      }
      psi.gest = NULL;
      dopt.gest = NULL;
      regret.gest = matrix(, nrow = n, ncol = K);
      for(j in K:1){
        if(j == K){
          Y = data$outcome;
        } else{
          Y = data$outcome + rowSums(as.matrix(regret.gest[,(j+1):K]));
        }
        Hpsi = model.matrix(full.blip.mod[[j]], data = data);
        estimating.equations = matrix(NA, nrow = ncol(Hpsi), ncol = 1 + ncol(Hpsi));
        cols = cbind(Y, -A[[j]]*Hpsi);
        rows = Hpsi;
        for(k in 1:ncol(rows)){
          for(l in 1:ncol(cols)){
            estimating.equations[k,l] = sum(((cols[,l]*(A[[j]]*rows[,k] - predA.star[[j]]*rows[,k]) - 
              (A[[j]]*rows[,k] - predA.star[[j]]*rows[,k])*XXXX[[j]]%*%cols[,l])));
          }
        }
        psi.gest[[j]] = solve(estimating.equations[,-1], -estimating.equations[,1]);
        dopt.gest[[j]] = 1*(Hpsi%*%psi.iptw.gest[[j]] > 0);
        regret.gest[,j] = dopt.gest[[j]]*Hpsi%*%psi.gest[[j]] - 
                          A[[j]]*Hpsi%*%psi.gest[[j]];
      }
    }


    ## IPTW + dWOLS
    if(any(method %in% c("iptw dwols"))){
      psi.iptw.dwols.star = NULL;
      dopt.iptw.dwols.star = NULL;
      regret.iptw.dwols.star = matrix(, nrow = n, ncol = K);
      for(j in K:1){
        if(j == K){
          Y = data$outcome;
        } else{
          Y = data$outcome + rowSums(as.matrix(regret.iptw.dwols.star[,(j+1):K]));
        }
        Hpsi = model.matrix(partial.blip.mod[[j]], data = data);
        X = cbind(model.matrix(tf.mod[[j]], data = data),
                  A[[j]]*Hpsi);
        W = diag(w.dWOLS.star[[j]]*w[[j]]);
        est = solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%Y;
        p2 = ncol(Hpsi);
        psi.iptw.dwols.star[[j]] = tail(est, p2);
        dopt.iptw.dwols.star[[j]] = 1*(Hpsi%*%psi.iptw.dwols.star[[j]] > 0);
        regret.iptw.dwols.star[,j] = dopt.iptw.dwols.star[[j]]*Hpsi%*%psi.iptw.dwols.star[[j]] - 
                                     A[[j]]*model.matrix(full.blip.mod[[j]], data = data)%*%psi.dwols[[j]];
        obj$psi[[j]] = psi.iptw.dwols.star[[j]];
        obj$opt.treat[[j]] = dopt.iptw.dwols.star[[j]];
      }
      obj$regret = rowSums(as.matrix(regret.iptw.dwols.star[,(j+1):K]));
      obj$opt.Y = data$outcome + obj$regret;
    }


    ## IPTW + G-est
    if(any(method %in% c("iptw gest"))){
      psi.iptw.gest.star = NULL;
      dopt.iptw.gest.star = NULL;
      regret.iptw.gest.star = matrix(, nrow = n, ncol = K);
      for(j in K:1){
        if(j == K){
          Y = data$outcome;
        } else{
          Y = data$outcome + rowSums(as.matrix(regret.iptw.gest.star[,(j+1):K]));
        }
        Hpsi = model.matrix(partial.blip.mod[[j]], data = data);
        estimating.equations = matrix(NA, nrow = ncol(Hpsi), ncol = 1 + ncol(Hpsi));
        cols = cbind(Y, -A[[j]]*Hpsi);
        rows = Hpsi;
        for(k in 1:ncol(rows)){
          for(l in 1:ncol(cols)){
            estimating.equations[k,l] = sum(w[[j]]*((cols[,l]*(A[[j]]*rows[,k] - predA.star[[j]]*rows[,k]) - 
              (A[[j]]*rows[,k] - predA.star[[j]]*rows[,k])*XXXX[[j]]%*%cols[,l])));
          }
        }
        psi.iptw.gest.star[[j]] = solve(estimating.equations[,-1], -estimating.equations[,1]);
        dopt.iptw.gest.star[[j]] = 1*(Hpsi%*%psi.iptw.gest.star[[j]] > 0);
        regret.iptw.gest.star[,j] = dopt.iptw.gest.star[[j]]*Hpsi%*%psi.iptw.gest.star[[j]] - 
                                    A[[j]]*model.matrix(full.blip.mod[[j]], data = data)%*%psi.gest[[j]];
        obj$psi[[j]] = psi.iptw.gest.star[[j]];
        obj$opt.treat[[j]] = dopt.iptw.gest.star[[j]];
      }
      obj$regret = rowSums(as.matrix(regret.iptw.gest.star));
      obj$opt.Y = data$outcome + obj$regret;
    }


    ## CE dWOLS
    if(any(method %in% c("ce dwols"))){
      psi.ce.dwols.star = NULL;
      dopt.ce.dwols.star = NULL;
      regret.ce.dwols.star = matrix(, nrow = n, ncol = K);
      for(j in K:1){
        if(j == K){
          Y = data$outcome;
        } else{
          Y = data$outcome + rowSums(as.matrix(regret.ce.dwols.star[,(j+1):K]));
        }
        Hpsi = model.matrix(partial.blip.mod[[j]], data = data);
        Y1_Y0 = model.matrix(full.blip.mod[[j]], data = data)%*%psi.dwols[[j]];
        mod = glm(Y1_Y0~-1+Hpsi);
        psi.ce.dwols.star[[j]] = mod$coefficients;
        names(psi.ce.dwols.star[[j]]) = colnames(Hpsi);
        dopt.ce.dwols.star[[j]] = 1*(model.matrix(partial.blip.mod[[j]], data = data)%*%psi.ce.dwols.star[[j]] > 0);
        regret.ce.dwols.star[,j] = dopt.ce.dwols.star[[j]]*model.matrix(partial.blip.mod[[j]], data = data)%*%psi.ce.dwols.star[[j]] - 
                                     A[[j]]*model.matrix(full.blip.mod[[j]], data = data)%*%psi.dwols[[j]];
        obj$psi[[j]] = psi.ce.dwols.star[[j]];
        obj$opt.treat[[j]] = dopt.ce.dwols.star[[j]];
      }
      obj$regret = rowSums(as.matrix(regret.ce.dwols.star));
      obj$opt.Y = data$outcome + obj$regret;
    }

    ## CE g-est
    if(any(method %in% c("CE G-est"))){
      psi.ce.gest.star = NULL;
      dopt.ce.gest.star = NULL;
      regret.ce.gest.star = matrix(, nrow = n, ncol = K);
      for(j in K:1){
        if(j == K){
          Y = data$outcome;
        } else{
          Y = data$outcome + rowSums(as.matrix(regret.ce.gest.star[,(j+1):K]));
        }
        Hpsi = model.matrix(partial.blip.mod[[j]], data = data);
        Y1_Y0 = model.matrix(full.blip.mod[[j]], data = data)%*%psi.gest[[j]];
        mod = glm(Y1_Y0~-1+Hpsi);
        psi.ce.gest.star[[j]] = coef(mod);
        names(psi.ce.gest.star[[j]]) = colnames(Hpsi);
        dopt.ce.gest.star[[j]] = 1*(model.matrix(partial.blip.mod[[j]], data = data)%*%psi.ce.gest.star[[j]] > 0);
        regret.ce.gest.star[,j] = dopt.ce.gest.star[[j]]*model.matrix(partial.blip.mod[[j]], data = data)%*%psi.ce.gest.star[[j]] - 
                                     A[[j]]*model.matrix(full.blip.mod[[j]], data = data)%*%psi.gest[[j]];
        obj$psi[[j]] = psi.ce.gest.star[[j]];
        obj$opt.treat[[j]] = dopt.ce.gest.star[[j]];
      }
      obj$regret = rowSums(as.matrix(regret.ce.gest.star));
      obj$opt.Y = data$outcome + obj$regret;
    }
    
    return(obj);
  }

  obj = analyze.f(data);
  if(verbose == TRUE){
    t1 = proc.time();
    cat("Estimating the PATS took", (t1-t0)[3], "seconds.\n");
    eta = (t1-t0)[3]*B/60;
    if(var.estim != "none" & interrupt == TRUE & eta > 10){
      cat("Performing the bootstrap will take approximately", eta, "minutes.\n");
      cont <- readline("Continue? y/n: ");
      if(cont == "n" | cont == "no" | cont == "NO") stop("Aborted.")
    }
  }      
  
  if(var.estim == "none") return(obj);

  if(var.estim == "bootstrap"){
    t0 = proc.time();
    if(M == 0) M = n;
    psi.boot = list();
    for(b in 1:B){
      index = sample(1:n, replace = TRUE, size = M); 
      psi.boot[[b]] = analyze.f(data[index,])$psi;
      cont = "u";
      if(verbose == TRUE & b >= 10){
        eta = (B - b)*((proc.time() - t0)[3]/b);
        if(b > 50 & eta > 600 & cont != "y" & interrupt == TRUE){
          cont <- readline(paste("Estimated run time", 
                    round(eta/60), "minutes. Continue? y/n: "))
          if(cont == "n" | cont == "no" | cont == "NO"){
            stop("Aborted.\n")
          } else{
            cont = "y";
          }
        }
        if(b == 10){
          last = eta + 31;
        }
        if(eta > 30 & eta < (last - 30)){
          cat("Approximately", round(eta), "seconds remaining.\n");
          last = eta;
        }
      }
    }
    psi.boot = do.call(function(...) mapply(rbind, ..., SIMPLIFY = FALSE), psi.boot);
    # covmat = lapply(psi.boot, function(x){var(x)*M/n});
    covmat = lapply(psi.boot, function(x){var(x)});
    obj$covmat = covmat;
    obj$psi.boot = psi.boot;

    if(K == 2){      
      L = model.matrix(partial.blip.mod[[2]], data = data);
      LB = L%*%obj$psi[[2]];
      se = sqrt(diag(L%*%obj$covmat[[2]]%*%t(L)));
      ll = LB - 1.96*se;
      ul = LB + 1.96*se;
      obj$nonreg = mean(ll < 0 & ul > 1);
    }
  }


  if(var.estim == "adapt"){

    # Estimate p - non regularity parameter
    index = matrix(,nrow = B, ncol = n);
    psi.boot = list();
    for(b in 1:B){
      index[b,] = sample(1:n, replace = TRUE, size = n); 
      psi.boot[[b]] = analyze.f(data[index[b,],])$psi;
    }
    psi.boot = do.call(function(...) mapply(rbind, ..., SIMPLIFY = FALSE), psi.boot);
    covmat = lapply(psi.boot, var);

    L = model.matrix(partial.blip.mod[[2]], data = data);
    LB = L%*%obj$psi[[2]];
    se = sqrt(diag(L%*%covmat[[2]]%*%t(L)));
    ll = LB - 1.96*se;
    ul = LB + 1.96*se;
    obj$nonreg = mean(ll < 0 & ul > 0);
    t2 = proc.time();
    if(obj$nonreg == 0){
      obj$covmat = covmat;
      obj$psi.boot = psi.boot;
    } else{
      if(verbose == TRUE){
        cat("First stage bootstrap took", (t2-t1)[3], "seconds.\n");
        eta = (t2-t1)[3]*B2/60;
        if(interrupt == TRUE & eta > 10){
          cat("Testing each alpha value will take approximately", eta, "minutes.\n");
          cont = readline(paste("Continue? y/n: "))
          if(cont == "n" | cont == "no" | cont == "NO") stop("Aborted.")
        }
      }
      alpha = 0.025;
      cover = 0;
      psi.boot2 = list();
      while(alpha < 1 & cover < 0.95){
        CI = numeric(B);
        for(b in 1:B){
          L = model.matrix(partial.blip.mod[[2]], data = data[index[b,],]);
          LB = L%*%obj$psi[[2]];
          se = sqrt(diag(L%*%covmat[[2]]%*%t(L)));
          ll = LB - 1.96*se;
          ul = LB + 1.96*se;
          p = mean(ll < 0 & ul > 0);
          M = floor(n**((1+alpha*(1-p))/(1+alpha)));
          for(b2 in 1:B2){
            index2 = sample(index[b,], replace = TRUE, size = M); 
            psi.boot2[[b2]] = analyze.f(data[index2,])$psi[[1]];
          }
          psi.boot2 = do.call(function(...) mapply(rbind, ..., SIMPLIFY = FALSE), psi.boot2);
          temp.CI = numeric(length(partial.blip.mod[[1]]));
          for(i in 1:length(temp.CI)){
            temp.CI[i] = psi.boot[[1]][b, i] - quantile(sqrt(M)*(psi.boot2[[i]] - psi.boot[[1]][b, i]), 0.975)/sqrt(M) < obj$psi[[1]][i] &
                         psi.boot[[1]][b, i] - quantile(sqrt(M)*(psi.boot2[[i]] - psi.boot[[1]][b, i]), 0.025)/sqrt(M) > obj$psi[[1]][i]
          }
          CI[b] = mean(temp.CI);
        }
        cover = mean(CI);
        if(cover < 0.95 & verbose == TRUE){
          t3 = proc.time();
          cat("alpha =", alpha, ". cover =", cover, ". The optimal value of M was not found yet.\n", 
              (t3-t0)[3]/60, "minutes elapsed so far.\n");
          if(interrupt == TRUE){
            cont = readline(paste("Continue? y/n: "));
            if(cont == "n" | cont == "no" | cont == "NO") stop("Aborted.");
          }
        }
        if(cover < 0.95) alpha = alpha + 0.025;
      }
      M = floor(n**((1+alpha*(1-p))/(1+alpha)));
      obj$M = M;
      obj$alpha = alpha;
      for(b in 1:B){
        index = sample(1:n, replace = TRUE, size = M); 
        psi.boot[[b]] = analyze.f(data[index,])$psi;
      }
      psi.boot = do.call(function(...) mapply(rbind, ..., SIMPLIFY = FALSE), psi.boot);
      covmat = lapply(psi.boot, function(x){var(x)});
      obj$covmat = covmat;
      obj$psi.boot = psi.boot;
    }
  }
  return(obj);
}


### Example
#set.seed(74917491);

#expit = plogis;

# Psi parameters
#psi1 = c(-0.5, 1, -1);
#psi2 = c(-0.2, 0.2, -0.2);

#EX22.0 = 0.4608152; #E[X22|X21 = 0]
#EX22.1 = 0.5390204; #E[X22|X21 = 1]
#psi2.true = c(psi2[1] + EX22.0*psi2[3], psi2[2]+psi2[3]*(EX22.1 - EX22.0));

#psi1.true = c(-1.014, 1.027);

#n = 300;

# Stage 1 data generation
#X11 = rbinom(n, 1, 0.5);
#X12 = rbinom(n, 1, 0.5);
#X1 = cbind(1, X11, X12);
#A1 = rbinom(n, 1, expit(-1 + X11 + X12));
#A1opt = 1*(X1%*% psi1 > 0);
#mu1 = (A1opt - A1)*X1%*%psi1;

# Stage 2 data generation
#X21 = rbinom(n, 1, expit(-1 + X11 + A1));
#X22 = rbinom(n, 1, expit(-1 + X12 + A1));
#X2 = cbind(1, X21, X22);
#A2 = rbinom(n, 1, expit(-1 + X21 + X22));
#A2opt = 1*(X2%*% psi2 > 0);
#mu2 = (A2opt - A2)*X2%*%psi2;

# Outcome generation
#Y = X11 + X12 - mu1 - mu2 + rnorm(n);

#outcome = Y;
#partial.blip.mod = list(~X11, ~X21);
#full.blip.mod = list(~X11+X12,~X21+X22+A1+X11+X12);
#partial.treat.mod = list(A1~X11, A2~X21);
#full.treat.mod = list(A1~X11+X12, A2~X21+X22+A1+X12+X12);
#tf.mod = list(~X11+X12, ~X21+X22+A1+X11+X12);
#method = "ce dwols";
#weight = "default";
#var.estim = "adapt";
#B = 200;
#B2 = 200;
#M = 0;
#verbose = TRUE; 
#interrupt = TRUE;
#data = NULL;


#res = PATS(outcome = outcome, partial.blip.mod = partial.blip.mod,
#           full.blip.mod = full.blip.mod, partial.treat.mod = partial.treat.mod,
#           full.treat.mod = full.treat.mod, tf.mod = tf.mod, data = NULL,
#           method = "ce dwols", weight = "default", var.estim = "adapt",
#           B = 500, B2 = 500, M = 0, verbose = TRUE, interrupt = FALSE);



