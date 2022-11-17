setwd("C:\\Users\\denis\\Dropbox\\Travail\\Recherche\\MSNMM\\Simulations\\Results");

nrep = 1000;
expit = plogis;
require(survival);
require(rms);


### Scernario 1 - Non-informative censoring
set.seed(37193179);
for(n in c(5000)){

  ## Initialize objects
  coefs.dWOLS.naive = coefs.dWOLS.IPTW = coefs.dWOLS.ICE = matrix(, nrow = nrep, ncol = 3);
  opt.dWOLS.naive = opt.dWOLS.IPTW = opt.dWOLS.ICE = opt = matrix(, nrow = n, ncol = nrep);
  val.dWOLS.naive = val.dWOLS.IPTW = val.dWOLS.ICE = val = matrix(, nrow = n, ncol = nrep);

  coefs.dWOLS.naive2 = coefs.dWOLS.IPTW2 = coefs.dWOLS.ICE2 = matrix(, nrow = nrep, ncol = 3);
  opt.dWOLS.naive2 = opt.dWOLS.IPTW2 = opt.dWOLS.ICE2 = matrix(, nrow = n, ncol = nrep);
  val.dWOLS.naive2 = val.dWOLS.IPTW2 = val.dWOLS.ICE2 = matrix(, nrow = n, ncol = nrep);

  for(i in 1:nrep){

    ## Generate data
    ER = rbinom(n, 1, p = 0.75);
    PR = rbinom(n, 1, p = 0.05 + 0.5*ER);
    AGE = rnorm(n, mean = 0, sd = 1);
    BMI = rbinom(n, 1, p = expit(-0.4 + AGE));
    YEAR = runif(n, min = 0, max = 25);
    A = rbinom(n, 1, p = expit(-3 + 0.5*AGE + 0.1*YEAR + 3*ER + 0.5*PR)); 
    Y = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
        A*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI) + rnorm(n, mean = 0, sd = 0.5);
    TCensor = 3.5 - 0.05*YEAR + 0.1*AGE + 0.2*A + rnorm(n, mean = 0, sd = 0.5);   
     # blip = -0.2 + 0.3*ER + 0.3*PR - 0.3*BMI
     # lm(blip ~ ER + BMI);
    Yobs = pmin(Y, TCensor);
    Censor = Y > TCensor;
    opt[,i] = 1*(-0.185 + 0.45*ER - 0.3*BMI) > 0;
    val[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
              opt[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);


    ## Propensity score models
    modA = glm(A ~ ER + PR + BMI + AGE + YEAR, family = "binomial");
    predA = predict(modA, type = "res");
    modA1 = glm(A ~ ER + BMI, family = "binomial");
    predA1 = predict(modA1, type = "res");
    modC = psm(Surv(Yobs, Censor) ~ A + ER + PR + BMI + AGE + YEAR, dist = "gaussian");
    predC = numeric(n);
    for(j in (1:n)[Censor == 0]){
      dati = data.frame(A = A[j], ER = ER[j], PR = PR[j], BMI = BMI[j],
                        AGE = AGE[j], YEAR = YEAR[j]);
      predC[j] = 1 - survest(modC, newdata = dati, times = Yobs[j], se.fit = FALSE)$surv;
    }

    ## Weights
    w = A*(predA1/predA) + (1 - A)*(1 - predA1)/(1 - predA);
    w.dWOLS = abs(A - predA1);
    w.dWOLS0 = abs(A - predA);
    wC = 1/(1 - predC);

    ## Naive dWOLS
    mod.dWOLS.naive = lm(Yobs ~ ER + PR + BMI + AGE + YEAR + A*(ER + BMI),
                         subset = Censor == 0,
                         weights = w.dWOLS0*wC);
    coefs.dWOLS.naive[i,] = tail(coef(mod.dWOLS.naive),3);
    opt.dWOLS.naive[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.naive[i,] > 0);
    val.dWOLS.naive[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                          opt.dWOLS.naive[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);

    mod.dWOLS.naive = survreg(Surv(exp(Yobs), 1 - Censor) ~ ER + PR + BMI + AGE + YEAR + A*(ER + BMI),
                              dist = "lognormal",
                              weights = w.dWOLS0);
    coefs.dWOLS.naive2[i,] = tail(coef(mod.dWOLS.naive),3);
    opt.dWOLS.naive2[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.naive2[i,] > 0);
    val.dWOLS.naive2[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                          opt.dWOLS.naive2[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);


    ## IPTW + dWOLS
    mod.dWOLS.IPTW = lm(Yobs ~ ER + PR + BMI + AGE + YEAR + A*(ER + BMI),
                        subset = Censor == 0,
                        weights = w.dWOLS*w*wC);
    coefs.dWOLS.IPTW[i,] = tail(coef(mod.dWOLS.IPTW), 3);
    opt.dWOLS.IPTW[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.IPTW[i,] > 0);
    val.dWOLS.IPTW[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                         opt.dWOLS.IPTW[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);

    mod.dWOLS.IPTW = survreg(Surv(exp(Yobs), 1 - Censor) ~ ER + PR + BMI + AGE + YEAR + A*(ER + BMI),
                             dist = "lognormal",
                             weights = w.dWOLS*w);
    coefs.dWOLS.IPTW2[i,] = tail(coef(mod.dWOLS.IPTW), 3);
    opt.dWOLS.IPTW2[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.IPTW2[i,] > 0);
    val.dWOLS.IPTW2[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                         opt.dWOLS.IPTW2[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);


    ## ICE dWOLS
    mod.dWOLS.integrate = lm(Yobs ~ ER + PR + BMI + AGE + YEAR + A*(ER + PR + BMI),
                             subset = Censor == 0,
                             weights = w.dWOLS0*wC);
    Y1_Y0.dWOLS = cbind(1, ER, PR, BMI)%*%tail(coef(mod.dWOLS.integrate), 4);
    mod.dWOLS.ICE = glm(Y1_Y0.dWOLS ~ ER + BMI);
    coefs.dWOLS.ICE[i,] = coef(mod.dWOLS.ICE);
    opt.dWOLS.ICE[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.ICE[i,] > 0);
    val.dWOLS.ICE[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                        opt.dWOLS.ICE[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);

    mod.dWOLS.integrate = survreg(Surv(exp(Yobs), 1 - Censor) ~ ER + PR + BMI + AGE + YEAR + A*(ER + PR + BMI),
                                  dist = "lognormal",
                                  weights = w.dWOLS0);
    Y1_Y0.dWOLS = cbind(1, ER, PR, BMI)%*%tail(coef(mod.dWOLS.integrate), 4);
    mod.dWOLS.ICE = glm(Y1_Y0.dWOLS ~ ER + BMI);
    coefs.dWOLS.ICE2[i,] = coef(mod.dWOLS.ICE);
    opt.dWOLS.ICE2[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.ICE2[i,] > 0);
    val.dWOLS.ICE2[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                        opt.dWOLS.ICE2[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);


    print(data.frame(n, i, time = Sys.time()));
  }
  
  file.name = paste("results_scenarioC1","_",gsub("-","",Sys.Date()),".Rdata", sep="");
  save(coefs.dWOLS.naive, coefs.dWOLS.IPTW, coefs.dWOLS.ICE,
       coefs.dWOLS.naive2, coefs.dWOLS.IPTW2, coefs.dWOLS.ICE2,
       val, opt, file = file.name);
} 







### Scernario 2 - Informative censoring
set.seed(37193179);
for(n in c(5000)){

  ## Initialize objects
  coefs.dWOLS.naive = coefs.dWOLS.IPTW = coefs.dWOLS.ICE = matrix(, nrow = nrep, ncol = 3);
  opt.dWOLS.naive = opt.dWOLS.IPTW = opt.dWOLS.ICE = opt = matrix(, nrow = n, ncol = nrep);
  val.dWOLS.naive = val.dWOLS.IPTW = val.dWOLS.ICE = val = matrix(, nrow = n, ncol = nrep);

  for(i in 1:nrep){

    ## Generate data
    ER = rbinom(n, 1, p = 0.75);
    PR = rbinom(n, 1, p = 0.05 + 0.5*ER);
    AGE = rnorm(n, mean = 0, sd = 1);
    BMI = rbinom(n, 1, p = expit(-0.4 + AGE));
    YEAR = runif(n, min = 0, max = 25);
    A = rbinom(n, 1, p = expit(-3 + 0.5*AGE + 0.1*YEAR + 3*ER + 0.5*PR)); 
    Y = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
        A*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI) + rnorm(n, mean = 0, sd = 0.5);
    Censor = 1*(exp(Y) + YEAR > 27);   
     # blip = -0.2 + 0.3*ER + 0.3*PR - 0.3*BMI
     # lm(blip ~ ER + BMI);
    Yobs = ifelse(Censor == 0, Y, log(27 - YEAR));
    opt[,i] = 1*(-0.185 + 0.45*ER - 0.3*BMI) > 0;
    val[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
              opt[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);


    ## Propensity score models
    modA = glm(A ~ ER + PR + BMI + AGE + YEAR, family = "binomial");
    predA = predict(modA, type = "res");
    modA1 = glm(A ~ ER + BMI, family = "binomial");
    predA1 = predict(modA1, type = "res");
    modC = cph(Surv(Yobs, Censor) ~ A + ER + PR + BMI + AGE + YEAR, surv = TRUE);
    modC = psm(Surv(Yobs, Censor) ~ A + ER + PR + BMI + AGE + YEAR, dist = "gaussian");
    predC = numeric(n);
    for(j in (1:n)[Censor == 0]){
      dati = data.frame(A = A[j], ER = ER[j], PR = PR[j], BMI = BMI[j],
                        AGE = AGE[j], YEAR = YEAR[j]);
      predC[j] = 1 - survest(modC, newdata = dati, times = Yobs[j], se.fit = FALSE)$surv;
    }

    ## Weights
    w = A*(predA1/predA) + (1 - A)*(1 - predA1)/(1 - predA);
    w.dWOLS = abs(A - predA1);
    w.dWOLS0 = abs(A - predA);
    wC = 1/(1 - predC);

    ## Naive dWOLS
    mod.dWOLS.naive = lm(Yobs ~ ER + PR + BMI + AGE + YEAR + A*(ER + BMI),
                         subset = Censor == 0,
                         weights = w.dWOLS0*wC);
    coefs.dWOLS.naive[i,] = tail(coef(mod.dWOLS.naive),3);
    opt.dWOLS.naive[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.naive[i,] > 0);
    val.dWOLS.naive[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                          opt.dWOLS.naive[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);

    mod.dWOLS.naive = survreg(Surv(exp(Yobs), 1 - Censor) ~ ER + PR + BMI + AGE + YEAR + A*(ER + BMI),
                              dist = "lognormal",
                              weights = w.dWOLS0);
    coefs.dWOLS.naive2[i,] = tail(coef(mod.dWOLS.naive),3);
    opt.dWOLS.naive2[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.naive2[i,] > 0);
    val.dWOLS.naive2[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                          opt.dWOLS.naive2[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);


    ## IPTW + dWOLS
    mod.dWOLS.IPTW = lm(Yobs ~ ER + PR + BMI + AGE + YEAR + A*(ER + BMI),
                        subset = Censor == 0,
                        weights = w.dWOLS*w*wC);
    coefs.dWOLS.IPTW[i,] = tail(coef(mod.dWOLS.IPTW), 3);
    opt.dWOLS.IPTW[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.IPTW[i,] > 0);
    val.dWOLS.IPTW[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                         opt.dWOLS.IPTW[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);

    mod.dWOLS.IPTW = survreg(Surv(exp(Yobs), 1 - Censor) ~ ER + PR + BMI + AGE + YEAR + A*(ER + BMI),
                             dist = "lognormal",
                             weights = w.dWOLS*w);
    coefs.dWOLS.IPTW2[i,] = tail(coef(mod.dWOLS.IPTW), 3);
    opt.dWOLS.IPTW2[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.IPTW2[i,] > 0);
    val.dWOLS.IPTW2[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                         opt.dWOLS.IPTW2[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);


    ## ICE dWOLS
    mod.dWOLS.integrate = lm(Yobs ~ ER + PR + BMI + AGE + YEAR + A*(ER + PR + BMI),
                             subset = Censor == 0,
                             weights = w.dWOLS0*wC);
    Y1_Y0.dWOLS = cbind(1, ER, PR, BMI)%*%tail(coef(mod.dWOLS.integrate), 4);
    mod.dWOLS.ICE = glm(Y1_Y0.dWOLS ~ ER + BMI);
    coefs.dWOLS.ICE[i,] = coef(mod.dWOLS.ICE);
    opt.dWOLS.ICE[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.ICE[i,] > 0);
    val.dWOLS.ICE[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                        opt.dWOLS.ICE[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);

    mod.dWOLS.integrate = survreg(Surv(exp(Yobs), 1 - Censor) ~ ER + PR + BMI + AGE + YEAR + A*(ER + PR + BMI),
                                  dist = "lognormal",
                                  weights = w.dWOLS0);
    Y1_Y0.dWOLS = cbind(1, ER, PR, BMI)%*%tail(coef(mod.dWOLS.integrate), 4);
    mod.dWOLS.ICE = glm(Y1_Y0.dWOLS ~ ER + BMI);
    coefs.dWOLS.ICE2[i,] = coef(mod.dWOLS.ICE);
    opt.dWOLS.ICE2[,i] = 1*(cbind(1, ER, BMI)%*%coefs.dWOLS.ICE2[i,] > 0);
    val.dWOLS.ICE2[,i] = 2.5 - 0.05*AGE - 0.1*ER - 0.1*PR - 0.1*BMI + 0.1*YEAR +
                        opt.dWOLS.ICE2[,i]*(-0.2 + 0.3*ER + 0.3*PR - 0.3*BMI);

    print(data.frame(n, i, time = Sys.time()));
  }
  
  file.name = paste("results_scenarioC2","_",gsub("-","",Sys.Date()),".Rdata", sep="");
  save(coefs.dWOLS.naive, coefs.dWOLS.IPTW, coefs.dWOLS.ICE,
       coefs.dWOLS.naive2, coefs.dWOLS.IPTW2, coefs.dWOLS.ICE2,
       val, opt, file = file.name);
} 
