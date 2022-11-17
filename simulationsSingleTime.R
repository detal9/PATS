# In this simulation, both X1 and X2 are confounders,
# but only X1 is used as a tailoring variable 

nrep = 1000;
expit = plogis;


### Scernario 1 - Correctly specified exposure model
###               Correctly specified outcome model

set.seed(37193179);
for(n in c(100, 1000, 10000)){

  ## Initialize objects
  coefs.dWOLS.naive = coefs.gest.naive =
  coefs.dWOLS.integrate = coefs.gest.integrate = 
  coefs.dWOLS.IPTW = coefs.gest.IPTW =
  coefs.dWOLS.ICE = coefs.gest.ICE = matrix(, nrow = nrep, ncol = 2);

  for(i in 1:nrep){

    ## Generate data
    X1 = rbinom(n, 1, p = 0.5);
    X2 = rbinom(n, 1, p = 0.5);
    A = rbinom(n, 1, p = expit(-0.5 + X1 + 0.5*X2));
    Y = rnorm(n, 0.25*X1 + X2 + A*(0.5 - 1*X1 + 1.5*X2));
     # True coefs = (0.5 + 1.5*E[X2] = 1.25, -1)  

  
    ## Propensity score models
    modA = glm(A ~ X1 + X2, family = "binomial");
    predA = predict(modA, type = "res");
    modA1 = glm(A ~ X1, family = "binomial");
    predA1 = predict(modA1, type = "res");


    ## Weights
    w = A*(predA1/predA) + (1 - A)*(1 - predA1)/(1 - predA);
    w.dWOLS = abs(A - predA1);
    w.dWOLS0 = abs(A - predA);


    ## X(X'X)^{-1}X'
    X = cbind(1, X1, X2);
    XXXX = X%*%solve((t(X)%*%X))%*%t(X);


    ## Naive dWOLS
    mod.dWOLS.naive = glm(Y ~ X1 + X2 + A*(1 + X1), weights = w.dWOLS0);
    coefs.dWOLS.naive[i,] = coef(mod.dWOLS.naive)[4:5];


    ## Naive g-est
    t.0 = rbind(sum(((A - predA)*Y - (A - predA)*XXXX%*%Y)),
                sum(((A*X1 - X1*predA)*Y - (A*X1 - X1*predA)*XXXX%*%Y)));

    t.psi0 = rbind(sum((-A*(A - predA) + (A - predA)*XXXX%*%A)),
                   sum((-A*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%A)));

    t.psi1 = rbind(sum((-A*X1*(A - predA) + (A - predA)*XXXX%*%(A*X1))),
                   sum((-A*X1*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%(A*X1))));
    coefs.gest.naive[i,] = solve(cbind(t.psi0, t.psi1), -t.0);


    ## Integrate dWOLS
    mod.dWOLS.integrate = glm(Y ~ X1 + X2 + A*(1 + X1 + X2), weights = w.dWOLS0);
    coefs.dWOLS.integrate[i,] = coef(mod.dWOLS.integrate)[4:5] +
                                c(mean(X2)*coef(mod.dWOLS.integrate)[6], 0);

    ## Integrate g-est
    t.0 = rbind(sum(((A - predA)*Y - (A - predA)*XXXX%*%Y)),
                sum(((A*X1 - X1*predA)*Y - (A*X1 - X1*predA)*XXXX%*%Y)),
                sum(((A*X2 - X2*predA)*Y - (A*X2 - X2*predA)*XXXX%*%Y)));

    t.psi0 = rbind(sum((-A*(A - predA) + (A - predA)*XXXX%*%A)),
                   sum((-A*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%A)),
                   sum((-A*(A*X2 - X2*predA) + (A*X2 - X2*predA)*XXXX%*%A)));

    t.psi1 = rbind(sum((-A*X1*(A - predA) + (A - predA)*XXXX%*%(A*X1))),
                   sum((-A*X1*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%(A*X1))),
                   sum((-A*X1*(A*X2 - X2*predA) + (A*X2 - X2*predA)*XXXX%*%(A*X1))));

    t.psi2 = rbind(sum((-A*X2*(A - predA) + (A - predA)*XXXX%*%(A*X2))),
                   sum((-A*X2*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%(A*X2))),
                   sum((-A*X2*(A*X2 - X2*predA) + (A*X2 - X2*predA)*XXXX%*%(A*X2))));

    gest = solve(cbind(t.psi0, t.psi1,t.psi2), -t.0);
    coefs.gest.integrate[i,] = c(gest[1] + mean(X2)*gest[3], gest[2]);


    ## IPTW + dWOLS
    mod.dWOLS.IPTW = glm(Y ~ X1 + X2 + A*(1 + X1), weights = w.dWOLS*w);
    coefs.dWOLS.IPTW[i,] = coef(mod.dWOLS.IPTW)[4:5];


    ## IPTW + gest
    t.0 = rbind(sum(w*((A - predA1)*Y - (A - predA1)*XXXX%*%Y)),
                sum(w*((A*X1 - X1*predA1)*Y - (A*X1 - X1*predA1)*XXXX%*%Y)));

    t.psi0 = rbind(sum(w*(-A*(A - predA1) + (A - predA1)*XXXX%*%A)),
                   sum(w*(-A*(A*X1 - X1*predA1) + (A*X1 - X1*predA1)*XXXX%*%A)));

    t.psi1 = rbind(sum(w*(-A*X1*(A - predA1) + (A - predA1)*XXXX%*%(A*X1))),
                   sum(w*(-A*X1*(A*X1 - X1*predA1) + (A*X1 - X1*predA1)*XXXX%*%(A*X1))));
    coefs.gest.IPTW[i,] = solve(cbind(t.psi0, t.psi1), -t.0);


    ## ICE dWOLS
    Y1_Y0.dWOLS = cbind(1, X1, X2)%*%coef(mod.dWOLS.integrate)[4:6];
    mod.dWOLS.ICE = glm(Y1_Y0.dWOLS~X1);
    coefs.dWOLS.ICE[i,] = coef(mod.dWOLS.ICE);


    ## ICE g-est
    Y1_Y0.gest = cbind(1, X1, X2)%*%gest;
    mod.gest.ICE = glm(Y1_Y0.gest~X1);
    coefs.gest.ICE[i,] = coef(mod.gest.ICE);
    

    print(data.frame(n, i, time = Sys.time()));
  }
  
  file.name = paste("results_scenario1n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
  save(coefs.dWOLS.naive, coefs.gest.naive,
       coefs.dWOLS.integrate, coefs.gest.integrate,
       coefs.dWOLS.IPTW, coefs.gest.IPTW,
       coefs.dWOLS.ICE, coefs.gest.ICE, file = file.name);
} 


### Scernario 2 - Incorrectly specified exposure model
###               Correctly specified outcome model

set.seed(37193179);
for(n in c(100, 1000, 10000)){

  ## Initialize objects
  coefs.dWOLS.naive = coefs.gest.naive =
  coefs.dWOLS.integrate = coefs.gest.integrate = 
  coefs.dWOLS.IPTW = coefs.gest.IPTW =
  coefs.dWOLS.ICE = coefs.gest.ICE = matrix(, nrow = nrep, ncol = 2);

  for(i in 1:nrep){

    ## Generate data
    X1 = rbinom(n, 1, p = 0.5);
    X2 = rbinom(n, 1, p = 0.5);
    A = rbinom(n, 1, p = expit(-0.5 + X1 + 0.5*X2 + X1*X2));
    Y = rnorm(n, 0.25*X1 + X2 + A*(0.5 - 1*X1 + 1.5*X2));
     # True coefs = (0.5 + 1.5*E[X2] = 1.25, -1)  

  
    ## Propensity score models
    modA = glm(A ~ X1 + X2, family = "binomial");
    predA = predict(modA, type = "res");
    modA1 = glm(A ~ X1, family = "binomial");
    predA1 = predict(modA1, type = "res");


    ## Weights
    w = A*(predA1/predA) + (1 - A)*(1 - predA1)/(1 - predA);
    w.dWOLS = abs(A - predA1);
    w.dWOLS0 = abs(A - predA);


    ## X(X'X)^{-1}X'
    X = cbind(1, X1, X2);
    XXXX = X%*%solve((t(X)%*%X))%*%t(X);


    ## Naive dWOLS
    mod.dWOLS.naive = glm(Y ~ X1 + X2 + A*(1 + X1), weights = w.dWOLS0);
    coefs.dWOLS.naive[i,] = coef(mod.dWOLS.naive)[4:5];


    ## Naive g-est
    t.0 = rbind(sum(((A - predA)*Y - (A - predA)*XXXX%*%Y)),
                sum(((A*X1 - X1*predA)*Y - (A*X1 - X1*predA)*XXXX%*%Y)));

    t.psi0 = rbind(sum((-A*(A - predA) + (A - predA)*XXXX%*%A)),
                   sum((-A*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%A)));

    t.psi1 = rbind(sum((-A*X1*(A - predA) + (A - predA)*XXXX%*%(A*X1))),
                   sum((-A*X1*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%(A*X1))));
    coefs.gest.naive[i,] = solve(cbind(t.psi0, t.psi1), -t.0);


    ## Integrate dWOLS
    mod.dWOLS.integrate = glm(Y ~ X1 + X2 + A*(1 + X1 + X2), weights = w.dWOLS0);
    coefs.dWOLS.integrate[i,] = coef(mod.dWOLS.integrate)[4:5] +
                                c(mean(X2)*coef(mod.dWOLS.integrate)[6], 0);

    ## Integrate g-est
    t.0 = rbind(sum(((A - predA)*Y - (A - predA)*XXXX%*%Y)),
                sum(((A*X1 - X1*predA)*Y - (A*X1 - X1*predA)*XXXX%*%Y)),
                sum(((A*X2 - X2*predA)*Y - (A*X2 - X2*predA)*XXXX%*%Y)));

    t.psi0 = rbind(sum((-A*(A - predA) + (A - predA)*XXXX%*%A)),
                   sum((-A*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%A)),
                   sum((-A*(A*X2 - X2*predA) + (A*X2 - X2*predA)*XXXX%*%A)));

    t.psi1 = rbind(sum((-A*X1*(A - predA) + (A - predA)*XXXX%*%(A*X1))),
                   sum((-A*X1*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%(A*X1))),
                   sum((-A*X1*(A*X2 - X2*predA) + (A*X2 - X2*predA)*XXXX%*%(A*X1))));

    t.psi2 = rbind(sum((-A*X2*(A - predA) + (A - predA)*XXXX%*%(A*X2))),
                   sum((-A*X2*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%(A*X2))),
                   sum((-A*X2*(A*X2 - X2*predA) + (A*X2 - X2*predA)*XXXX%*%(A*X2))));

    gest = solve(cbind(t.psi0, t.psi1,t.psi2), -t.0);
    coefs.gest.integrate[i,] = c(gest[1] + mean(X2)*gest[3], gest[2]);

  
    ## IPTW + dWOLS
    mod.dWOLS.IPTW = glm(Y ~ X1 + X2 + A*(1 + X1), weights = w.dWOLS*w);
    coefs.dWOLS.IPTW[i,] = coef(mod.dWOLS.IPTW)[4:5];


    ## IPTW + gest
    t.0 = rbind(sum(w*((A - predA1)*Y - (A - predA1)*XXXX%*%Y)),
                sum(w*((A*X1 - X1*predA1)*Y - (A*X1 - X1*predA1)*XXXX%*%Y)));

    t.psi0 = rbind(sum(w*(-A*(A - predA1) + (A - predA1)*XXXX%*%A)),
                   sum(w*(-A*(A*X1 - X1*predA1) + (A*X1 - X1*predA1)*XXXX%*%A)));

    t.psi1 = rbind(sum(w*(-A*X1*(A - predA1) + (A - predA1)*XXXX%*%(A*X1))),
                   sum(w*(-A*X1*(A*X1 - X1*predA1) + (A*X1 - X1*predA1)*XXXX%*%(A*X1))));
    coefs.gest.IPTW[i,] = solve(cbind(t.psi0, t.psi1), -t.0);


    ## ICE dWOLS
    Y1_Y0.dWOLS = cbind(1, X1, X2)%*%coef(mod.dWOLS.integrate)[4:6];
    mod.dWOLS.ICE = glm(Y1_Y0.dWOLS~X1);
    coefs.dWOLS.ICE[i,] = coef(mod.dWOLS.ICE);


    ## ICE g-est
    Y1_Y0.gest = cbind(1, X1, X2)%*%gest;
    mod.gest.ICE = glm(Y1_Y0.gest~X1);
    coefs.gest.ICE[i,] = coef(mod.gest.ICE);


    print(data.frame(n, i, time = Sys.time()));
  }
  
  file.name = paste("results_scenario2n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
  save(coefs.dWOLS.naive, coefs.gest.naive,
       coefs.dWOLS.integrate, coefs.gest.integrate,
       coefs.dWOLS.IPTW, coefs.gest.IPTW,
       coefs.dWOLS.ICE, coefs.gest.ICE, file = file.name);
} 


### Scernario 3 - Correctly specified exposure model
###               Incorrectly specified outcome model

set.seed(37193179);
for(n in c(100, 1000, 10000)){

  ## Initialize objects
  coefs.dWOLS.naive = coefs.gest.naive =
  coefs.dWOLS.integrate = coefs.gest.integrate = 
  coefs.dWOLS.IPTW = coefs.gest.IPTW =
  coefs.dWOLS.ICE = coefs.gest.ICE = matrix(, nrow = nrep, ncol = 2);

  for(i in 1:nrep){

    ## Generate data
    X1 = rbinom(n, 1, p = 0.5);
    X2 = rbinom(n, 1, p = 0.5);
    A = rbinom(n, 1, p = expit(-0.5 + X1 + 0.5*X2));
    Y = rnorm(n, 0.25*X1 + X2 + X1*X2 + A*(0.5 - 1*X1 + 1.5*X2));
     # True coefs = (0.5 + 1.5*E[X2] = 1.25, -1)  

  
    ## Propensity score models
    modA = glm(A ~ X1 + X2, family = "binomial");
    predA = predict(modA, type = "res");
    modA1 = glm(A ~ X1, family = "binomial");
    predA1 = predict(modA1, type = "res");


    ## Weights
    w = A*(predA1/predA) + (1 - A)*(1 - predA1)/(1 - predA);
    w.dWOLS = abs(A - predA1);
    w.dWOLS0 = abs(A - predA);


    ## X(X'X)^{-1}X'
    X = cbind(1, X1, X2);
    XXXX = X%*%solve((t(X)%*%X))%*%t(X);


    ## Naive dWOLS
    mod.dWOLS.naive = glm(Y ~ X1 + X2 + A*(1 + X1), weights = w.dWOLS0);
    coefs.dWOLS.naive[i,] = coef(mod.dWOLS.naive)[4:5];


    ## Naive g-est
    t.0 = rbind(sum(((A - predA)*Y - (A - predA)*XXXX%*%Y)),
                sum(((A*X1 - X1*predA)*Y - (A*X1 - X1*predA)*XXXX%*%Y)));

    t.psi0 = rbind(sum((-A*(A - predA) + (A - predA)*XXXX%*%A)),
                   sum((-A*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%A)));

    t.psi1 = rbind(sum((-A*X1*(A - predA) + (A - predA)*XXXX%*%(A*X1))),
                   sum((-A*X1*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%(A*X1))));
    coefs.gest.naive[i,] = solve(cbind(t.psi0, t.psi1), -t.0);


    ## Integrate dWOLS
    mod.dWOLS.integrate = glm(Y ~ X1 + X2 + A*(1 + X1 + X2), weights = w.dWOLS0);
    coefs.dWOLS.integrate[i,] = coef(mod.dWOLS.integrate)[4:5] +
                                c(mean(X2)*coef(mod.dWOLS.integrate)[6], 0);

    ## Integrate g-est
    t.0 = rbind(sum(((A - predA)*Y - (A - predA)*XXXX%*%Y)),
                sum(((A*X1 - X1*predA)*Y - (A*X1 - X1*predA)*XXXX%*%Y)),
                sum(((A*X2 - X2*predA)*Y - (A*X2 - X2*predA)*XXXX%*%Y)));

    t.psi0 = rbind(sum((-A*(A - predA) + (A - predA)*XXXX%*%A)),
                   sum((-A*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%A)),
                   sum((-A*(A*X2 - X2*predA) + (A*X2 - X2*predA)*XXXX%*%A)));

    t.psi1 = rbind(sum((-A*X1*(A - predA) + (A - predA)*XXXX%*%(A*X1))),
                   sum((-A*X1*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%(A*X1))),
                   sum((-A*X1*(A*X2 - X2*predA) + (A*X2 - X2*predA)*XXXX%*%(A*X1))));

    t.psi2 = rbind(sum((-A*X2*(A - predA) + (A - predA)*XXXX%*%(A*X2))),
                   sum((-A*X2*(A*X1 - X1*predA) + (A*X1 - X1*predA)*XXXX%*%(A*X2))),
                   sum((-A*X2*(A*X2 - X2*predA) + (A*X2 - X2*predA)*XXXX%*%(A*X2))));

    gest = solve(cbind(t.psi0, t.psi1,t.psi2), -t.0);
    coefs.gest.integrate[i,] = c(gest[1] + mean(X2)*gest[3], gest[2]);

  
    ## IPTW + dWOLS
    mod.dWOLS.IPTW = glm(Y ~ X1 + X2 + A*(1 + X1), weights = w.dWOLS*w);
    coefs.dWOLS.IPTW[i,] = coef(mod.dWOLS.IPTW)[4:5];


    ## IPTW + gest
    t.0 = rbind(sum(w*((A - predA1)*Y - (A - predA1)*XXXX%*%Y)),
                sum(w*((A*X1 - X1*predA1)*Y - (A*X1 - X1*predA1)*XXXX%*%Y)));

    t.psi0 = rbind(sum(w*(-A*(A - predA1) + (A - predA1)*XXXX%*%A)),
                   sum(w*(-A*(A*X1 - X1*predA1) + (A*X1 - X1*predA1)*XXXX%*%A)));

    t.psi1 = rbind(sum(w*(-A*X1*(A - predA1) + (A - predA1)*XXXX%*%(A*X1))),
                   sum(w*(-A*X1*(A*X1 - X1*predA1) + (A*X1 - X1*predA1)*XXXX%*%(A*X1))));
    coefs.gest.IPTW[i,] = solve(cbind(t.psi0, t.psi1), -t.0);


    ## ICE dWOLS
    Y1_Y0.dWOLS = cbind(1, X1, X2)%*%coef(mod.dWOLS.integrate)[4:6];
    mod.dWOLS.ICE = glm(Y1_Y0.dWOLS~X1);
    coefs.dWOLS.ICE[i,] = coef(mod.dWOLS.ICE);


    ## ICE g-est
    Y1_Y0.gest = cbind(1, X1, X2)%*%gest;
    mod.gest.ICE = glm(Y1_Y0.gest~X1);
    coefs.gest.ICE[i,] = coef(mod.gest.ICE);


    print(data.frame(n, i, time = Sys.time()));
  }
  
  file.name = paste("results_scenario3n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
  save(coefs.dWOLS.naive, coefs.gest.naive,
       coefs.dWOLS.integrate, coefs.gest.integrate,
       coefs.dWOLS.IPTW, coefs.gest.IPTW,
       coefs.dWOLS.ICE, coefs.gest.ICE, file = file.name);
} 



