library(doParallel)
ncores <- detectCores()
ncores
detectCores(logical = FALSE)
cl <- makePSOCKcluster(ncores/2)
registerDoParallel(cl)

### Scernario 1 - Correctly specified exposure model
###               Correctly specified outcome model


# In this simulation, both X1 and X2 are confounders,
# but only X1 is used as a tailoring variable 

nrep = 1000;
expit = plogis;
source("PATS_functions_v2.R");



# Psi parameters
psi1 = c(-0.5, 1, -1);
psi2 = c(-0.5, 1, -1);

### Scernario 1 - Correctly specified exposure model
###               Correctly specified outcome model
#
#    # Stage 1 data generation
#    X11 = rbinom(n, 1, 0.5);
#    X12 = rbinom(n, 1, 0.5);
#    X1 = cbind(1, X11, X12);
#    A1 = rbinom(n, 1, expit(-1 + X11 + X12));
#    A1opt = 1*(X1%*% psi1 > 0);
#    mu1 = (A1opt - A1)*X1%*%psi1;
#
#    # Stage 2 data generation
#    X21 = rbinom(n, 1, expit(-1 + X11 + A1));
#    X22 = rbinom(n, 1, expit(-1 + X12 + A1));
#    X2 = cbind(1, X21, X22);
#    A2 = rbinom(n, 1, expit(-1 + X21 + X22));
#    A2opt = 1*(X2%*% psi2 > 0);
# H^C_2 = {X_22, X_11, X_12, A1}
# H*_2 = X_21
#    mu2 = (A2opt - A2)*X2%*%psi2;
#
# gamma2* = \int{X_22, X_11, X_12, A1} A(psi02 + psi12*X_21 + psi22*X_22) dF_{H^C_2|H*_2}  
#         = A(psi02 + psi12*X_21 + psi22*\int{X_22, X_11, X_12, A1} X_22 dF_{H^C_2|H*_2}
#         = A(pis02 + psi12*X_21 + psi22*E[X_22|X_21])
#

# Monte Carlo estimation of E[X22|X21]
# set.seed(37193179);
# n = 20000000;
# X11 = rbinom(n, 1, 0.5);
# X12 = rbinom(n, 1, 0.5);
# A1 = rbinom(n, 1, expit(-1 + X11 + X12));
# X21 = rbinom(n, 1, expit(-1 + X11 + A1));
# X22 = rbinom(n, 1, expit(-1 + X12 + A1));
EX22.0 = 0.4608152; #E[X22|X21 = 0]
EX22.1 = 0.5390204; #E[X22|X21 = 1]
psi2.true = c(psi2[1] + EX22.0*psi2[3], psi2[2]+psi2[3]*(EX22.1 - EX22.0));

# Monte Carlo estimation of gamma1*
# set.seed(37193179);
# n = 1000000;
# X11 = rbinom(n, 1, 0.5);
# X12 = rbinom(n, 1, 0.5);
# A1.0 = rep(0, n);
# A1.1 = rep(1, n);
# A1opt = 1*(cbind(1, X11, X12)%*% psi1 > 0);
# mu1.0 = (A1opt - A1.0)*cbind(1, X11, X12)%*%psi1;
# mu1.1 = (A1opt - A1.1)*cbind(1, X11, X12)%*%psi1;
# X21.0 = rbinom(n, 1, expit(-1 + X11 + 0));
# X22.0 = rbinom(n, 1, expit(-1 + X12 + 0));
# X21.1 = rbinom(n, 1, expit(-1 + X11 + 1));
# X22.1 = rbinom(n, 1, expit(-1 + X12 + 1));
# A2.0 = 1*(cbind(1, X21.0)%*%psi2.true > 0);
# A2.1 = 1*(cbind(1, X21.1)%*%psi2.true > 0);
# A2opt.0 = 1*(cbind(1, X21.0, X22.0)%*% psi2 > 0);
# A2opt.1 = 1*(cbind(1, X21.1, X22.1)%*% psi2 > 0);
# mu2.0 = (A2opt.0 - A2.0)*cbind(1, X21.0, X22.0)%*%psi2;
# mu2.1 = (A2opt.1 - A2.1)*cbind(1, X21.1, X22.1)%*%psi2;
# Y.0 = X11 + X12 - mu1.0 - mu2.0;
# Y.1 = X11 + X12 - mu1.1 - mu2.1;
# glm(Y.1 - Y.0 ~ X11);
psi1.true = c(-1.014, 1.027);


file.name = paste("results_CIPartTwoC",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");

set.seed(841);
n = 300;
{

  ## Initialize objects
  results.sim = list();

  for(i in 1:nrep){

    # Stage 1 data generation
    X11 = rbinom(n, 1, 0.5);
    X12 = rbinom(n, 1, 0.5);
    X1 = cbind(1, X11, X12);
    A1 = rbinom(n, 1, expit(-1 + X11 + X12));
    A1opt = 1*(X1%*% psi1 > 0);
    mu1 = (A1opt - A1)*X1%*%psi1;

    # Stage 2 data generation
    X21 = rbinom(n, 1, expit(-1 + X11 + A1));
    X22 = rbinom(n, 1, expit(-1 + X12 + A1));
    X2 = cbind(1, X21, X22);
    A2 = rbinom(n, 1, expit(-1 + X21 + X22));
    A2opt = 1*(X2%*% psi2 > 0);
    mu2 = (A2opt - A2)*X2%*%psi2;

    # Outcome generation
    Y = X11 + X12 - mu1 - mu2 + rnorm(n);

    results.sim[[i]] = PATS(outcome = Y,
                            partial.blip.mod = list(~X11, ~X21),
                            full.blip.mod = list(~X11+X12, ~X21 + X22),
                            full.treat.mod = list(A1~X11+X12, A2~X21+X22+A1+X12+X12),
                            tf.mod = list(~X11+X12, ~X21+X22+A1+X11+X12),
                            method = "ce dwols", weight = "default", var.estim = "adapt",
                            B = 500, B2 = 500, verbose = TRUE);

    print(data.frame(n, i, time = Sys.time()));
    save(results.sim, file = file.name);
  }
} 



