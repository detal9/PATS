### Scernario 1 - Correctly specified exposure model
###               Correctly specified outcome model


# In this simulation, both X1 and X2 are confounders,
# but only X1 is used as a tailoring variable 

nrep = 1000;
expit = plogis;

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



set.seed(37193179);
for(n in c(100,1000,10000)){

  ## Initialize objects
  coefs.dWOLS.naive = coefs.gest.naive =
  coefs.dWOLS.integrate = coefs.gest.integrate = 
  coefs.dWOLS.IPTW = coefs.gest.IPTW =
  coefs.dWOLS.ICE = coefs.gest.ICE = matrix(, nrow = nrep, ncol = 4);

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

  
    ## Propensity score models
    modA1 = glm(A1 ~ X11 + X12, family = "binomial");
    predA1 = predict(modA1, type = "res");

    modA2 = glm(A2 ~ X21 + X22 + X11 + X12 + A1, family = "binomial");
    predA2 = predict(modA2, type = "res");

    modA1S = glm(A1 ~ X11, family = "binomial");
    predA1S = predict(modA1S, type = "res");

    modA2S = glm(A2 ~ X21, family = "binomial");
    predA2S = predict(modA2S, type = "res");


    ## Weights
    w2 = (A2*(predA2S/predA2) + (1 - A2)*(1 - predA2S)/(1 - predA2));
    w1 = A1*(predA1S/predA1) + (1 - A1)*(1 - predA1S)/(1 - predA1);
    w1.dWOLS = abs(A1 - predA1S);
    w1.dWOLS0 = abs(A1 - predA1);
    w2.dWOLS = abs(A2 - predA2S);
    w2.dWOLS0 = abs(A2 - predA2);


    ## Naive dWOLS
    # Stage 2
    m2 = summary(glm(Y~ X11 + X12 + X21 + X22 + A1*(1 + X11 + X12 + X21 + X21) 
              + A2*(1 + X21), weights = w2.dWOLS0, family = "gaussian"))$coef;
    psi2.hat = m2[c("A2", "X21:A2"),1];

    # Stage 1
    H2 = cbind(1, X21);
    A2opt.hat = 1*(H2%*%psi2.hat > 0);
    mu2.hat = (A2opt.hat - A2)*H2%*%psi2.hat;
    Yopt2.hat = Y + mu2.hat;
    m1 = summary(glm(Yopt2.hat ~ X11 + X12 
                     + A1*(1 + X11), weights = w1.dWOLS0))$coef;
    coefs.dWOLS.naive[i,] = c(m1[c("A1", "X11:A1"),1], psi2.hat);


    ## Naive g-est

    # Stage 2

    # X(X'X)^{-1}X'
    H2 = cbind(1, X11, X12, X21, X22, A1, A1*X11, A1*X12, A1*X21, A1*X22);
    XXXX = H2%*%solve(t(H2)%*%H2)%*%t(H2);
    XXXXY = XXXX%*%Y;
    XXXXA2 = XXXX%*%A2;
    XXXXA2X21 = XXXX%*%(A2*X21);

    t.0 = rbind(sum(((A2 - predA2)*Y - (A2 - predA2)*XXXXY)),
                sum(((A2*X21 - X21*predA2)*Y - (A2*X21 - X21*predA2)*XXXXY)));

    t.psi0 = rbind(sum((-A2*(A2 - predA2) + (A2 - predA2)*XXXXA2)),
                   sum((-A2*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2)));

    t.psi1 = rbind(sum((-A2*X21*(A2 - predA2) + (A2 - predA2)*XXXXA2X21)),
               sum((-A2*X21*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2X21)));

    psi2.hat = solve(cbind(t.psi0, t.psi1), -t.0);

    # stage 1
    H2 = cbind(1, X21);
    A2opt.hat = 1*(H2%*%psi2.hat > 0);
    mu2.hat = (A2opt.hat - A2)*H2%*%psi2.hat;
    Yopt2.hat = Y + mu2.hat;

    H1 = cbind(1, X11, X12);
    XXXX = H1%*%solve(t(H1)%*%H1)%*%t(H1);
    XXXXY = XXXX%*%Yopt2.hat;
    XXXXA1 = XXXX%*%A1;
    XXXXA1X11 = XXXX%*%(A1*X11);

    t.0 = rbind(sum(((A1 - predA1)*Yopt2.hat - (A1 - predA1)*XXXXY)),
                sum(((A1*X11 - X11*predA1)*Yopt2.hat - (A1*X11 - X11*predA1)*XXXXY)));

    t.psi0 = rbind(sum((-A1*(A1 - predA1) + (A1 - predA1)*XXXXA1)),
                   sum((-A1*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1)));

    t.psi1 = rbind(sum((-A1*X11*(A1 - predA1) + (A1 - predA1)*XXXXA1X11)),
                   sum((-A1*X11*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1X11)));

    psi1.hat = solve(cbind(t.psi0, t.psi1), -t.0);
    coefs.gest.naive[i,] = c(psi1.hat, psi2.hat);


    ## Integrate dWOLS
    # Stage 2
    m2 = summary(glm(Y~ X11 + X12 + X21 + X22 + A1*(1 + X11 + X12 + X21 + X21) 
              + A2*(1 + X21 + X22), weights = w2.dWOLS0, family = "gaussian"))$coef;
    psi2.hat = m2[c("A2", "X21:A2", "X22:A2"),1];
    psi2.hat.dWOLS.integrate = psi2.hat;
    EX12 = prop.table(table(X12, X11), 2)[2,];
    EX22 = prop.table(table(X22, X21), 2)[2,];
    psi2.star.hat = c(psi2.hat[1] + psi2.hat[3]*EX22[1],
                      psi2.hat[2] + psi2.hat[3]*(EX22[2] - EX22[1]));

    # Stage 1
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    m1 = summary(glm(Yopt2.hat ~ X11 + X12 
                     + A1*(1 + X11 + X12), weights = w1.dWOLS0))$coef;
    psi1.hat = m1[c("A1", "X11:A1", "X12:A1"),1];
    psi1.hat.dWOLS.integrate = psi1.hat;
    psi1.star.hat = c(psi1.hat[1] + psi1.hat[3]*EX12[1],
                      psi1.hat[2] + psi1.hat[3]*(EX12[2] - EX12[1]));

    coefs.dWOLS.integrate[i,] = c(psi1.star.hat, psi2.star.hat);


    ## Integrate g-est
    # Stage 2
    # X(X'X)^{-1}X'
    H2 = cbind(1, X11, X12, X21, X22, A1, A1*X11, A1*X12, A1*X21, A1*X22);
    XXXX = H2%*%solve(t(H2)%*%H2)%*%t(H2);
    XXXXY = XXXX%*%Y;
    XXXXA2 = XXXX%*%A2;
    XXXXA2X21 = XXXX%*%(A2*X21);
    XXXXA2X22 = XXXX%*%(A2*X22);

    t.0 = rbind(sum(((A2 - predA2)*Y - (A2 - predA2)*XXXXY)),
                sum(((A2*X21 - X21*predA2)*Y - (A2*X21 - X21*predA2)*XXXXY)),
                sum(((A2*X22 - X22*predA2)*Y - (A2*X22 - X22*predA2)*XXXXY)));

    t.psi0 = rbind(sum((-A2*(A2 - predA2) + (A2 - predA2)*XXXXA2)),
                   sum((-A2*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2)),
                   sum((-A2*(A2*X22 - X22*predA2) + (A2*X22 - X22*predA2)*XXXXA2)));

    t.psi1 = rbind(sum((-A2*X21*(A2 - predA2) + (A2 - predA2)*XXXXA2X21)),
               sum((-A2*X21*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2X21)),
               sum((-A2*X21*(A2*X22 - X22*predA2) + (A2*X22 - X22*predA2)*XXXXA2X21)));

    t.psi2 = rbind(sum((-A2*X22*(A2 - predA2) + (A2 - predA2)*XXXXA2X22)),
               sum((-A2*X22*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2X22)),
               sum((-A2*X22*(A2*X22 - X22*predA2) + (A2*X22 - X22*predA2)*XXXXA2X22)));

    psi2.hat = solve(cbind(t.psi0, t.psi1, t.psi2), -t.0);
    psi2.hat.gest.integrate = psi2.hat;
    psi2.star.hat = c(psi2.hat[1] + psi2.hat[3]*EX22[1],
                      psi2.hat[2] + psi2.hat[3]*(EX22[2] - EX22[1]));

    # stage 1
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    H1 = cbind(1, X11, X12);
    XXXX = H1%*%solve(t(H1)%*%H1)%*%t(H1);
    XXXXY = XXXX%*%Yopt2.hat;
    XXXXA1 = XXXX%*%A1;
    XXXXA1X11 = XXXX%*%(A1*X11);
    XXXXA1X12 = XXXX%*%(A1*X12);

    t.0 = rbind(sum(((A1 - predA1)*Yopt2.hat - (A1 - predA1)*XXXXY)),
                sum(((A1*X11 - X11*predA1)*Yopt2.hat - (A1*X11 - X11*predA1)*XXXXY)),
                sum(((A1*X12 - X12*predA1)*Yopt2.hat - (A1*X12 - X12*predA1)*XXXXY)));

    t.psi0 = rbind(sum((-A1*(A1 - predA1) + (A1 - predA1)*XXXXA1)),
                   sum((-A1*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1)),
                   sum((-A1*(A1*X12 - X12*predA1) + (A1*X12 - X12*predA1)*XXXXA1)));

    t.psi1 = rbind(sum((-A1*X11*(A1 - predA1) + (A1 - predA1)*XXXXA1X11)),
                   sum((-A1*X11*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1X11)),
                   sum((-A1*X11*(A1*X12 - X12*predA1) + (A1*X12 - X12*predA1)*XXXXA1X11)));

    t.psi2 = rbind(sum((-A1*X12*(A1 - predA1) + (A1 - predA1)*XXXXA1X12)),
                   sum((-A1*X12*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1X12)),
                   sum((-A1*X12*(A1*X12 - X12*predA1) + (A1*X12 - X12*predA1)*XXXXA1X12)));

    psi1.hat = solve(cbind(t.psi0, t.psi1, t.psi2), -t.0);
    psi1.hat.gest.integrate = psi1.hat;
    psi1.star.hat = c(psi1.hat[1] + psi1.hat[3]*EX12[1],
                      psi1.hat[2] + psi1.hat[3]*(EX12[2] - EX12[1]));
    coefs.gest.integrate[i,] =  c(psi1.star.hat, psi2.star.hat);


    ## IPTW + dWOLS
    # Stage 2
    m2 = summary(glm(Y~ X11 + X12 + X21 + X22 + A1*(1 + X11 + X12 + X21 + X21) 
              + A2*(1 + X21), weights = w2.dWOLS*w2, family = "gaussian"))$coef;
    psi2.star.hat = m2[c("A2", "X21:A2"),1];

    # Stage 1
    psi2.hat = psi2.hat.dWOLS.integrate;
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    m1 = summary(glm(Yopt2.hat ~ X11 + X12 
                     + A1*(1 + X11), weights = w1*w1.dWOLS))$coef;
    coefs.dWOLS.IPTW[i,] = c(m1[c("A1", "X11:A1"),1], psi2.star.hat);


    ## IPTW + gest
    # Stage 2

    # X(X'X)^{-1}X'
    H2 = cbind(1, X11, X12, X21, X22, A1, A1*X11, A1*X12, A1*X21, A1*X22);
    XXXX = H2%*%solve(t(H2)%*%H2)%*%t(H2);
    XXXXY = XXXX%*%Y;
    XXXXA2 = XXXX%*%A2;
    XXXXA2X21 = XXXX%*%(A2*X21);

    t.0 = rbind(sum(w2*((A2 - predA2S)*Y - (A2 - predA2S)*XXXXY)),
                sum(w2*((A2*X21 - X21*predA2S)*Y - (A2*X21 - X21*predA2S)*XXXXY)));

    t.psi0 = rbind(sum(w2*(-A2*(A2 - predA2S) + (A2 - predA2S)*XXXXA2)),
                   sum(w2*(-A2*(A2*X21 - X21*predA2S) + (A2*X21 - X21*predA2S)*XXXXA2)));

    t.psi1 = rbind(sum(w2*(-A2*X21*(A2 - predA2S) + (A2 - predA2S)*XXXXA2X21)),
               sum(w2*(-A2*X21*(A2*X21 - X21*predA2S) + (A2*X21 - X21*predA2S)*XXXXA2X21)));

    psi2.star.hat = solve(cbind(t.psi0, t.psi1), -t.0);

    # stage 1
    psi2.hat = psi2.hat.gest.integrate;
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    H1 = cbind(1, X11, X12);
    XXXX = H1%*%solve(t(H1)%*%H1)%*%t(H1);
    XXXXY = XXXX%*%Yopt2.hat;
    XXXXA1 = XXXX%*%A1;
    XXXXA1X11 = XXXX%*%(A1*X11);

    t.0 = rbind(sum(w1*((A1 - predA1S)*Yopt2.hat - (A1 - predA1S)*XXXXY)),
                sum(w1*((A1*X11 - X11*predA1S)*Yopt2.hat - (A1*X11 - X11*predA1S)*XXXXY)));

    t.psi0 = rbind(sum(w1*(-A1*(A1 - predA1S) + (A1 - predA1S)*XXXXA1)),
                   sum(w1*(-A1*(A1*X11 - X11*predA1S) + (A1*X11 - X11*predA1S)*XXXXA1)));

    t.psi1 = rbind(sum(w1*(-A1*X11*(A1 - predA1S) + (A1 - predA1S)*XXXXA1X11)),
                   sum(w1*(-A1*X11*(A1*X11 - X11*predA1S) + (A1*X11 - X11*predA1S)*XXXXA1X11)));

    psi1.star.hat = solve(cbind(t.psi0, t.psi1), -t.0);
    coefs.gest.IPTW[i,] = c(psi1.star.hat, psi2.star.hat);


    ## ICE dWOLS
    Y1_Y0.dWOLS = cbind(1, X21, X22)%*%psi2.hat.dWOLS.integrate;
    mod.dWOLS2.ICE = glm(Y1_Y0.dWOLS~X21);
    Y1_Y0.dWOLS = cbind(1, X11, X12)%*%psi1.hat.dWOLS.integrate;
    mod.dWOLS1.ICE = glm(Y1_Y0.dWOLS~X11);
    coefs.dWOLS.ICE[i,] = c(coef(mod.dWOLS1.ICE), coef(mod.dWOLS2.ICE));


    ## ICE g-est
    Y1_Y0.gest = cbind(1, X21, X22)%*%psi2.hat.gest.integrate;
    mod.gest2.ICE = glm(Y1_Y0.gest~X21);
    Y1_Y0.gest = cbind(1, X11, X12)%*%psi1.hat.gest.integrate;
    mod.gest1.ICE = glm(Y1_Y0.gest~X11);
    coefs.gest.ICE[i,] = c(coef(mod.gest1.ICE), coef(mod.gest2.ICE));


    print(data.frame(n, i, time = Sys.time()));
  }
  
  file.name = paste("results_2TP_scenario1n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
  save(coefs.dWOLS.naive, coefs.gest.naive,
       coefs.dWOLS.integrate, coefs.gest.integrate,
       coefs.dWOLS.IPTW, coefs.gest.IPTW,
       coefs.dWOLS.ICE, coefs.gest.ICE, file = file.name);
} 




### Scernario 2 - Incorrectly specified exposure model
###               Correctly specified outcome model


# In this simulation, both X1 and X2 are confounders,
# but only X1 is used as a tailoring variable 

setwd("C:\\Users\\Denis Talbot\\Dropbox\\Travail\\Recherche\\MSNMM\\Simulations\\Results");
nrep = 1000;
expit = plogis;

# Psi parameters
psi1 = c(-0.5, 1, -1);
psi2 = c(-0.5, 1, -1);



# Monte Carlo estimation of E[X22|X21]
# set.seed(37193179);
# n = 20000000;
# X11 = rbinom(n, 1, 0.5);
# X12 = rbinom(n, 1, 0.5);
# A1 = rbinom(n, 1, expit(-1 + X11 + X12 + X11*X12));
# X21 = rbinom(n, 1, expit(-1 + X11 + A1));
# X22 = rbinom(n, 1, expit(-1 + X12 + A1));
EX22.0 = 0.4648520; #E[X22|X21 = 0]
EX22.1 = 0.5507178; #E[X22|X21 = 1]
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
psi1.true = c(-1.014, 1.028);



set.seed(37193179);
#for(n in c(100,1000,10000)){
for(n in c(100,1000,10000)){

  ## Initialize objects
  coefs.dWOLS.naive = coefs.gest.naive =
  coefs.dWOLS.integrate = coefs.gest.integrate = 
  coefs.dWOLS.IPTW = coefs.gest.IPTW =
  coefs.dWOLS.ICE = coefs.gest.ICE = matrix(, nrow = nrep, ncol = 4);

  for(i in 1:nrep){

    # Stage 1 data generation
    X11 = rbinom(n, 1, 0.5);
    X12 = rbinom(n, 1, 0.5);
    X1 = cbind(1, X11, X12);
    A1 = rbinom(n, 1, expit(-1 + X11 + X12 + X11*X12));
    A1opt = 1*(X1%*% psi1 > 0);
    mu1 = (A1opt - A1)*X1%*%psi1;

    # Stage 2 data generation
    X21 = rbinom(n, 1, expit(-1 + X11 + A1));
    X22 = rbinom(n, 1, expit(-1 + X12 + A1));
    X2 = cbind(1, X21, X22);
    A2 = rbinom(n, 1, expit(-1 + X21 + X22 + X21*X22));
    A2opt = 1*(X2%*% psi2 > 0);
    mu2 = (A2opt - A2)*X2%*%psi2;

    # Outcome generation
    Y = X11 + X12 - mu1 - mu2 + rnorm(n);

  
    ## Propensity score models
    modA1 = glm(A1 ~ X11 + X12, family = "binomial");
    predA1 = predict(modA1, type = "res");

    modA2 = glm(A2 ~ X21 + X22 + X11 + X12 + A1, family = "binomial");
    predA2 = predict(modA2, type = "res");

    modA1S = glm(A1 ~ X11, family = "binomial");
    predA1S = predict(modA1S, type = "res");

    modA2S = glm(A2 ~ X21, family = "binomial");
    predA2S = predict(modA2S, type = "res");


    ## Weights
    w2 = (A2*(predA2S/predA2) + (1 - A2)*(1 - predA2S)/(1 - predA2));
    w1 = A1*(predA1S/predA1) + (1 - A1)*(1 - predA1S)/(1 - predA1);
    w1.dWOLS = abs(A1 - predA1S);
    w1.dWOLS0 = abs(A1 - predA1);
    w2.dWOLS = abs(A2 - predA2S);
    w2.dWOLS0 = abs(A2 - predA2);


    ## Naive dWOLS
    # Stage 2
    m2 = summary(glm(Y~ X11 + X12 + X21 + X22 + A1*(1 + X11 + X12 + X21 + X21) 
              + A2*(1 + X21), weights = w2.dWOLS0, family = "gaussian"))$coef;
    psi2.hat = m2[c("A2", "X21:A2"),1];

    # Stage 1
    H2 = cbind(1, X21);
    A2opt.hat = 1*(H2%*%psi2.hat > 0);
    mu2.hat = (A2opt.hat - A2)*H2%*%psi2.hat;
    Yopt2.hat = Y + mu2.hat;
    m1 = summary(glm(Yopt2.hat ~ X11 + X12 
                     + A1*(1 + X11), weights = w1.dWOLS0))$coef;
    coefs.dWOLS.naive[i,] = c(m1[c("A1", "X11:A1"),1], psi2.hat);


    ## Naive g-est

    # Stage 2

    # X(X'X)^{-1}X'
    H2 = cbind(1, X11, X12, X21, X22, A1, A1*X11, A1*X12, A1*X21, A1*X22);
    XXXX = H2%*%solve(t(H2)%*%H2)%*%t(H2);
    XXXXY = XXXX%*%Y;
    XXXXA2 = XXXX%*%A2;
    XXXXA2X21 = XXXX%*%(A2*X21);

    t.0 = rbind(sum(((A2 - predA2)*Y - (A2 - predA2)*XXXXY)),
                sum(((A2*X21 - X21*predA2)*Y - (A2*X21 - X21*predA2)*XXXXY)));

    t.psi0 = rbind(sum((-A2*(A2 - predA2) + (A2 - predA2)*XXXXA2)),
                   sum((-A2*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2)));

    t.psi1 = rbind(sum((-A2*X21*(A2 - predA2) + (A2 - predA2)*XXXXA2X21)),
               sum((-A2*X21*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2X21)));

    psi2.hat = solve(cbind(t.psi0, t.psi1), -t.0);

    # stage 1
    H2 = cbind(1, X21);
    A2opt.hat = 1*(H2%*%psi2.hat > 0);
    mu2.hat = (A2opt.hat - A2)*H2%*%psi2.hat;
    Yopt2.hat = Y + mu2.hat;

    H1 = cbind(1, X11, X12);
    XXXX = H1%*%solve(t(H1)%*%H1)%*%t(H1);
    XXXXY = XXXX%*%Yopt2.hat;
    XXXXA1 = XXXX%*%A1;
    XXXXA1X11 = XXXX%*%(A1*X11);

    t.0 = rbind(sum(((A1 - predA1)*Yopt2.hat - (A1 - predA1)*XXXXY)),
                sum(((A1*X11 - X11*predA1)*Yopt2.hat - (A1*X11 - X11*predA1)*XXXXY)));

    t.psi0 = rbind(sum((-A1*(A1 - predA1) + (A1 - predA1)*XXXXA1)),
                   sum((-A1*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1)));

    t.psi1 = rbind(sum((-A1*X11*(A1 - predA1) + (A1 - predA1)*XXXXA1X11)),
                   sum((-A1*X11*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1X11)));

    psi1.hat = solve(cbind(t.psi0, t.psi1), -t.0);
    coefs.gest.naive[i,] = c(psi1.hat, psi2.hat);


    ## Integrate dWOLS
    # Stage 2
    m2 = summary(glm(Y~ X11 + X12 + X21 + X22 + A1*(1 + X11 + X12 + X21 + X21) 
              + A2*(1 + X21 + X22), weights = w2.dWOLS0, family = "gaussian"))$coef;
    psi2.hat = m2[c("A2", "X21:A2", "X22:A2"),1];
    psi2.hat.dWOLS.integrate = psi2.hat;
    EX12 = prop.table(table(X12, X11), 2)[2,];
    EX22 = prop.table(table(X22, X21), 2)[2,];
    psi2.star.hat = c(psi2.hat[1] + psi2.hat[3]*EX22[1],
                      psi2.hat[2] + psi2.hat[3]*(EX22[2] - EX22[1]));

    # Stage 1
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    m1 = summary(glm(Yopt2.hat ~ X11 + X12 
                     + A1*(1 + X11 + X12), weights = w1.dWOLS0))$coef;
    psi1.hat = m1[c("A1", "X11:A1", "X12:A1"),1];
    psi1.hat.dWOLS.integrate = psi1.hat;
    psi1.star.hat = c(psi1.hat[1] + psi1.hat[3]*EX12[1],
                      psi1.hat[2] + psi1.hat[3]*(EX12[2] - EX12[1]));

    coefs.dWOLS.integrate[i,] = c(psi1.star.hat, psi2.star.hat);


    ## Integrate g-est
    # Stage 2
    # X(X'X)^{-1}X'
    H2 = cbind(1, X11, X12, X21, X22, A1, A1*X11, A1*X12, A1*X21, A1*X22);
    XXXX = H2%*%solve(t(H2)%*%H2)%*%t(H2);
    XXXXY = XXXX%*%Y;
    XXXXA2 = XXXX%*%A2;
    XXXXA2X21 = XXXX%*%(A2*X21);
    XXXXA2X22 = XXXX%*%(A2*X22);

    t.0 = rbind(sum(((A2 - predA2)*Y - (A2 - predA2)*XXXXY)),
                sum(((A2*X21 - X21*predA2)*Y - (A2*X21 - X21*predA2)*XXXXY)),
                sum(((A2*X22 - X22*predA2)*Y - (A2*X22 - X22*predA2)*XXXXY)));

    t.psi0 = rbind(sum((-A2*(A2 - predA2) + (A2 - predA2)*XXXXA2)),
                   sum((-A2*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2)),
                   sum((-A2*(A2*X22 - X22*predA2) + (A2*X22 - X22*predA2)*XXXXA2)));

    t.psi1 = rbind(sum((-A2*X21*(A2 - predA2) + (A2 - predA2)*XXXXA2X21)),
               sum((-A2*X21*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2X21)),
               sum((-A2*X21*(A2*X22 - X22*predA2) + (A2*X22 - X22*predA2)*XXXXA2X21)));

    t.psi2 = rbind(sum((-A2*X22*(A2 - predA2) + (A2 - predA2)*XXXXA2X22)),
               sum((-A2*X22*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2X22)),
               sum((-A2*X22*(A2*X22 - X22*predA2) + (A2*X22 - X22*predA2)*XXXXA2X22)));

    psi2.hat = solve(cbind(t.psi0, t.psi1, t.psi2), -t.0);
    psi2.hat.gest.integrate = psi2.hat;
    psi2.star.hat = c(psi2.hat[1] + psi2.hat[3]*EX22[1],
                      psi2.hat[2] + psi2.hat[3]*(EX22[2] - EX22[1]));

    # stage 1
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    H1 = cbind(1, X11, X12);
    XXXX = H1%*%solve(t(H1)%*%H1)%*%t(H1);
    XXXXY = XXXX%*%Yopt2.hat;
    XXXXA1 = XXXX%*%A1;
    XXXXA1X11 = XXXX%*%(A1*X11);
    XXXXA1X12 = XXXX%*%(A1*X12);

    t.0 = rbind(sum(((A1 - predA1)*Yopt2.hat - (A1 - predA1)*XXXXY)),
                sum(((A1*X11 - X11*predA1)*Yopt2.hat - (A1*X11 - X11*predA1)*XXXXY)),
                sum(((A1*X12 - X12*predA1)*Yopt2.hat - (A1*X12 - X12*predA1)*XXXXY)));

    t.psi0 = rbind(sum((-A1*(A1 - predA1) + (A1 - predA1)*XXXXA1)),
                   sum((-A1*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1)),
                   sum((-A1*(A1*X12 - X12*predA1) + (A1*X12 - X12*predA1)*XXXXA1)));

    t.psi1 = rbind(sum((-A1*X11*(A1 - predA1) + (A1 - predA1)*XXXXA1X11)),
                   sum((-A1*X11*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1X11)),
                   sum((-A1*X11*(A1*X12 - X12*predA1) + (A1*X12 - X12*predA1)*XXXXA1X11)));

    t.psi2 = rbind(sum((-A1*X12*(A1 - predA1) + (A1 - predA1)*XXXXA1X12)),
                   sum((-A1*X12*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1X12)),
                   sum((-A1*X12*(A1*X12 - X12*predA1) + (A1*X12 - X12*predA1)*XXXXA1X12)));

    psi1.hat = solve(cbind(t.psi0, t.psi1, t.psi2), -t.0);
    psi1.hat.gest.integrate = psi1.hat;
    psi1.star.hat = c(psi1.hat[1] + psi1.hat[3]*EX12[1],
                      psi1.hat[2] + psi1.hat[3]*(EX12[2] - EX12[1]));
    coefs.gest.integrate[i,] =  c(psi1.star.hat, psi2.star.hat);


    ## IPTW + dWOLS
    # Stage 2
    m2 = summary(glm(Y~ X11 + X12 + X21 + X22 + A1*(1 + X11 + X12 + X21 + X21) 
              + A2*(1 + X21), weights = w2.dWOLS*w2, family = "gaussian"))$coef;
    psi2.star.hat = m2[c("A2", "X21:A2"),1];

    # Stage 1
    psi2.hat = psi2.hat.dWOLS.integrate;
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    m1 = summary(glm(Yopt2.hat ~ X11 + X12 
                     + A1*(1 + X11), weights = w1*w1.dWOLS))$coef;
    coefs.dWOLS.IPTW[i,] = c(m1[c("A1", "X11:A1"),1], psi2.star.hat);


    ## IPTW + gest
    # Stage 2

    # X(X'X)^{-1}X'
    H2 = cbind(1, X11, X12, X21, X22, A1, A1*X11, A1*X12, A1*X21, A1*X22);
    XXXX = H2%*%solve(t(H2)%*%H2)%*%t(H2);
    XXXXY = XXXX%*%Y;
    XXXXA2 = XXXX%*%A2;
    XXXXA2X21 = XXXX%*%(A2*X21);

    t.0 = rbind(sum(w2*((A2 - predA2S)*Y - (A2 - predA2S)*XXXXY)),
                sum(w2*((A2*X21 - X21*predA2S)*Y - (A2*X21 - X21*predA2S)*XXXXY)));

    t.psi0 = rbind(sum(w2*(-A2*(A2 - predA2S) + (A2 - predA2S)*XXXXA2)),
                   sum(w2*(-A2*(A2*X21 - X21*predA2S) + (A2*X21 - X21*predA2S)*XXXXA2)));

    t.psi1 = rbind(sum(w2*(-A2*X21*(A2 - predA2S) + (A2 - predA2S)*XXXXA2X21)),
               sum(w2*(-A2*X21*(A2*X21 - X21*predA2S) + (A2*X21 - X21*predA2S)*XXXXA2X21)));

    psi2.star.hat = solve(cbind(t.psi0, t.psi1), -t.0);

    # stage 1
    psi2.hat = psi2.hat.gest.integrate;
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    H1 = cbind(1, X11, X12);
    XXXX = H1%*%solve(t(H1)%*%H1)%*%t(H1);
    XXXXY = XXXX%*%Yopt2.hat;
    XXXXA1 = XXXX%*%A1;
    XXXXA1X11 = XXXX%*%(A1*X11);

    t.0 = rbind(sum(w1*((A1 - predA1S)*Yopt2.hat - (A1 - predA1S)*XXXXY)),
                sum(w1*((A1*X11 - X11*predA1S)*Yopt2.hat - (A1*X11 - X11*predA1S)*XXXXY)));

    t.psi0 = rbind(sum(w1*(-A1*(A1 - predA1S) + (A1 - predA1S)*XXXXA1)),
                   sum(w1*(-A1*(A1*X11 - X11*predA1S) + (A1*X11 - X11*predA1S)*XXXXA1)));

    t.psi1 = rbind(sum(w1*(-A1*X11*(A1 - predA1S) + (A1 - predA1S)*XXXXA1X11)),
                   sum(w1*(-A1*X11*(A1*X11 - X11*predA1S) + (A1*X11 - X11*predA1S)*XXXXA1X11)));

    psi1.star.hat = solve(cbind(t.psi0, t.psi1), -t.0);
    coefs.gest.IPTW[i,] = c(psi1.star.hat, psi2.star.hat);


    ## ICE dWOLS
    Y1_Y0.dWOLS = cbind(1, X21, X22)%*%psi2.hat.dWOLS.integrate;
    mod.dWOLS2.ICE = glm(Y1_Y0.dWOLS~X21);
    Y1_Y0.dWOLS = cbind(1, X11, X12)%*%psi1.hat.dWOLS.integrate;
    mod.dWOLS1.ICE = glm(Y1_Y0.dWOLS~X11);
    coefs.dWOLS.ICE[i,] = c(coef(mod.dWOLS1.ICE), coef(mod.dWOLS2.ICE));


    ## ICE g-est
    Y1_Y0.gest = cbind(1, X21, X22)%*%psi2.hat.gest.integrate;
    mod.gest2.ICE = glm(Y1_Y0.gest~X21);
    Y1_Y0.gest = cbind(1, X11, X12)%*%psi1.hat.gest.integrate;
    mod.gest1.ICE = glm(Y1_Y0.gest~X11);
    coefs.gest.ICE[i,] = c(coef(mod.gest1.ICE), coef(mod.gest2.ICE));


    print(data.frame(n, i, time = Sys.time()));
  }
  
  file.name = paste("results_2TP_scenario2n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
  save(coefs.dWOLS.naive, coefs.gest.naive,
       coefs.dWOLS.integrate, coefs.gest.integrate,
       coefs.dWOLS.IPTW, coefs.gest.IPTW,
       coefs.dWOLS.ICE, coefs.gest.ICE, file = file.name);
} 




### Scernario 3 - Correctly specified exposure model
###               Correctly specified outcome model


# In this simulation, both X1 and X2 are confounders,
# but only X1 is used as a tailoring variable 

setwd("C:\\Users\\Denis Talbot\\Dropbox\\Travail\\Recherche\\MSNMM\\Simulations\\Results");
nrep = 1000;
expit = plogis;

# Psi parameters
psi1 = c(-0.5, 1, -1);
psi2 = c(-0.5, 1, -1);

### Scernario 1 - Correctly specified exposure model
###               Incorrectly specified outcome model
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
# Y.0 = X11 + X12 + X11*X12 - mu1.0 - mu2.0;
# Y.1 = X11 + X12 + X11*X12 - mu1.1 - mu2.1;
# glm(Y.1 - Y.0 ~ X11);
psi1.true = c(-1.014, 1.028);



set.seed(37193179);
for(n in c(100,1000,10000)){

  ## Initialize objects
  coefs.dWOLS.naive = coefs.gest.naive =
  coefs.dWOLS.integrate = coefs.gest.integrate = 
  coefs.dWOLS.IPTW = coefs.gest.IPTW =
  coefs.dWOLS.ICE = coefs.gest.ICE = matrix(, nrow = nrep, ncol = 4);

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
    Y = X11 + X12 + X11*X12 - mu1 - mu2 + rnorm(n);

  
    ## Propensity score models
    modA1 = glm(A1 ~ X11 + X12, family = "binomial");
    predA1 = predict(modA1, type = "res");

    modA2 = glm(A2 ~ X21 + X22 + X11 + X12 + A1, family = "binomial");
    predA2 = predict(modA2, type = "res");

    modA1S = glm(A1 ~ X11, family = "binomial");
    predA1S = predict(modA1S, type = "res");

    modA2S = glm(A2 ~ X21, family = "binomial");
    predA2S = predict(modA2S, type = "res");


    ## Weights
    w2 = (A2*(predA2S/predA2) + (1 - A2)*(1 - predA2S)/(1 - predA2));
    w1 = A1*(predA1S/predA1) + (1 - A1)*(1 - predA1S)/(1 - predA1);
    w1.dWOLS = abs(A1 - predA1S);
    w1.dWOLS0 = abs(A1 - predA1);
    w2.dWOLS = abs(A2 - predA2S);
    w2.dWOLS0 = abs(A2 - predA2);


    ## Naive dWOLS
    # Stage 2
    m2 = summary(glm(Y~ X11 + X12 + X21 + X22 + A1*(1 + X11 + X12 + X21 + X21) 
              + A2*(1 + X21), weights = w2.dWOLS0, family = "gaussian"))$coef;
    psi2.hat = m2[c("A2", "X21:A2"),1];

    # Stage 1
    H2 = cbind(1, X21);
    A2opt.hat = 1*(H2%*%psi2.hat > 0);
    mu2.hat = (A2opt.hat - A2)*H2%*%psi2.hat;
    Yopt2.hat = Y + mu2.hat;
    m1 = summary(glm(Yopt2.hat ~ X11 + X12 
                     + A1*(1 + X11), weights = w1.dWOLS0))$coef;
    coefs.dWOLS.naive[i,] = c(m1[c("A1", "X11:A1"),1], psi2.hat);


    ## Naive g-est

    # Stage 2

    # X(X'X)^{-1}X'
    H2 = cbind(1, X11, X12, X21, X22, A1, A1*X11, A1*X12, A1*X21, A1*X22);
    XXXX = H2%*%solve(t(H2)%*%H2)%*%t(H2);
    XXXXY = XXXX%*%Y;
    XXXXA2 = XXXX%*%A2;
    XXXXA2X21 = XXXX%*%(A2*X21);

    t.0 = rbind(sum(((A2 - predA2)*Y - (A2 - predA2)*XXXXY)),
                sum(((A2*X21 - X21*predA2)*Y - (A2*X21 - X21*predA2)*XXXXY)));

    t.psi0 = rbind(sum((-A2*(A2 - predA2) + (A2 - predA2)*XXXXA2)),
                   sum((-A2*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2)));

    t.psi1 = rbind(sum((-A2*X21*(A2 - predA2) + (A2 - predA2)*XXXXA2X21)),
               sum((-A2*X21*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2X21)));

    psi2.hat = solve(cbind(t.psi0, t.psi1), -t.0);

    # stage 1
    H2 = cbind(1, X21);
    A2opt.hat = 1*(H2%*%psi2.hat > 0);
    mu2.hat = (A2opt.hat - A2)*H2%*%psi2.hat;
    Yopt2.hat = Y + mu2.hat;

    H1 = cbind(1, X11, X12);
    XXXX = H1%*%solve(t(H1)%*%H1)%*%t(H1);
    XXXXY = XXXX%*%Yopt2.hat;
    XXXXA1 = XXXX%*%A1;
    XXXXA1X11 = XXXX%*%(A1*X11);

    t.0 = rbind(sum(((A1 - predA1)*Yopt2.hat - (A1 - predA1)*XXXXY)),
                sum(((A1*X11 - X11*predA1)*Yopt2.hat - (A1*X11 - X11*predA1)*XXXXY)));

    t.psi0 = rbind(sum((-A1*(A1 - predA1) + (A1 - predA1)*XXXXA1)),
                   sum((-A1*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1)));

    t.psi1 = rbind(sum((-A1*X11*(A1 - predA1) + (A1 - predA1)*XXXXA1X11)),
                   sum((-A1*X11*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1X11)));

    psi1.hat = solve(cbind(t.psi0, t.psi1), -t.0);
    coefs.gest.naive[i,] = c(psi1.hat, psi2.hat);


    ## Integrate dWOLS
    # Stage 2
    m2 = summary(glm(Y~ X11 + X12 + X21 + X22 + A1*(1 + X11 + X12 + X21 + X21) 
              + A2*(1 + X21 + X22), weights = w2.dWOLS0, family = "gaussian"))$coef;
    psi2.hat = m2[c("A2", "X21:A2", "X22:A2"),1];
    psi2.hat.dWOLS.integrate = psi2.hat;
    EX12 = prop.table(table(X12, X11), 2)[2,];
    EX22 = prop.table(table(X22, X21), 2)[2,];
    psi2.star.hat = c(psi2.hat[1] + psi2.hat[3]*EX22[1],
                      psi2.hat[2] + psi2.hat[3]*(EX22[2] - EX22[1]));

    # Stage 1
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    m1 = summary(glm(Yopt2.hat ~ X11 + X12 
                     + A1*(1 + X11 + X12), weights = w1.dWOLS0))$coef;
    psi1.hat = m1[c("A1", "X11:A1", "X12:A1"),1];
    psi1.hat.dWOLS.integrate = psi1.hat;
    psi1.star.hat = c(psi1.hat[1] + psi1.hat[3]*EX12[1],
                      psi1.hat[2] + psi1.hat[3]*(EX12[2] - EX12[1]));

    coefs.dWOLS.integrate[i,] = c(psi1.star.hat, psi2.star.hat);


    ## Integrate g-est
    # Stage 2
    # X(X'X)^{-1}X'
    H2 = cbind(1, X11, X12, X21, X22, A1, A1*X11, A1*X12, A1*X21, A1*X22);
    XXXX = H2%*%solve(t(H2)%*%H2)%*%t(H2);
    XXXXY = XXXX%*%Y;
    XXXXA2 = XXXX%*%A2;
    XXXXA2X21 = XXXX%*%(A2*X21);
    XXXXA2X22 = XXXX%*%(A2*X22);

    t.0 = rbind(sum(((A2 - predA2)*Y - (A2 - predA2)*XXXXY)),
                sum(((A2*X21 - X21*predA2)*Y - (A2*X21 - X21*predA2)*XXXXY)),
                sum(((A2*X22 - X22*predA2)*Y - (A2*X22 - X22*predA2)*XXXXY)));

    t.psi0 = rbind(sum((-A2*(A2 - predA2) + (A2 - predA2)*XXXXA2)),
                   sum((-A2*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2)),
                   sum((-A2*(A2*X22 - X22*predA2) + (A2*X22 - X22*predA2)*XXXXA2)));

    t.psi1 = rbind(sum((-A2*X21*(A2 - predA2) + (A2 - predA2)*XXXXA2X21)),
               sum((-A2*X21*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2X21)),
               sum((-A2*X21*(A2*X22 - X22*predA2) + (A2*X22 - X22*predA2)*XXXXA2X21)));

    t.psi2 = rbind(sum((-A2*X22*(A2 - predA2) + (A2 - predA2)*XXXXA2X22)),
               sum((-A2*X22*(A2*X21 - X21*predA2) + (A2*X21 - X21*predA2)*XXXXA2X22)),
               sum((-A2*X22*(A2*X22 - X22*predA2) + (A2*X22 - X22*predA2)*XXXXA2X22)));

    psi2.hat = solve(cbind(t.psi0, t.psi1, t.psi2), -t.0);
    psi2.hat.gest.integrate = psi2.hat;
    psi2.star.hat = c(psi2.hat[1] + psi2.hat[3]*EX22[1],
                      psi2.hat[2] + psi2.hat[3]*(EX22[2] - EX22[1]));

    # stage 1
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    H1 = cbind(1, X11, X12);
    XXXX = H1%*%solve(t(H1)%*%H1)%*%t(H1);
    XXXXY = XXXX%*%Yopt2.hat;
    XXXXA1 = XXXX%*%A1;
    XXXXA1X11 = XXXX%*%(A1*X11);
    XXXXA1X12 = XXXX%*%(A1*X12);

    t.0 = rbind(sum(((A1 - predA1)*Yopt2.hat - (A1 - predA1)*XXXXY)),
                sum(((A1*X11 - X11*predA1)*Yopt2.hat - (A1*X11 - X11*predA1)*XXXXY)),
                sum(((A1*X12 - X12*predA1)*Yopt2.hat - (A1*X12 - X12*predA1)*XXXXY)));

    t.psi0 = rbind(sum((-A1*(A1 - predA1) + (A1 - predA1)*XXXXA1)),
                   sum((-A1*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1)),
                   sum((-A1*(A1*X12 - X12*predA1) + (A1*X12 - X12*predA1)*XXXXA1)));

    t.psi1 = rbind(sum((-A1*X11*(A1 - predA1) + (A1 - predA1)*XXXXA1X11)),
                   sum((-A1*X11*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1X11)),
                   sum((-A1*X11*(A1*X12 - X12*predA1) + (A1*X12 - X12*predA1)*XXXXA1X11)));

    t.psi2 = rbind(sum((-A1*X12*(A1 - predA1) + (A1 - predA1)*XXXXA1X12)),
                   sum((-A1*X12*(A1*X11 - X11*predA1) + (A1*X11 - X11*predA1)*XXXXA1X12)),
                   sum((-A1*X12*(A1*X12 - X12*predA1) + (A1*X12 - X12*predA1)*XXXXA1X12)));

    psi1.hat = solve(cbind(t.psi0, t.psi1, t.psi2), -t.0);
    psi1.hat.gest.integrate = psi1.hat;
    psi1.star.hat = c(psi1.hat[1] + psi1.hat[3]*EX12[1],
                      psi1.hat[2] + psi1.hat[3]*(EX12[2] - EX12[1]));
    coefs.gest.integrate[i,] =  c(psi1.star.hat, psi2.star.hat);


    ## IPTW + dWOLS
    # Stage 2
    m2 = summary(glm(Y~ X11 + X12 + X21 + X22 + A1*(1 + X11 + X12 + X21 + X21) 
              + A2*(1 + X21), weights = w2.dWOLS*w2, family = "gaussian"))$coef;
    psi2.star.hat = m2[c("A2", "X21:A2"),1];

    # Stage 1
    psi2.hat = psi2.hat.dWOLS.integrate;
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    m1 = summary(glm(Yopt2.hat ~ X11 + X12 
                     + A1*(1 + X11), weights = w1*w1.dWOLS))$coef;
    coefs.dWOLS.IPTW[i,] = c(m1[c("A1", "X11:A1"),1], psi2.star.hat);


    ## IPTW + gest
    # Stage 2

    # X(X'X)^{-1}X'
    H2 = cbind(1, X11, X12, X21, X22, A1, A1*X11, A1*X12, A1*X21, A1*X22);
    XXXX = H2%*%solve(t(H2)%*%H2)%*%t(H2);
    XXXXY = XXXX%*%Y;
    XXXXA2 = XXXX%*%A2;
    XXXXA2X21 = XXXX%*%(A2*X21);

    t.0 = rbind(sum(w2*((A2 - predA2S)*Y - (A2 - predA2S)*XXXXY)),
                sum(w2*((A2*X21 - X21*predA2S)*Y - (A2*X21 - X21*predA2S)*XXXXY)));

    t.psi0 = rbind(sum(w2*(-A2*(A2 - predA2S) + (A2 - predA2S)*XXXXA2)),
                   sum(w2*(-A2*(A2*X21 - X21*predA2S) + (A2*X21 - X21*predA2S)*XXXXA2)));

    t.psi1 = rbind(sum(w2*(-A2*X21*(A2 - predA2S) + (A2 - predA2S)*XXXXA2X21)),
               sum(w2*(-A2*X21*(A2*X21 - X21*predA2S) + (A2*X21 - X21*predA2S)*XXXXA2X21)));

    psi2.star.hat = solve(cbind(t.psi0, t.psi1), -t.0);

    # stage 1
    psi2.hat = psi2.hat.gest.integrate;
    H2 = cbind(1, X21, X22);
    H2.star = cbind(1, X21);
    A2opt.star.hat = 1*(H2.star%*%psi2.star.hat > 0);
    blip.A2.hat = A2*(H2%*%psi2.hat);
    blip.A2opt.star.hat = A2opt.star.hat*(H2.star%*%psi2.star.hat);
    Yopt2.hat = Y - blip.A2.hat + blip.A2opt.star.hat;

    H1 = cbind(1, X11, X12);
    XXXX = H1%*%solve(t(H1)%*%H1)%*%t(H1);
    XXXXY = XXXX%*%Yopt2.hat;
    XXXXA1 = XXXX%*%A1;
    XXXXA1X11 = XXXX%*%(A1*X11);

    t.0 = rbind(sum(w1*((A1 - predA1S)*Yopt2.hat - (A1 - predA1S)*XXXXY)),
                sum(w1*((A1*X11 - X11*predA1S)*Yopt2.hat - (A1*X11 - X11*predA1S)*XXXXY)));

    t.psi0 = rbind(sum(w1*(-A1*(A1 - predA1S) + (A1 - predA1S)*XXXXA1)),
                   sum(w1*(-A1*(A1*X11 - X11*predA1S) + (A1*X11 - X11*predA1S)*XXXXA1)));

    t.psi1 = rbind(sum(w1*(-A1*X11*(A1 - predA1S) + (A1 - predA1S)*XXXXA1X11)),
                   sum(w1*(-A1*X11*(A1*X11 - X11*predA1S) + (A1*X11 - X11*predA1S)*XXXXA1X11)));

    psi1.star.hat = solve(cbind(t.psi0, t.psi1), -t.0);
    coefs.gest.IPTW[i,] = c(psi1.star.hat, psi2.star.hat);


    ## ICE dWOLS
    Y1_Y0.dWOLS = cbind(1, X21, X22)%*%psi2.hat.dWOLS.integrate;
    mod.dWOLS2.ICE = glm(Y1_Y0.dWOLS~X21);
    Y1_Y0.dWOLS = cbind(1, X11, X12)%*%psi1.hat.dWOLS.integrate;
    mod.dWOLS1.ICE = glm(Y1_Y0.dWOLS~X11);
    coefs.dWOLS.ICE[i,] = c(coef(mod.dWOLS1.ICE), coef(mod.dWOLS2.ICE));


    ## ICE g-est
    Y1_Y0.gest = cbind(1, X21, X22)%*%psi2.hat.gest.integrate;
    mod.gest2.ICE = glm(Y1_Y0.gest~X21);
    Y1_Y0.gest = cbind(1, X11, X12)%*%psi1.hat.gest.integrate;
    mod.gest1.ICE = glm(Y1_Y0.gest~X11);
    coefs.gest.ICE[i,] = c(coef(mod.gest1.ICE), coef(mod.gest2.ICE));


    print(data.frame(n, i, time = Sys.time()));
  }
  
  file.name = paste("results_2TP_scenario3n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
  save(coefs.dWOLS.naive, coefs.gest.naive,
       coefs.dWOLS.integrate, coefs.gest.integrate,
       coefs.dWOLS.IPTW, coefs.gest.IPTW,
       coefs.dWOLS.ICE, coefs.gest.ICE, file = file.name);
} 





