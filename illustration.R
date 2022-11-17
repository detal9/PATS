#### OR+BMI as tailoring, PR as neglected interaction

#### Loading required libraries
require(tableone);
require(stargazer);
require(survey);
require(boot);
require(rms);
require(survival);


#### Set working directory
setwd("P:\\SPPOS\\Equipe_Talbot_D\\Projet_Diorio\\version2");


#### Read dataset
ds = read.csv("cms_er_bmi.csv", sep = ";", na.strings = c("."));
head(ds);
summary(ds);


#### Change dates to date format
ds$d_diag = as.Date(ds$d_diag, "%d-%m-%Y");
ds$date_bx_sein = as.Date(ds$date_bx_sein, "%d-%m-%Y");
ds$date_chx_sein = as.Date(ds$date_chx_sein, "%d-%m-%Y");
ds$date_recep = as.Date(ds$date_recep, "%d-%m-%Y");
ds$dfin_2011 = as.Date(ds$dfin_2011, "%d-%m-%Y");


### Restrict attention to certain participants
ds2 = ds[ds$d_diag >= as.Date("01-01-1987", "%d-%m-%Y") & ds$d_diag < as.Date("01-01-2010", "%d-%m-%Y"),];


#### Data manipulation
## Outcome and censoring
ds2$censor = (ds2$statut_2011 == 1)*1;
ds2$Y = ds2$suivi_2011a;

## Treatment 
ds2$A = ds2$ind_hor;

## Potential confounders
ds2$cbmi = ifelse(ds2$bmi > 25, 1, 0);
ds2$meno = ds2$meno_statut; 
ds2$meno[ds2$meno_statut == 9] = ifelse(ds2$cage[ds2$meno_statut == 9] >= 3, 2, 1);
ds2$staging_rvalue[ds2$staging_rvalue == 9] = NA;
ds2$ctabac[ds2$ctabac %in% c(2,3,4)] = 2;
ds2$ctabac = factor(ds2$ctabac);
ds2$histfm_1degre[ds2$histfm_1degre == 9] = NA;
ds2$prise_hts[ds2$prise_hts %in% c(2,3,4)] = 2;
ds2$prise_hts[ds2$prise_hts == 9] = NA;
ds2$ind_chi[ds2$ind_chi == 9] = 0;
ds2$ind_rad[ds2$ind_rad == 9] = 1;
ds2$ind_herceptin[ds2$ind_herceptin == 9] = NA;
ds2$year = as.numeric(format(ds2$d_diag, '%Y'));
ds2$cage = factor(ds2$cage);
ds2$cbmi = factor(ds2$cbmi);
ds2$meno = factor(ds2$meno);
ds2$grade = factor(ds2$grade);
ds2$staging_rvalue = factor(ds2$staging_rvalue);
ds2$recod_oest = factor(ds2$recod_oest);
ds2$recod_prog = factor(ds2$recod_prog);
ds2$type_chx = factor(ds2$type_chx);
ds2$histfm_1degre = factor(ds2$histfm_1degre);
ds2$prise_hts = factor(ds2$prise_hts);
ds2$ind_chi = factor(ds2$ind_chi);
ds2$ind_rad = factor(ds2$ind_rad);
ds2$ind_herceptin = factor(ds2$ind_herceptin);

ds2$yearcat = ifelse(ds2$d_diag >= as.Date("01-01-1985", "%d-%m-%Y") & ds2$d_diag <= as.Date("31-12-1989", "%d-%m-%Y"), "85-89",
               ifelse(ds2$d_diag >= as.Date("01-01-1990", "%d-%m-%Y") & ds2$d_diag <= as.Date("31-12-1994", "%d-%m-%Y"), "90-94", 
                ifelse(ds2$d_diag >= as.Date("01-01-1995", "%d-%m-%Y") & ds2$d_diag <= as.Date("31-12-1999", "%d-%m-%Y"), "95-99",
                 ifelse(ds2$d_diag >= as.Date("01-01-2000", "%d-%m-%Y") & ds2$d_diag <= as.Date("31-12-2004", "%d-%m-%Y"), "00-04", "05-09"))));


#### Remove obs w/ missing data
cov_names = c("cage", "ctabac", "bmi", "cbmi", "meno", "grade", "staging_rvalue", "recod_oest", "recod_prog",
              "type_chx", "histfm_1degre", "prise_hts", "ind_chi", "ind_rad", "ind_herceptin",
              "cd_diag", "year", "yearcat")
ds3 = ds2[, c("A", "Y", "censor", cov_names)];
ds4 = ds3[complete.cases(ds3[,-2]),];


#### Descriptive statistics
print(CreateTableOne(vars = cov_names, strata = "A", data = ds4), test = FALSE, smd = TRUE);
stargazer(print(CreateTableOne(vars = cov_names, strata = "A", data = ds4), test = FALSE, smd = TRUE));

print(CreateTableOne(vars = cov_names, strata = "censor", data = ds4), test = FALSE, smd = TRUE);
stargazer(print(CreateTableOne(vars = cov_names, strata = "censor", data = ds4), test = FALSE, smd = TRUE));

by(ds4$Y + ds4$year, ds4$censor, summary);
plot(ds4$Y, ds4$year);

# ds4$censor: 0
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1988    1998    2004    2003    2008    2012 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
# ds4$censor: 1
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1994    2011    2011    2011    2012    2012 



#### Inverse probability of treatment weighting
mod.iptw = glm(A ~ cage + cbmi + ctabac + meno + grade + staging_rvalue +
               recod_oest + recod_prog + type_chx + histfm_1degre + prise_hts +
               ind_chi + ind_rad + ind_herceptin + rcs(year),
               data = ds4, family = "binomial");
summary(mod.iptw);
iptw = abs(ds4$A - predict(mod.iptw, type = "res"));
summary(iptw);



#### Naive dWOLS
ds4$oes_bmi = paste0("oest=",ds4$recod_oest,",cbmi=",ds4$cbmi);

naive = survreg(Surv(Y, 1 - censor) ~ cage + oes_bmi + ctabac + meno + grade + staging_rvalue +
                                      recod_prog + type_chx + histfm_1degre + prise_hts +
                                      ind_chi + ind_rad + ind_herceptin + rcs(year) +
                                      A*(1 + oes_bmi),
                weights = iptw, data = ds4, dist = "loglogistic");
summary(naive);

# log-norm     : 3501.74
# Weibull      : 3512.165
# log-logistic : 3493.276
# exp          : 3530.706



#### CE dWOLS
PATS1 = survreg(Surv(Y, 1 - censor) ~ cage + oes_bmi + ctabac + meno + grade + staging_rvalue +
                                      recod_prog + type_chx + histfm_1degre + prise_hts +
                                      ind_chi + ind_rad + ind_herceptin + rcs(year) +
                                      A*(1 + oes_bmi + recod_prog),
                weights = iptw, data = ds4, dist = "loglogistic");

ds4$Q = with(ds4, model.matrix(Y ~ oes_bmi + recod_prog))%*%tail(coef(PATS1), 6);
lm(Q ~ oes_bmi, data = ds4);


#### Confidence intervals w/bootstrap
bootf = function(ds, index){
  ds4 = ds[index,];


  # Inverse probability of treatment weighting
  mod.iptw = glm(A ~ cage + cbmi + ctabac + meno + grade + staging_rvalue +
                 recod_oest + recod_prog + type_chx + histfm_1degre + prise_hts +
                 ind_chi + ind_rad + ind_herceptin + rcs(year),
                 data = ds4, family = "binomial");
  iptw = abs(ds4$A - predict(mod.iptw, type = "res"));


  # Naive dWOLS
  ds4$oes_bmi = paste0("oest=",ds4$recod_oest,",cbmi=",ds4$cbmi);

  naive = survreg(Surv(Y, 1 - censor) ~ cage + oes_bmi + ctabac + meno + grade + staging_rvalue +
                                        recod_prog + type_chx + histfm_1degre + prise_hts +
                                        ind_chi + ind_rad + ind_herceptin + rcs(year) +
                                        A*(1 + oes_bmi),
                  weights = iptw, data = ds4, dist = "loglogistic");
  psi.naive = tail(coef(naive),4);


  # CE dWOLS
  PATS1 = survreg(Surv(Y, 1 - censor) ~ cage + oes_bmi + ctabac + meno + grade + staging_rvalue +
                                        recod_prog + type_chx + histfm_1degre + prise_hts +
                                        ind_chi + ind_rad + ind_herceptin + rcs(year) +
                                        A*(1 + oes_bmi + recod_prog),
                  weights = iptw, data = ds4, dist = "loglogistic");
  ds4$Q = with(ds4, model.matrix(Y ~ oes_bmi + recod_prog))%*%tail(coef(PATS1), 6);
  psi.PATS1 = coef(lm(Q ~ oes_bmi, data = ds4));
  results = c(psi.naive, psi.PATS1);
  return(results);
}
set.seed(7481971);
bootres = boot(ds4, bootf, 1000);

cbind(bootres$t0[1], bootres$t0[1] + bootres$t0[2],  bootres$t0[1] + bootres$t0[3],  bootres$t0[1] + bootres$t0[4],
      bootres$t0[5], bootres$t0[5] + bootres$t0[6],  bootres$t0[5] + bootres$t0[7],  bootres$t0[5] + bootres$t0[8]);
#       [,1]      [,2]        [,3]       [,4]      [,5]      [,6]        [,7]       [,8]
# A 0.154664 0.1981747 -0.07453349 -0.3263892 0.1775836 0.2202928 -0.06584198 -0.3313612



psi = cbind(bootres$t[,1], bootres$t[,1] + bootres$t[,2],  bootres$t[,1] + bootres$t[,3],  bootres$t[,1] + bootres$t[,4],
            bootres$t[,5], bootres$t[,5] + bootres$t[,6],  bootres$t[,5] + bootres$t[,7],  bootres$t[,5] + bootres$t[,8]);

### BMI > 25
# [1]  0.15331432  0.19658338 -0.05832835 -0.32350693
#      0.17615475  0.21820802 -0.04889002 -0.32753882

# > apply(psi, 2, quantile, c(0.025, 0.975));
#              [,1]       [,2]       [,3]       [,4]     
# 2.5%  -0.01316625 0.005094372 -0.4545959 -0.7845331 
# 97.5%  0.30977457 0.389947583  0.3382736  0.1095211 
# 2.5%   0.001869733 0.01746494 -0.4510364 -0.78336887
# 97.5%  0.340064849 0.41900361  0.3574609  0.09506687






#### Inverse probability of censoring weighting
ds4$year = ds4$year - 1987;
mod.ipcw = psm(Surv(Y, censor) ~ A + cage + cbmi + ctabac + meno + grade + staging_rvalue +
                        recod_oest + recod_prog + type_chx + histfm_1degre + prise_hts +
                        ind_chi + ind_rad + ind_herceptin + rcs(year), data = ds4, dist = "weibull");
# BIC(mod.ipcw);
# gaussian = 8625.866
# loglogistic = 5690.438
# weibull =  5376.942
# exponential = 26627.41
n = nrow(ds4);
predC = numeric(n);
for(j in (1:n)[ds4$censor == 0]){
  dati = ds4[j,];
  predC[j] = 1 - survest(mod.ipcw, newdata = dati, times = ds4$Y[j], se.fit = FALSE)$surv;
}
ipcw = 1/(1 - predC);

#### Naive dWOLS
naive = lm(log(Y) ~ cage + oes_bmi + ctabac + meno + grade + staging_rvalue +
                    recod_prog + type_chx + histfm_1degre + prise_hts +
                    ind_chi + ind_rad + ind_herceptin + rcs(year) +
                    A*(1 + oes_bmi),
           weights = iptw*ipcw, data = ds4, subset = censor == 0);
summary(naive);
# A                       0.078460   0.070440   1.114  0.26555    
# oes_bmioest=1,cbmi=1:A  0.050455   0.102631   0.492  0.62308    
# oes_bmioest=2,cbmi=0:A -0.313148   0.156782  -1.997  0.04600 *  
# oes_bmioest=2,cbmi=1:A -0.419895   0.158523  -2.649  0.00818 ** 

#### CE dWOLS
PATS1 = lm(log(Y) ~ cage + oes_bmi + ctabac + meno + grade + staging_rvalue +
                    recod_prog + type_chx + histfm_1degre + prise_hts +
                    ind_chi + ind_rad + ind_herceptin + rcs(year) +
                    A*(1 + oes_bmi + recod_prog),
           weights = iptw*ipcw, data = ds4, subset = censor == 0);
ds4$Q = with(ds4, model.matrix(Y ~ oes_bmi + recod_prog))%*%tail(coef(PATS1), 6);
lm(Q ~ oes_bmi, data = ds4);
# Coefficients:
#          (Intercept)  oes_bmioest=1,cbmi=1  oes_bmioest=2,cbmi=0  oes_bmioest=2,cbmi=1  
#              0.06708               0.02997              -0.27357              -0.39651  



#### Confidence intervals w/bootstrap
bootf = function(ds, index){
  ds4 = ds[index,];


  # Inverse probability of treatment weighting
  mod.iptw = glm(A ~ cage + cbmi + ctabac + meno + grade + staging_rvalue +
                 recod_oest + recod_prog + type_chx + histfm_1degre + prise_hts +
                 ind_chi + ind_rad + ind_herceptin + rcs(year),
                 data = ds4, family = "binomial");
  iptw = abs(ds4$A - predict(mod.iptw, type = "res"));

  mod.ipcw = psm(Surv(Y, censor) ~ A + cage + cbmi + ctabac + meno + grade + staging_rvalue +
                                   recod_oest + recod_prog + type_chx + histfm_1degre + prise_hts +
                                   ind_chi + ind_rad + ind_herceptin + rcs(year), data = ds4, dist = "weibull");
  n = nrow(ds4);
  predC = numeric(n);
  for(j in (1:n)[ds4$censor == 0]){
    dati = ds4[j,];
    predC[j] = 1 - survest(mod.ipcw, newdata = dati, times = ds4$Y[j], se.fit = FALSE)$surv;
  }
  ipcw = 1/(1 - predC);



  # Naive dWOLS
  ds4$oes_bmi = paste0("oest=",ds4$recod_oest,",cbmi=",ds4$cbmi);

  naive = lm(log(Y) ~ cage + oes_bmi + ctabac + meno + grade + staging_rvalue +
                      recod_prog + type_chx + histfm_1degre + prise_hts +
                      ind_chi + ind_rad + ind_herceptin + rcs(year) +
                      A*(1 + oes_bmi),
             weights = iptw*ipcw, data = ds4, subset = censor == 0);
  psi.naive = tail(coef(naive),4);


  # CE dWOLS
  PATS1 = lm(log(Y) ~ cage + oes_bmi + ctabac + meno + grade + staging_rvalue +
                      recod_prog + type_chx + histfm_1degre + prise_hts +
                      ind_chi + ind_rad + ind_herceptin + rcs(year) +
                      A*(1 + oes_bmi + recod_prog),
             weights = iptw*ipcw, data = ds4, subset = censor == 0);
  ds4$Q = with(ds4, model.matrix(Y ~ oes_bmi + recod_prog))%*%tail(coef(PATS1), 6);
  psi.PATS1 = coef(lm(Q ~ oes_bmi, data = ds4));
  results = c(psi.naive, psi.PATS1);
  return(results);
}
set.seed(7481971);
bootres = boot(ds4, bootf, 1000);

psi = cbind(bootres$t[,1], bootres$t[,1] + bootres$t[,2],  bootres$t[,1] + bootres$t[,3],  bootres$t[,1] + bootres$t[,4],
            bootres$t[,5], bootres$t[,5] + bootres$t[,6],  bootres$t[,5] + bootres$t[,7],  bootres$t[,5] + bootres$t[,8]);

cbind(bootres$t0[1], bootres$t0[1] + bootres$t0[2],  bootres$t0[1] + bootres$t0[3],  bootres$t0[1] + bootres$t0[4],
      bootres$t0[5], bootres$t0[5] + bootres$t0[6],  bootres$t0[5] + bootres$t0[7],  bootres$t0[5] + bootres$t0[8]);
apply(psi, 2, quantile, c(0.025, 0.975));


#              [,1]        [,2]       [,3]       [,4]       [,5]       [,6]       [,7]       [,8]
# A      0.07845994   0.1289148 -0.2346885 -0.3414352 0.06708429 0.09705687 -0.2064854 -0.3294241

#              [,1]        [,2]       [,3]        [,4]        [,5]       [,6]       [,7]        [,8]
# 2.5%  -0.07267342 -0.08220434 -0.6418125 -0.67374376 -0.08579715 -0.1045635 -0.6088102 -0.65269214
# 97.5%  0.24231790  0.34398695  0.1348076  0.02691033  0.22121062  0.2970119  0.1884321  0.03437435

















