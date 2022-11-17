# PATS
The PATS_function.R file contains an R implementation of the optimal partially adaptive
treatment strategies estimators proposed in Talbot et al (2022) Double robust estimation 
of optimal partially adaptive treatment strategies: An pplication to breast cancer 
treatment using hormonal therapy. Statistics in Medicine. More specifically, the function
implements the IPTW+dWOLS, IPTW+G-est, CE+dWOLS and CE+G-est estimators with an arbitrary
number of time points. However, note that the adaptive m-out-of-n bootstrap is only
available for the two time points setting.

The other files contain the code required to replicate the simulation studies performed
in the aforementioned article.
