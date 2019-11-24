# MEV {Error}	Package for R 
Regression in model with errors

# Description
MEV is used to adjust models with errors in the variables. It can be used to perform regressions on variables that have measurement errors.

# Usage
```
MEV(x, y, lambda_x, Correction = FALSE, method = "Bartlett",
  beta_til = 0, conf.level = 0.95, ...)
```

# Arguments
- **x**  Values vector
- **y** Error vector
- **lambda_x** Ratio of variances between values and error
- **Correction** Small sample correction, use TRUE or FALSE
- **method** The correction method that will be considered "Bartlett", "B1" or "B2"
- **beta_til** Parameter value to be tested
- **conf.level** Significance level for the hypothesis test
- **...** Not always use all arguments

# Details
Only use the conf.level, beta_til e method arguments when running equals "Y", otherwise these arguments will not affect the results.

# Author(s)
Kaíque Ferreira Henrique de Souza

Dra. Tatiane Ferreira do N. M. da Silva

# References
Aoki, R., Bolfarine, H. e Singer, J.M. (2001). Null intercept measurement error regressionmodels.TestA, 10, 441-457.

Aubin, E.C.Q. e Cordeiro, G.M. (2000). Bartlett-corrected tests for normal linear modelswhen the error covariance matrix is nonscalar.Communications in Statistics, Theory andMethods,29, 2405-2426.

Barndorff–Nielsen, O.E. (1986), Inference on full or partial parameters, based on the stan-dardized signed log likelihood ratio,Biometrika, 73, 307-322.

Bartlett, M.S. (1937), Properties of sufficiency and statistical tests,Proceedings of RoyalSociety of London A, 160, 268-282.

Bayer, F.M., Cribari–Neto, F. (2013). Bartlett corrections in beta regression models.Journalof Statistical Planning and Inference, 143, 531–547.

Botter, D.A. e Cordeiro G.M. (1997). Bartlett corrections for generalized linear models withdispersion covariates.Communications in Statistics, Theory and Methods, 26, 279-307.

Chan, L.K.&Mak, T.K., On the maximum likelihood estimation of a linear structuralrelashionship when the intercept is known,Journal of Multivariate Analysis, 9, 304-313(1979).

Fuller, S. (1987),Measurement Error Models. Wiley, New York.

Melo, T.F.N., Vasconcellos, K.L.P., Lemonte, A.J. (2009). Some restriction tests in a newclass of regression models for proportions.Computational Statistics and Data Analysis, 53,3972–3979.

# See Also
lm

# Examples

```
x = seq(1,100,length.out = 100)
y = sort(runif(100,1,10))

MEV(x, y, lambda_x = 2 )
MEV(x, y, lambda_x = 1, beta_til = 4 ,conf.level = 0.95 ,Correction = TRUE)
```
