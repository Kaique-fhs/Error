#' @title Regression in model with errors
#' @name MEM
#'
#' @description MEM is used to adjust models with errors in the variables.    It can be used to perform regressions on variables that have measurement errors.
#'
#' @param x Values vector
#' @param y Values vector
#' @param lambda_x Ratio of variances between values and error
#' @param conf.level Significance level for the hypothesis test
#' @param Correction Small sample correction, use TRUE or FALSE
#' @param method The correction method that will be considered "Bartlett", "B1" or "B2"
#' @param ... Not always use all arguments
#'
#' @details Only use the conf.level, beta_til e method arguments when running equals "Y", otherwise these arguments will not affect the results.
#'    
#' @return Estimated Parameters \code{x} e \code{y}.
#'
#' @author Kaíque Ferreira Henrique de Souza
#' @author Tatiane Ferreira do N. M. da Silva
#' 
#' @seealso \code{\link[stats]{lm}}
#'
#' @examples
#' x = seq(1,100,length.out = 100)
#' y = sort(runif(100,1,10))
#' 
#' MEM(x, y, lambda_x = 2 )
#' MEM(x, y, lambda_x = 1, conf.level = 0.95 ,Correction = TRUE)
#' 
#' @importFrom stats qchisq
#' 
#' @references Bartlett, M.S. (1937), Properties of sufficiency and statistical tests,Proceedings of RoyalSociety of London A, 160, 268-282.
#' @references Fuller, S. (1987),Measurement Error Models. Wiley, New York.
#' @references Melo, T.F.N., Vasconcellos, K.L.P., Lemonte, A.J. (2009). Some restriction tests in a newclass of regression models for proportions.Computational Statistics and Data Analysis, 53,3972–3979.
#'
#' @export
MEM = function(x, y, lambda_x , Correction = FALSE, method = "Bartlett", conf.level = 0.95 , ...){
 
  
  
  # Operating Conditions
  if(Correction == TRUE){
  
  if(is.vector(x) & is.numeric(x)){
    if(is.vector(y) & is.numeric(y)){
      if(is.numeric(lambda_x) & (length(lambda_x)==1)){
        
        
        ##### Calculating the sums of squares
        x_barra = mean(x)
        y_barra = mean(y)
        
        beta_hat      = c()
        alfa_hat      = c()
        mu_x_hat      = c()
        sigma2_e_hat  = c()
        sigma2_u_hat  = c()
        Det_Sigma_hat = c()
        
        alfa_til      = c()
        mu_x_til      = c()
        sigma2_e_til  = c()
        sigma2_u_til  = c()
        Det_Sigma_til = c()
        beta_til      = 0
        
        mu_hat = array(NA,c(2,1))
        mu_til = array(NA,c(2,1))
        
        Sigma_hat     = array(NA,c(2,2))
        Inv_Sigma_hat = array(NA,c(2,2))
        Det_Sigma_hat = c()
        
        Sigma_til     = array(NA,c(2,2))
        Inv_Sigma_til = array(NA,c(2,2))
        Det_Sigma_til = c()
        
        L_teta_hat = c()
        L_teta_til = c()
        
        n = 1
        Z             = array(NA,c(2,1,n))
        Z_mean        = matrix(NA,2,1)
        diferenca_hat = array(NA,c(2,1,n))
        diferenca_til = array(NA,c(2,1,n))
        Y_um_dois     = matrix(NA,c(2,1))
        Y_zero        = c()
        x_til         = c()
        v_hat = c()
        
        
        Sigma      = list()
        Sigma[[1]] = array(0,c(2,2,5))
        Sigma[[2]] = array(0,c(2,2,5,5))
        Sigma[[3]] = array(0,c(2,2,5,5,5))
        Sigma[[4]] = array(0,c(2,2,5,5,5,5))
        
        Mu      = list()
        Mu[[1]] = array(0,c(2,1,5))
        Mu[[2]] = array(0,c(2,1,5,5))
        Mu[[3]] = array(0,c(2,1,5,5,5))
        Mu[[4]] = array(0,c(2,1,5,5,5,5))
        
        D_Inv_Sigma  = list()
        D_Inv_Sigma[[1]] = array(0,c(2,2,5))
        D_Inv_Sigma[[2]] = array(0,c(2,2,5,5))
        D_Inv_Sigma[[3]] = array(0,c(2,2,5,5,5))
        D_Inv_Sigma[[4]] = array(0,c(2,2,5,5,5,5)) 
        
        K_twrs  = array(0,c(5,5,5,5))
        K_twr   = array(0,c(5,5,5,5))
        K_twr_s = array(0,c(5,5,5,5))
        K_tr_ws = array(0,c(5,5,5,5))
        K_ws_t  = array(0,c(5,5,5,5))
        
        A = array(0,c(5,5,5,5))
        P = array(0,c(5,5,5))
        Q = array(0,c(5,5,5))
        
        K     = array(NA,c(5,5))
        L     = array(0,c(5,5))
        M     = array(0,c(5,5))
        N     = array(0,c(5,5))
        Inv_K = array(NA,c(5,5))
        
        Epsilon_p   = c()
        Epsilon_p_1 = c()
        
        K_p_1     = array(NA,c(4,4))
        L_p_1     = array(NA,c(4,4))
        M_p_1     = array(NA,c(4,4))
        N_p_1     = array(NA,c(4,4))
        Inv_K_p_1 = array(NA,c(4,4))
        
        CB1 = c()
        CB2 = c()
        CB3 = c()
        
        Produto_hat = matrix(0,nrow=n,ncol=1)
        Produto_til = matrix(0,nrow=n,ncol=1) 
        
        ##### Calculando as somas de quadrados
        S_xx    = mean((x[] - x_barra)^2)
        S_yy    = mean((y[] - y_barra)^2)
        S_xy    = mean((x[] - x_barra)*(y[] - y_barra))
        M_zz    = ((1/(length(x)-1))*sum(((y[] - y_barra)^2)+((x[] - x_barra)^2)))
        M_zx    = ((1/(length(x)-1))*sum(((y[] - y_barra)*(x[] - x_barra))+((x[] - x_barra)^2)))
        M_yy    = ((1/(length(x)-1))*sum((y[]-y_barra)^2))
        M_x    = ((1/(length(x)-1))*sum((x[]-x_barra)^2))
        M_xy    = ((1/(length(x)-1))*sum((x[]-x_barra)*(y[]-y_barra)))
        
        
        ## Maximum likelihood estimates (e.m.v.)
        beta_hat     = ((lambda_x + 1)*S_xy)/(lambda_x*S_xx)
      #  beta_um_hat  = (M_xx - )
        mu_x_hat     = x_barra
        alfa_hat     = y_barra - (beta_hat*x_barra)
        sigma2_u_hat = S_xx/(lambda_x + 1)
        sigma2_e_hat = (lambda_x*(S_yy*S_xx - S_xy^2) - S_xy^2)/(lambda_x*S_xx)
    #    sigma_xx     = (M_xx - sigma_uu)
    #    sigma_uu     = (M_xx - sigma_xx)
        Y_um_dois    = ((1/M_zz)*(M_zx-(matrix(c(0,sigma2_u_hat),2,1))))
        Y_zero       = (((1-Y_um_dois[2,1])*x_barra)-(Y_um_dois[1,1]*y_barra))
        i = 1:length(x)
        x_til        = (Y_zero + (Y_um_dois[1,1]*y[i]) + (Y_um_dois[2,1]*x[i])) 
        v_hat        = ( y[i] - 67.56 -(x[i]*0.4232) )
        x_hat        = (((beta_hat*sigma2_u_hat*(y[]-alfa_hat)) + (sigma2_e_hat*x[])) / ((beta_hat*beta_hat*sigma2_u_hat) + sigma2_e_hat))
        
        ##Estimativas de máxima verossimilhança (e.m.v.) sob H0
        alfa_til     = y_barra - (beta_til *x_barra)
        mu_x_til     = x_barra;
        sigma2_u_til = S_xx/(lambda_x + 1);
        sigma2_e_til = ((lambda_x-1)*lambda_x*beta_til^2*S_xx - 2*(lambda_x+1)*lambda_x*beta_til*S_xy +	(lambda_x+1)^2*S_yy)/((lambda_x+1)^2)
        
        ##### Matrix de Médias
        mu_hat[,] = matrix(c(alfa_hat + (beta_hat*mu_x_hat), mu_x_hat), ncol=1, nrow=2)
        mu_til[,] = matrix(c(alfa_til + (beta_til    *mu_x_til), mu_x_til), ncol=1, nrow=2)
        
        
        ##### Matrix de Covariãncias
        Sigma_hat[,] = matrix(c((sigma2_e_hat)+(beta_hat^2*lambda_x*sigma2_u_hat),
                                beta_hat*lambda_x*(sigma2_u_hat),
                                beta_hat*lambda_x*(sigma2_u_hat),
                                ((lambda_x+1)*sigma2_u_hat)), ncol=2,nrow=2)
        
        
        Inv_Sigma_hat[,] = solve(Sigma_hat[,])
        Det_Sigma_hat= det(Sigma_hat[,])
        
        Sigma_til[,] = matrix(c((sigma2_e_til)+(beta_til^2*lambda_x*sigma2_u_til),
                                beta_til *lambda_x*(sigma2_u_til),
                                beta_til *lambda_x*(sigma2_u_til),
                                ((lambda_x+1)*sigma2_u_til)), ncol=2,nrow=2)
        
        Inv_Sigma_til[,] = solve(Sigma_til[,])
        Det_Sigma_til= det(Sigma_til[,])
        
        ##### Preenchendo os vetores Z, diferença e produto
        for (i in 1:n)
        {
          Z[,,i]             = matrix(c(y[][i],x[][i]),ncol = 1,nrow = 2)
          
          diferenca_hat[,,i] = matrix(Z[,,i] - mu_hat[,],ncol = 1)       
          diferenca_til[,,i] = matrix(Z[,,i] - mu_til[,],ncol = 1)     
          
          Produto_hat[i]     = (t(matrix(diferenca_hat[,,i],ncol = 1)) %*% Inv_Sigma_hat[,] %*% matrix(diferenca_hat[,,i],ncol = 1)) 
          Produto_til[i]     = (t(matrix(diferenca_til[,,i],ncol = 1)) %*% Inv_Sigma_til[,] %*% matrix(diferenca_til[,,i],ncol = 1)) 
        }  
        
        
        ##### Log da verossimilhança avaliados nas estimativas de máxima 
        ##### verossimilhança irrestrita e restrita
        
        L_teta_hat = - n*log(2*pi) - ((n/2)*log(Det_Sigma_hat)) - ((1/2)*sum(Produto_hat))    
        L_teta_til = - n*log(2*pi) - ((n/2)*log(Det_Sigma_til)) - ((1/2)*sum(Produto_til)) 
        
        ##### Primeira derivada do vetor mu      
        Mu[[1]][1,1,1] = 1                  
        Mu[[1]][1,1,2] = mu_x_til       
        Mu[[1]][1,1,3] = beta_til           
        Mu[[1]][2,1,3] = 1                  
        
        
        ##### Segunda derivada do vetor mu 
        Mu[[2]][1,1,2,3] = Mu[[2]][1,1,3,2] = 1
        
        ##### Primeira derivada da matriz Sigma
        Sigma[[1]][1,1,2] = 2*lambda_x*sigma2_u_til*beta_til
        Sigma[[1]][1,2,1] = Sigma[[1]][2,1,1] = lambda_x*sigma2_u_til
        Sigma[[1]][1,1,4] = 1
        Sigma[[1]][1,1,5] = beta_til^2*lambda_x
        Sigma[[1]][2,1,5] = Sigma[[1]][1,2,5] = beta_til*lambda_x
        Sigma[[1]][2,2,5] = lambda_x+1
        
        ##### Segunda derivada da matriz Sigma
        Sigma[[2]][1,1,2,2]   = 2*lambda_x*sigma2_u_til
        Sigma[[2]][1,1,5,2]   = Sigma[[2]][1,1,2,5] = 2*lambda_x*beta_til
        Sigma[[2]][1,2,5,2]   = Sigma[[2]][2,1,5,2] =   Sigma[[2]][1,2,2,5]   = Sigma[[2]][2,1,2,5] = lambda_x
        
        ##### Terceira derivada da matriz Sigma
        Sigma[[3]][1,1,2,2,5] = Sigma[[3]][1,1,2,5,2] = Sigma[[3]][1,1,5,2,2] = 2*lambda_x
        
        ##### Primeira derivada da inversa de Sigma
        ##D_Inv_Sigma[[1]]
        for (r  in c(1:5)) 
        {D_Inv_Sigma[[1]][,,r]=(-(Inv_Sigma_til[,] %*% Sigma[[1]][,,r] %*% Inv_Sigma_til[,]))}
        
        
        ##### Segunda derivada da inversa de Sigma
        ##D_Inv_Sigma[[2]]
        for (r in 1:5)  {
          for (s in 1:5)  {
            D_Inv_Sigma[[2]][,,r,s] = (-(D_Inv_Sigma[[1]][,,s] %*% Sigma[[1]][,,r]   %*% Inv_Sigma_til[,])
                                       -(Inv_Sigma_til[,]      %*% Sigma[[2]][,,r,s] %*% Inv_Sigma_til[,])
                                       -(Inv_Sigma_til[,]      %*% Sigma[[1]][,,r]   %*% D_Inv_Sigma[[1]][,,s]))
          }
        }
        
        ##### Terceira derivada da inversa de Sigma
        ##D_Inv_Sigma[[3]]
        for (t in 1:5) 
        {
          for (s in 1:5)
          {
            for (r in 1:5) 
            {
              D_Inv_Sigma[[3]][,,r,s,t] = (-(D_Inv_Sigma[[2]][,,s,t]  %*% Sigma[[1]][,,r]    %*% Inv_Sigma_til[,])
                                           -(D_Inv_Sigma[[1]][,,s]    %*% Sigma[[2]][,,r,t]  %*% Inv_Sigma_til[,])
                                           -(D_Inv_Sigma[[1]][,,s]    %*% Sigma[[1]][,,r]    %*% D_Inv_Sigma[[1]][,,t])
                                           -(D_Inv_Sigma[[1]][,,t]    %*% Sigma[[2]][,,r,s]  %*% Inv_Sigma_til[,])
                                           -(Inv_Sigma_til[,]         %*% Sigma[[3]][,,r,s,t]%*% Inv_Sigma_til[,])
                                           -(Inv_Sigma_til[,]         %*% Sigma[[2]][,,r,s]  %*% D_Inv_Sigma[[1]][,,t])
                                           -(D_Inv_Sigma[[1]][,,t]    %*% Sigma[[1]][,,r]    %*% D_Inv_Sigma[[1]][,,s])
                                           -(Inv_Sigma_til[,]         %*% Sigma[[2]][,,r,t]  %*% D_Inv_Sigma[[1]][,,s])
                                           -(Inv_Sigma_til[,]         %*% Sigma[[1]][,,r]    %*% D_Inv_Sigma[[2]][,,s,t]))
            }
          }
        }
        
        
        ##### Quarta derivada da inversa de Sigma
        ##D_Inv_Sigma[[4]]
        for (r in 1:5)  {
          for (s in 1:5)   {
            for (t in 1:5)  {
              for (w in 1:5)   {
                D_Inv_Sigma[[4]][,,r,s,t,w]=(-(D_Inv_Sigma[[3]][,,s,t,w] %*% Sigma[[1]][,,r]       %*% Inv_Sigma_til[,])
                                             -(D_Inv_Sigma[[2]][,,s,t]   %*% Sigma[[2]][,,r,w]     %*% Inv_Sigma_til[,])
                                             -(D_Inv_Sigma[[2]][,,s,t]   %*% Sigma[[1]][,,r]       %*% D_Inv_Sigma[[1]][,,w])
                                             -(D_Inv_Sigma[[2]][,,s,w]   %*% Sigma[[2]][,,r,t]     %*% Inv_Sigma_til[,])
                                             -(D_Inv_Sigma[[1]][,,s]     %*% Sigma[[3]][,,r,t,w]   %*% Inv_Sigma_til[,])
                                             -(D_Inv_Sigma[[1]][,,s]     %*% Sigma[[2]][,,r,t]     %*% D_Inv_Sigma[[1]][,,w])
                                             -(D_Inv_Sigma[[2]][,,s,w]   %*% Sigma[[1]][,,r]       %*% D_Inv_Sigma[[1]][,,t])
                                             -(D_Inv_Sigma[[1]][,,s]     %*% Sigma[[2]][,,r,w]     %*% D_Inv_Sigma[[1]][,,t])
                                             -(D_Inv_Sigma[[1]][,,s]     %*% Sigma[[1]][,,r]       %*% D_Inv_Sigma[[2]][,,t,w])
                                             -(D_Inv_Sigma[[2]][,,t,w]   %*% Sigma[[2]][,,r,s]     %*% Inv_Sigma_til[,])
                                             -(D_Inv_Sigma[[1]][,,t]     %*% Sigma[[3]][,,r,s,w]   %*% Inv_Sigma_til[,])
                                             -(D_Inv_Sigma[[1]][,,t]     %*% Sigma[[2]][,,r,s]     %*% D_Inv_Sigma[[1]][,,w])
                                             -(D_Inv_Sigma[[1]][,,w]     %*% Sigma[[3]][,,r,s,t]   %*% Inv_Sigma_til[,])
                                             -(Inv_Sigma_til[,]          %*% Sigma[[3]][,,r,s,t]   %*% D_Inv_Sigma[[1]][,,w])
                                             -(D_Inv_Sigma[[1]][,,w]     %*% Sigma[[2]][,,r,s]     %*% D_Inv_Sigma[[1]][,,t])
                                             -(Inv_Sigma_til[,]          %*% Sigma[[3]][,,r,s,w]   %*% D_Inv_Sigma[[1]][,,t])
                                             -(Inv_Sigma_til[,]          %*% Sigma[[2]][,,r,s]     %*% D_Inv_Sigma[[2]][,,t,w])
                                             -(D_Inv_Sigma[[2]][,,t,w]   %*% Sigma[[1]][,,r]       %*% D_Inv_Sigma[[1]][,,s])
                                             -(D_Inv_Sigma[[1]][,,t]     %*% Sigma[[2]][,,r,w]     %*% D_Inv_Sigma[[1]][,,s])
                                             -(D_Inv_Sigma[[1]][,,t]     %*% Sigma[[1]][,,r]       %*% D_Inv_Sigma[[2]][,,s,w])
                                             -(D_Inv_Sigma[[1]][,,w]     %*% Sigma[[2]][,,r,t]     %*% D_Inv_Sigma[[1]][,,s])
                                             -(Inv_Sigma_til[,]          %*% Sigma[[3]][,,r,t,w]   %*% D_Inv_Sigma[[1]][,,s])
                                             -(Inv_Sigma_til[,]          %*% Sigma[[2]][,,r,t]     %*% D_Inv_Sigma[[2]][,,s,w])
                                             -(D_Inv_Sigma[[1]][,,w]     %*% Sigma[[1]][,,r]       %*% D_Inv_Sigma[[2]][,,s,t])
                                             -(Inv_Sigma_til[,]          %*% Sigma[[2]][,,r,w]     %*% D_Inv_Sigma[[2]][,,s,t])
                                             -(Inv_Sigma_til[,]          %*% Sigma[[1]][,,r]       %*% D_Inv_Sigma[[3]][,,s,t,w]))
              }
            }
          }
        }
        
        ### Função para calcular o traço de uma matriz
        traco = function(M){sum(diag(M))}  
        
        for (t in 1:5)  {
          for (w in 1:5)  {
            ##### Matriz de Informação de Fisher K
            
            K[t,w] = (-(n/2)*traco(Sigma[[1]][,,t] %*% D_Inv_Sigma[[1]][,,w])+n*t(Mu[[1]][,,w]) %*% Inv_Sigma_til[,] %*% Mu[[1]][,,t])
            
            
            for (s in 1:5)  {
              
              ### Q^s = {K_{ws}^t} 
              Q[t,w,s] = (
                (n/2)*traco(
                  Sigma[[2]][,,s,t] %*% D_Inv_Sigma[[1]][,,w] + Sigma[[1]][,,s] %*% D_Inv_Sigma[[2]][,,w,t]
                )
                -n*
                  (
                    t(Mu[[2]][,,s,t]) %*% Inv_Sigma_til[,]      %*% Mu[[1]][,,w]
                    + t(Mu[[1]][,,s])   %*% D_Inv_Sigma[[1]][,,t] %*% Mu[[1]][,,w]
                    + t(Mu[[1]][,,s])   %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,w,t]
                  )
              )
              
              
              ### P^s = {K_{tws}}
              P[t,w,s] = (
                -(n/2)*traco(
                  D_Inv_Sigma[[2]][,,w,t]   %*% Sigma[[1]][,,s]
                  + D_Inv_Sigma[[1]][,,w]     %*% Sigma[[2]][,,s,t]
                  + D_Inv_Sigma[[1]][,,t]     %*% Sigma[[2]][,,s,w]
                  + Inv_Sigma_til[,]          %*% Sigma[[3]][,,s,w,t]
                  + D_Inv_Sigma[[3]][,,s,w,t] %*% Sigma_til[,]
                )
                -n*( 
                  t(Mu[[2]][,,s,w]) %*% Inv_Sigma_til[,]      %*% Mu[[1]][,,t]
                  + t(Mu[[2]][,,s,t]) %*% Inv_Sigma_til[,]      %*% Mu[[1]][,,w]
                  + t(Mu[[1]][,,s])   %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,w,t]
                  + t(Mu[[1]][,,s])   %*% D_Inv_Sigma[[1]][,,w] %*% Mu[[1]][,,t]
                  + t(Mu[[1]][,,w])   %*% D_Inv_Sigma[[1]][,,s] %*% Mu[[1]][,,t]
                  + t(Mu[[1]][,,s])   %*% D_Inv_Sigma[[1]][,,t] %*% Mu[[1]][,,w]
                )
              )
              
              
              for (r in 1:5)  {
                ### A^{rs} = {A_{tw}^{rs}}
                A[t,w,r,s] = (
                  -(n/8)*traco(
                    D_Inv_Sigma[[3]][,,r,w,t]   %*% Sigma[[1]][,,s]
                    +   D_Inv_Sigma[[2]][,,r,w]     %*% Sigma[[2]][,,s,t]
                    +   D_Inv_Sigma[[2]][,,r,t]     %*% Sigma[[2]][,,s,w]
                    +   D_Inv_Sigma[[2]][,,w,t]     %*% Sigma[[2]][,,s,r]
                    + 4*D_Inv_Sigma[[2]][,,w,s]     %*% Sigma[[2]][,,r,t]
                    + 5*D_Inv_Sigma[[1]][,,w]       %*% Sigma[[3]][,,s,r,t]
                    +   D_Inv_Sigma[[1]][,,t]       %*% Sigma[[3]][,,s,r,w]
                    +   D_Inv_Sigma[[1]][,,r]       %*% Sigma[[3]][,,s,w,t]
                    +   D_Inv_Sigma[[4]][,,s,r,w,t] %*% Sigma_til[,]
                    + 8*D_Inv_Sigma[[2]][,,w,s]     %*% Sigma[[1]][,,r]   %*% D_Inv_Sigma[[1]][,,t]   %*% Sigma_til[,]
                    + 8*D_Inv_Sigma[[1]][,,w]       %*% Sigma[[2]][,,r,s] %*% D_Inv_Sigma[[1]][,,t]   %*% Sigma_til[,]
                    + 8*D_Inv_Sigma[[1]][,,w]       %*% Sigma[[1]][,,r]   %*% D_Inv_Sigma[[2]][,,t,s] %*% Sigma_til[,]
                    + 8*D_Inv_Sigma[[1]][,,w]       %*% Sigma[[1]][,,r]   %*% D_Inv_Sigma[[1]][,,t]   %*% Sigma[[1]][,,r]
                  )
                  -(n/4)*(
                    t(Mu[[2]][,,r,w]) %*% D_Inv_Sigma[[1]][,,s]   %*% Mu[[1]][,,t]
                    +   t(Mu[[2]][,,r,w]) %*% D_Inv_Sigma[[1]][,,t]   %*% Mu[[1]][,,s]
                    - 3*t(Mu[[2]][,,r,t]) %*% D_Inv_Sigma[[1]][,,s]   %*% Mu[[1]][,,w]
                    +   t(Mu[[2]][,,r,t]) %*% D_Inv_Sigma[[1]][,,w]   %*% Mu[[1]][,,s]
                    +   t(Mu[[2]][,,s,r]) %*% D_Inv_Sigma[[1]][,,w]   %*% Mu[[1]][,,t]
                    - 3*t(Mu[[2]][,,s,r]) %*% D_Inv_Sigma[[1]][,,t]   %*% Mu[[1]][,,w]
                    - 3*t(Mu[[2]][,,s,w]) %*% D_Inv_Sigma[[1]][,,t]   %*% Mu[[1]][,,r]
                    - 3*t(Mu[[2]][,,s,w]) %*% D_Inv_Sigma[[1]][,,r]   %*% Mu[[1]][,,t]
                    - 3*t(Mu[[2]][,,s,t]) %*% D_Inv_Sigma[[1]][,,r]   %*% Mu[[1]][,,w]
                    +   t(Mu[[2]][,,s,t]) %*% D_Inv_Sigma[[1]][,,w]   %*% Mu[[1]][,,r]
                    +   t(Mu[[2]][,,w,t]) %*% D_Inv_Sigma[[1]][,,s]   %*% Mu[[1]][,,r]
                    +   t(Mu[[2]][,,w,t]) %*% D_Inv_Sigma[[1]][,,r]   %*% Mu[[1]][,,s]
                    - 3*t(Mu[[1]][,,w])   %*% D_Inv_Sigma[[2]][,,s,r] %*% Mu[[1]][,,t]
                    +   t(Mu[[1]][,,r])   %*% D_Inv_Sigma[[2]][,,s,w] %*% Mu[[1]][,,t]
                    - 3*t(Mu[[1]][,,r])   %*% D_Inv_Sigma[[2]][,,s,t] %*% Mu[[1]][,,w]
                    +   t(Mu[[1]][,,s])   %*% D_Inv_Sigma[[2]][,,r,t] %*% Mu[[1]][,,w]
                    +   t(Mu[[1]][,,s])   %*% D_Inv_Sigma[[2]][,,r,w] %*% Mu[[1]][,,t]
                    +   t(Mu[[1]][,,s])   %*% D_Inv_Sigma[[2]][,,w,t] %*% Mu[[1]][,,r]
                    +   t(Mu[[2]][,,s,r]) %*% Inv_Sigma_til[,]        %*% Mu[[2]][,,w,t]
                    - 3*t(Mu[[2]][,,s,w]) %*% Inv_Sigma_til[,]        %*% Mu[[2]][,,r,t]
                    +   t(Mu[[2]][,,s,t]) %*% Inv_Sigma_til[,]        %*% Mu[[2]][,,r,w]
                  )
                )
                
              }## Fim do loop do r            
            }## Fim do loop do s
          }## Fim do loop do w
        } ## Fim do loop do t
        
        Inv_K[,] = (solve(K[,]))
        
        
        for (r in 1:5) {
          for (s in 1:5)  {
            ##### L
            L[r,s] = traco(Inv_K[,] %*% A[,,r,s])
            
            ##### M
            M[r,s] = (
              -(1/6)*traco(Inv_K[,] %*% P[,,r] %*% Inv_K[,] %*%   P[,,s])
              +traco(Inv_K[,] %*% P[,,r] %*% Inv_K[,] %*% t(Q[,,s]))
              -traco(Inv_K[,] %*% Q[,,r] %*% Inv_K[,] %*%   Q[,,s])
            )
            
            ##### N
            N[r,s] = (
              -(1/4)*traco(P[,,r] %*% Inv_K[,])*traco(P[,,s] %*% Inv_K[,])
              +traco(P[,,r] %*% Inv_K[,])*traco(Q[,,s] %*% Inv_K[,])
              -traco(Q[,,r] %*% Inv_K[,])*traco(Q[,,s] %*% Inv_K[,])
            )
          }
        }
        
        ##### Epsilon(p)
        Epsilon_p = traco(Inv_K[,] %*% (L[,] - M[,] - N[,]))
        
        ##### Cálculo de Epsilon(p-1): 
        ##### Retiramos a 2ª coluna e 2ª linha das matrizes, pois teta = (beta,alfa,mu_x,sigma2_e,sigma2_u)
        a = 2
        K_p_1[,] = K[-a,-a]
        L_p_1[,] = L[-a,-a]
        M_p_1[,] = M[-a,-a]
        N_p_1[,] = N[-a,-a]
        
        Inv_K_p_1[,] = solve(K_p_1[,])
        
        Epsilon_p_1 = traco(Inv_K_p_1[,] %*% (L_p_1[,] - M_p_1[,] - N_p_1[,]))
        
        ##### Constante C (Correção de Bartlett)
        CB1 = 1 + (Epsilon_p - Epsilon_p_1)
        CB2 = exp(1 - CB1)
        CB3 = 1 - (Epsilon_p - Epsilon_p_1)
        
        ##### Estatística da razão de verossimilhanças
        LR = 2*(L_teta_hat - L_teta_til)
        
        LR_B1 = LR/CB1 # Method Bartlett
        LR_B2 = LR*CB2 # Method B2
        LR_B3 = LR*CB3 # Method B3
        
        
        if(          method == "Bartlett"){ Lrc = LR_B1 
           }else if( method == "B2"      ){ Lrc = LR_B2
            }else if(method == "B3"      ){ Lrc = LR_B3
             }else{ warning("Argument only accepts 'Bartlett' or 'B2' or 'B3' options")}
        
        ##### Impressão de resultados - Parámetros estimados
        Parametros  = c("Beta","Alfa","Mu_x","Sigma2_e","Sigma2_u")
        EMV         = c(round(beta_hat,2),round(alfa_hat,2),round(mu_x_hat,2),round(sigma2_e_hat,2),round(sigma2_u_hat,2))
        EMV_H0      = c(round(beta_til,2),round(alfa_til,2),round(mu_x_til,2),round(sigma2_e_til,2),round(sigma2_u_til,2))
        
        # Function exit
        Saida = function(){
          cat("Call: \n")
          cat("MEV( x, y ", lambda_x," \n")
          cat("\n")
          cat("Residuals: \n")
          print( summary(x_hat) )
          cat("\n")
          cat("Coefficients: \n")
          
          cat(  " XXXXX "  )
          cat("\n")
          cat("--- \n")
          cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 \n")
          cat("\n")
          cat("Residual standard error:", " XXXXX " ," on", " XXXXX " ,"degrees of freedom \n")
          cat("Multiple R-squared:", " XXXXX ",	"Adjusted R-squared:",  " XXXXX ", "\n")
          cat("F-statistic:", " XXXXX " ,"on "," XXXXX "," and", " XXXXX "," DF", " p-value: "," XXXXX ", " \n")
        } #End of function exit list
        
        Saida()
        
        
      }else{ warning("lambda_x must be a single numbers") }
    }else{ warning("Y must be a vector of numbers") }
  }else{ warning("X must be a vector of numbers") }
    
  }else if(Correction == FALSE){
    
    # Operating Conditions
    if(is.vector(x) & is.numeric(x)){
      if(is.vector(y) & is.numeric(y)){
        if(is.numeric(lambda_x) & (length(lambda_x)==1)){
          
          ##### Calculating the sums of squares
          x_barra = mean(x)
          y_barra = mean(y)
          S_xx    = mean((x - x_barra)^2)
          S_yy    = mean((y - y_barra)^2)
          S_xy    = mean((x - x_barra)*(y - y_barra))
          
          ## Maximum likelihood estimates (e.m.v.)
          beta_hat     = ((lambda_x + 1)*S_xy)/(lambda_x*S_xx)
          alfa_hat     = y_barra - (beta_hat*x_barra)
          sigma2_u_hat = S_xx/(lambda_x + 1)
          sigma2_e_hat = (lambda_x*(S_yy*S_xx - S_xy^2) - S_xy^2)/(lambda_x*S_xx)
          x_hat        = (((beta_hat*sigma2_u_hat*(y[]-alfa_hat)) + (sigma2_e_hat*x[])) / ((beta_hat*beta_hat*sigma2_u_hat) + sigma2_e_hat))
          
        
          # Function exit
          Saida = function(){
            cat("Call: \n")
            cat("MEV( x, y ", lambda_x," \n")
            cat("\n")
            cat("Residuals: \n")
            print( summary(x_hat) )
            cat("\n")
            cat("Coefficients: \n")
            
            cat(  " XXXXX "  )
            cat("\n")
            cat("--- \n")
            cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 \n")
            cat("\n")
            cat("Residual standard error:", " XXXXX " ," on", " XXXXX " ,"degrees of freedom \n")
            cat("Multiple R-squared:", " XXXXX ",	"Adjusted R-squared:",  " XXXXX ", "\n")
            cat("F-statistic:", " XXXXX " ,"on "," XXXXX "," and", " XXXXX "," DF", " p-value: "," XXXXX ", " \n")
          } #End of function exit list
          
          Saida()
          
          
        }else{ warning("lambda_x must be a single numbers") }
      }else{ warning("Y must be a vector of numbers") }
    }else{ warning("X must be a vector of numbers") }
    
    # Else do Correction
  }else{ warning("The argument thus accepts the options TRUE or FALSE") }
  
} # End of MEM Function
