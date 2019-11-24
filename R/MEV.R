# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#Caso Lambda_x conhecido


MEV = function(x, y, lambda_x){

#

##### Calculando as somas de quadrados
x_barra = mean(x)
y_barra = mean(y)
S_xx    = mean((x - x_barra)^2)
S_yy    = mean((y - y_barra)^2)
S_xy    = mean((x - x_barra)*(y - y_barra))

##Estimativas de máxima verossimilhança (e.m.v.)
beta_hat     = ((lambda_x + 1)*S_xy)/(lambda_x*S_xx)
alfa_hat     = y_barra - (beta_hat*x_barra)
sigma2_u_hat = S_xx/(lambda_x + 1)
sigma2_e_hat = (lambda_x*(S_yy*S_xx - S_xy^2) - S_xy^2)/(lambda_x*S_xx)

# Saida da função

Saida = list(
  Coefi= data.frame(Beta     = beta_hat,
                        Alpha    = alfa_hat,
                        Sigma2_u = sigma2_u_hat,
                        Sigma2_e = sigma2_e_hat)
Saida[[2]] = x_barra

 ) Fim da lista de saida da função
} # Fim da função MEV

r=MEV(a,b,3)
str(r)
r
w=lm(b~a)
str(w)
list(a= 3)
