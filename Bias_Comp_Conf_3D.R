###################################################
### Programm to derive the bias with competing confounding and plot it in a 3D Plot
### Author: Jost Viebrock
### 04/04/2025
############################

library(scatterplot3d)
library(plotly)
library(akima) 
library(reshape2)
library(latex2exp)
logit <- function(x){
  log(x/(1-x))
}

inv.logit <- function(x){ 
  exp(x)/(1+exp(x))
}




### Plotting function for different prevalence of the competing event
plot_bias_curves_3D_comp.rate <- function( a=0,
                                  c,
                                  f, 
                                  g,
                                  h,
                                  p=1/2,
                                  comp.rate=1,
                                  haz, #baseline hazard for non-treated  
                                  max=50
){
  Stat.time<-Sys.time()
  CIF <- function(k, a, lambda_Y, comp.rate, zeta, gamma, theta, eta, alpha, p) {
    lambda_D=lambda_Y*comp.rate
    # First term calculations
    term1 <- (1 - exp((-lambda_Y * exp(zeta * a) - lambda_D * exp(gamma * a)) * k)) * 
      (lambda_Y * exp(zeta * a) / (lambda_Y * exp(zeta * a) + lambda_D * exp(gamma * a)))
    
    weight1 <- ((1 - p) / 2) / (p*( inv.logit(alpha))^a * ((1 - inv.logit(alpha)))^(1 - a) + (1 - p) / 2)
    
    # Second term calculations
    term2 <- (1 - exp((-lambda_Y * exp(zeta * a + theta) - lambda_D * exp(gamma * a + eta)) * k)) * 
      (lambda_Y * exp(zeta * a + theta) / (lambda_Y * exp(zeta * a + theta) + lambda_D * exp(gamma * a + eta)))
    
    weight2 <- (inv.logit(alpha)^a * (1 - inv.logit(alpha))^(1 - a) * p) / 
      ( p*( inv.logit(alpha))^a * ((1 - inv.logit(alpha)))^(1 - a) + (1 - p) / 2)
    
    # Final result
    result <- term1 * (weight1) + term2 * (weight2)
    return(result)
  }
  
  # Alle combinations from max and alpha 
  combinations <- expand.grid(max = 1:max, comp.rate = exp((1:(999+1)-(999+1)/2)/50))
  

  num_columns <- nrow(combinations)
  
 #matrix 
  plot_matrix <- matrix(0, nrow = 3, ncol = num_columns)
  
  plot_matrix[1, ] <- combinations$max
  plot_matrix[2, ] <- combinations$comp.rate
  

  
  
  for (i in 1:length(plot_matrix[1,])){
    y1 <- CIF(a=1,
              alpha = a,
              lambda_Y = haz,
              comp.rate = plot_matrix[2, i],
              zeta = f,
              gamma = c,
              theta= h,
              eta = g,
              k= plot_matrix[1,i ],
              p= p)
    y0 <- CIF(a=0,
              alpha = a, 
              lambda_Y = haz,
              comp.rate = plot_matrix[2, i],
              zeta = f,
              gamma = c,
              theta= h,
              eta = g,
              k= plot_matrix[1, i],
              p= p)
    # true is randomised (a=0)
    y1true <- CIF(a=1,
                  alpha = 0, 
                  lambda_Y = haz,
                  comp.rate = plot_matrix[2, i],
                  zeta = f,
                  gamma = c,
                  theta= h,
                  eta = g,
                  k= plot_matrix[1, i],
                  p= p)
    y0true <- CIF(a=0,
                  alpha = 0, 
                  lambda_Y = haz,
                  comp.rate = plot_matrix[2, i],
                  zeta = f,
                  gamma = c,
                  theta= h,
                  eta = g,
                  k= plot_matrix[1, i],
                  p= p)
    
    plot_matrix[3,i]<- (y1-y0)-(y1true-y0true)
    
  }
  
 
  grid <- expand.grid(time = unique(plot_matrix[1, ]), alpha = unique(plot_matrix[2, ]))
  grid$bias <- with(grid, apply(grid, 1, function(row) {
    plot_matrix[3, which(plot_matrix[1, ] == row[1] & plot_matrix[2, ] == row[2])]
  }))
  grid <- acast(grid, alpha ~ time, value.var = "bias")
  
  print(plot_ly(
    x = unique(plot_matrix[1, ]),
    y = log(unique(plot_matrix[2, ])),
    z = as.matrix(grid),
    type = "surface"
  ) %>%
    layout(
      title = "3D-Plot",
      scene = list(
        xaxis = list(title = "Time"),
        yaxis = list(title = list(text = "log(λ_D/λ_Y)")),
        zaxis = list(title = "Bias")
      )
    ))
  print(Stat.time-Sys.time()) 
  print(max)
}
plot_bias_curves_3D_comp.rate(a=-2, #U0 -> A
                      c=log(1), #A -> D
                      f=log(1), #A -> Y
                      g=2, #U0 -> D
                      h=log(1), #U0 -> Y
                      p=.2, #P(U_0=1)
                      comp.rate=25, #\lambdaD/lambdaY
                      haz=1/1000, #hazard of the event of interest and the competing event, if all parameters are 1
                      max=200)


################## Bias for different values of p #################
plot_bias_curves_3D_p <- function( a=0,
                                           c,
                                           f, 
                                           g,
                                           h,
                                           p=.2,
                                           comp.rate=25,
                                           haz, #baseline hazard for non-treated  
                                           max=50
){
  Stat.time<-Sys.time()
  CIF <- function(k, a, lambda_Y, comp.rate, zeta, gamma, theta, eta, alpha, p) {
    lambda_D=lambda_Y*comp.rate
    # First term calculations
    term1 <- (1 - exp((-lambda_Y * exp(zeta * a) - lambda_D * exp(gamma * a)) * k)) * 
      (lambda_Y * exp(zeta * a) / (lambda_Y * exp(zeta * a) + lambda_D * exp(gamma * a)))
    
    weight1 <- ((1 - p) / 2) / (p*( inv.logit(alpha))^a * ((1 - inv.logit(alpha)))^(1 - a) + (1 - p) / 2)
    
    # Second term calculations
    term2 <- (1 - exp((-lambda_Y * exp(zeta * a + theta) - lambda_D * exp(gamma * a + eta)) * k)) * 
      (lambda_Y * exp(zeta * a + theta) / (lambda_Y * exp(zeta * a + theta) + lambda_D * exp(gamma * a + eta)))
    
    weight2 <- (inv.logit(alpha)^a * (1 - inv.logit(alpha))^(1 - a) * p) / 
      ( p*( inv.logit(alpha))^a * ((1 - inv.logit(alpha)))^(1 - a) + (1 - p) / 2)
    
    # Final result
    result <- term1 * (weight1) + term2 * (weight2)
    return(result)
  }
  
  # All combinations of max and alpha 
  combinations <- expand.grid(max = 1:max, p = 1:99/100)
  
  num_columns <- nrow(combinations)
  
  # Matrix 
  plot_matrix <- matrix(0, nrow = 3, ncol = num_columns)
  

  plot_matrix[1, ] <- combinations$max
  plot_matrix[2, ] <- combinations$p
  

  
  for (i in 1:length(plot_matrix[1,])){
    y1 <- CIF(a=1,
              alpha = a, # überprüfen
              lambda_Y = haz,
              comp.rate = comp.rate,
              zeta = f,
              gamma = c,
              theta= h,
              eta = g,
              k= plot_matrix[1,i ],
              p= plot_matrix[2, i])
    y0 <- CIF(a=0,
              alpha = a, # überprüfen
              lambda_Y = haz,
              comp.rate = comp.rate,
              zeta = f,
              gamma = c,
              theta= h,
              eta = g,
              k= plot_matrix[1, i],
              p= plot_matrix[2, i])
    y1true <- CIF(a=1,
                  alpha = 0, # überprüfen
                  lambda_Y = haz,
                  comp.rate = comp.rate,
                  zeta = f,
                  gamma = c,
                  theta= h,
                  eta = g,
                  k= plot_matrix[1, i],
                  p= plot_matrix[2, i])
    y0true <- CIF(a=0,
                  alpha = 0, # überprüfen
                  lambda_Y = haz,
                  comp.rate = comp.rate,
                  zeta = f,
                  gamma = c,
                  theta= h,
                  eta = g,
                  k= plot_matrix[1, i],
                  p= plot_matrix[2, i])
    
    plot_matrix[3,i]<- (y1-y0)-(y1true-y0true)
    
  }

  grid <- expand.grid(time = unique(plot_matrix[1, ]), alpha = unique(plot_matrix[2, ]))
  grid$bias <- with(grid, apply(grid, 1, function(row) {
    plot_matrix[3, which(plot_matrix[1, ] == row[1] & plot_matrix[2, ] == row[2])]
  }))
  grid <- acast(grid, alpha ~ time, value.var = "bias")
  
  print(plot_ly(
    x = unique(plot_matrix[1, ]),
    y = unique(plot_matrix[2, ]),
    z = as.matrix(grid),
    type = "surface"
  ) %>%
    layout(
      title = "3D-Plot",
      scene = list(
        xaxis = list(title = "Time"),
        yaxis = list(title = list(text = "p")),
        zaxis = list(title = "Bias")
      )
    ))
  print(Stat.time-Sys.time()) 
  print(max)
}
plot_bias_curves_3D_p(a=-2, #U0 -> A
                              c=log(1), #A -> D
                              f=log(1), #A -> Y
                              g=2, #U0 -> D
                              h=log(1), #U0 -> Y
                              p=.5, #P(U_0=1)
                              comp.rate=25, #lambdaD/lambdaY
                              haz=1/1000, #hazard of the event of interest and the competing event, if all parameters are 1
                              max=200)
