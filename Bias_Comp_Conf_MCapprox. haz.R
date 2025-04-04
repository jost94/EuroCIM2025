################################################
### Monte Carlo approximation of the bias
### Author: Jost Viebrock
### 04/04/2025
##############################################

### in.logit function 
inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

### Generation of the Data
gen_data <- function(a=0, #U0 -> A 
                     c, #A -> D 
                     f, #A -> Y
                     g, #U0 -> D
                     h, # U0 -> Y
                     p=1/2, # prevalence of comp.ev
                     k_comp_ev=1, #shape parameter of weibull distribution for the competing event
                     k_outcome=1, #shape parameter of weibull distribution for the competing event
                     cens.rate = 1/4, #regulates propotrion of censored data
                     comp.rate = 1, #regulates propotrion of competed data
                     haz, #baseline hazard for non-treated  
                     n = 10^5, #number of samples
                     max= 500){ #expected time, where the time dependent confounder does or does not occur
  Start.Time <- Sys.time() 
  ##set all parameters, which are needed for the simulation
  #one dimensional parameters
  
  cf <- 1/cens.rate-1 # censoring factor for exponential distribution
  U0 <- rbinom(n,1,p) #runif(n) #baseline confounder
  treat.prob <- inv.logit(a*U0)*rep(1,n) #inv.logit(U0*a*rep(1,n)) #probability of treatment; 1/2 if a=0
  #later (un)observed variables
  comp.ev <- rep(0,n) #competing event {0,1}
  outcome <- rep(0,n) #outcome {0,1}
  censoring <- rep(0,n) #censoring {0,1}
  event.type <- rep(0,n) #event type {0,1,2}
  Y <- rep(0,n) #observed time of outcome
  #baseline variables
  
  treated <- rbinom(n,1,treat.prob) #treatment of any person {0,1}
  #simulating the variables
  for (i in 1:n){
    comp.ev[i] <- rweibull(n=1, scale =1/(haz*exp(c*treated[i])*exp(g*U0[i])*comp.rate), shape = k_comp_ev)
    outcome[i] <- rweibull(n=1, scale =1/(haz*exp(f*treated[i])*exp(h*U0[i])), shape =k_outcome)
    Y[i] <- min(comp.ev[i],outcome[i])
    #keeping only the one, which occurs first
    if(comp.ev[i] < outcome[i]){
      outcome[i] <- 0
      comp.ev[i] <- 1
      event.type[i] <- 2
    }
    else if(comp.ev[i] > outcome[i]){
      outcome[i] <- 1
      comp.ev[i] <- 0
      event.type[i] <- 1
    }
    censoring[i] <- 10^5#rweibull(n=1, scale = 1/ (haz*exp(f*treated[i])*exp(h*U0[i])*cens.rate), shape =k_outcome )
    if(censoring[i] < Y[i]){
      outcome[i] <- 0
      comp.ev[i] <- 0
      Y[i] <- censoring[i]
      censoring[i] <- 1 
      event.type[i] <- 0
    }
    else{
      censoring[i] <- 0
    }
    
  }
  #print(sum(censoring)/length(censoring)) #checking if the censoring rate is correct. 
  #Y <- ceiling(Y)
  DATA_1 <- cbind(Y,outcome,treated,comp.ev,U0,rep(1,n),censoring,event.type) #combining all information in a matrix
  DATA <- as.data.frame(DATA_1)
  DATA[DATA$Y> (max+3),]$outcome <- 0
  DATA[DATA$Y> (max+3),]$comp.ev <- 0
  DATA[DATA$Y> (max+3),]$censoring <- 1
  DATA[DATA$Y> (max+3),]$event.type <- 0
  DATA[DATA$Y> (max+3),]$Y <- (max+3)
  DATA
}
### Estimation

model_simu <- function(a=0, 
                       c,
                       f, 
                       g,
                       h,
                       p=1/2,
                       k_comp_ev=1, #shape parameter of weibull distribution for comp_ev
                       k_outcome=1, #shape parameter of weibull distribution for outcome
                       cens.rate = 1/4, #rough proportion of censored data
                       comp.rate=1,
                       haz, #baseline hazard for non-treated  
                       n = 10^5, #number of individuals
                       max= 500,
                       ymax=1){
  start.time <- Sys.time()
  DATA_1 <- gen_data(a=a,c=c,f=f,g=g,h=h,p=p,k_comp_ev=k_comp_ev,k_outcome=k_outcome,cens.rate = cens.rate,comp.rate=comp.rate,haz = haz,n=n,max=max) # generate data from previous functions
  model_tot_comp <- cuminc(ftime = DATA_1[,1],
                           fstatus = DATA_1[,8],
                           group = DATA_1[,3])
  return(model_tot_comp)
}

total_eff <- function(model, max_points = 10000){ # Funktion to calculate the difference
  curve_1_est <- model$`0 1`$est
  curve_1_time <- model$`0 1`$time
  curve_1 <- data.frame(time = curve_1_time, estimate = curve_1_est)
  
  curve_2_est <- model$`1 1`$est  
  curve_2_time <- model$`1 1`$time
  curve_2 <- data.frame(time = curve_2_time, estimate = curve_2_est)
  
  index_nearest <- findInterval(curve_1$time, curve_2$time)
  
  # Extrahieren der Schätzwerte aus curve_2 für die entsprechenden Zeitpunkte
  nearest_estimates_curve_2 <- curve_2$estimate[index_nearest]
  
  # Berechnen der Differenz zwischen curve_1 und den zugewiesenen Werten von curve_2
  curve_diff <- data.frame(
    time = curve_1$time,
    estimate_diff = nearest_estimates_curve_2 - curve_1$estimate
  )
  # Reduzieren auf maximal max_points
  if (nrow(curve_diff) > max_points) {
    indices <- seq(1, nrow(curve_diff), length.out = max_points) # Gleichmäßige Auswahl
    curve_diff <- curve_diff[round(indices), ]
  }
  return(curve_diff)
}


## for plotting curves 
plot_curves_k_comp_ev <- function(a=0, 
                                  c,
                                  f, 
                                  g,
                                  h,
                                  p=1/2,
                                  k_comp_ev=1, #shape parameter of weibull distribution for comp_ev
                                  k_comp_ev1=2,
                                  k_comp_ev2=0.5,
                                  k_outcome=1, #shape parameter of weibull distribution for outcome
                                  cens.rate = 1/4, #rough proportion of censored data
                                  comp.rate=1,
                                  haz, #baseline hazard for non-treated  
                                  n = 10^7, #number of individuals
                                  max=500,
                                  ymax=1){
  Start.Time <- Sys.time()
  
  # Simulating  with different values of k
  model_biased <- model_simu(a=a, c=c, f=f, g=g, h=h, p=p, 
                             k_comp_ev=k_comp_ev, 
                             k_outcome=k_outcome, 
                             cens.rate=cens.rate, comp.rate=comp.rate, 
                             haz=haz, n=n, max=max, ymax=ymax)
  
  model_biased_1 <- model_simu(a=a, c=c, f=f, g=g, h=h, p=p, 
                               k_comp_ev=k_comp_ev1, 
                               k_outcome=k_outcome, 
                               cens.rate=cens.rate, comp.rate=comp.rate/2/gamma(1+1/k_comp_ev1), # to make the expactation values comparable
                               haz=haz, n=n, max=max, ymax=ymax)
  
  model_biased_2 <- model_simu(a=a, c=c, f=f, g=g, h=h, p=p, 
                               k_comp_ev=k_comp_ev2, 
                               k_outcome=k_outcome, 
                               cens.rate=cens.rate, comp.rate=comp.rate/2/gamma(1+1/k_comp_ev2), # to make the expactation values comparable
                               haz=haz, n=n, max=max, ymax=ymax)
  
  model_unbiased <- model_simu(a=0, c=c, f=f, g=g, h=h, p=p, 
                               k_comp_ev=k_comp_ev, 
                               k_outcome=k_outcome, 
                               cens.rate=cens.rate, comp.rate=comp.rate, 
                               haz=haz, n=n, max=max, ymax=ymax)
  model_unbiased_1 <- model_simu(a=0, c=c, f=f, g=g, h=h, p=p, 
                               k_comp_ev=k_comp_ev1, 
                               k_outcome=k_outcome, 
                               cens.rate=cens.rate, comp.rate=comp.rate/2/gamma(1+1/k_comp_ev1), # to make the expactation values comparable
                               haz=haz, n=n, max=max, ymax=ymax)
  model_unbiased_2 <- model_simu(a=0, c=c, f=f, g=g, h=h, p=p, 
                               k_comp_ev=k_comp_ev2, 
                               k_outcome=k_outcome, 
                               cens.rate=cens.rate, comp.rate=comp.rate/2/gamma(1+1/k_comp_ev2), # to make the expactation values comparable
                               haz=haz, n=n, max=max, ymax=ymax)
  
  
  interpol<- function(Name, max=200) { #Interpolation for calculation of the difference of the CIFs

    new_time <- seq(0, max, length.out = max * 5)
    
    interpolated_values <- approx(x = Name$time, y = Name$estimate_diff, xout = new_time, rule = 2)
    
    new_data <- data.frame(time = interpolated_values$x, estimate_diff = interpolated_values$y)
    
    return(new_data)
  }
  
  
  total_eff_biased <- interpol(total_eff(model_biased)) - interpol(total_eff(model_unbiased))
  total_eff_biased_1 <- interpol(total_eff(model_biased_1)) - interpol(total_eff(model_unbiased_1))
  total_eff_biased_2 <- interpol(total_eff(model_biased_2)) - interpol(total_eff(model_unbiased_2))
  total_eff_unbiased <- total_eff(model_unbiased) 
  total_eff_biased$time <- total_eff_biased_1$time <- total_eff_biased_2$time <- interpol(total_eff(model_biased))$time
  print(total_eff_biased)
  # Maximum  y-Axes
  ymax <- max(abs(total_eff_biased$estimate_diff), 
              abs(total_eff_biased_1$estimate_diff), 
              abs(total_eff_biased_2$estimate_diff))


  # Plot 
  plot1 <- plot(total_eff_biased$time, total_eff_biased$estimate_diff, type="l", col="blue", lty=1,
       ylim=c(-ymax, ymax), lwd=4,
       xlab="Time", ylab="Bias" 
       , cex.lab=1.5, cex.axis=1.5, cex.main=2)  
  

  lines(total_eff_biased_1$time, total_eff_biased_1$estimate_diff, type="l", col="green", lwd=4, lty=1)
  lines(total_eff_biased_2$time, total_eff_biased_2$estimate_diff, type="l", col="red", lwd=4, lty=1)
  abline(h=0)
  
  
  legend("topleft",  # Verschiebt die Legende von oben links nach unten rechts
         legend = c(
           substitute(paste( s[D], " = ", k_comp_ev), list(k_comp_ev = k_comp_ev)),
           substitute(paste( s[D], " = ", k_comp_ev1), list(k_comp_ev1 = k_comp_ev1)),
           substitute(paste( s[D], " = ", k_comp_ev2), list(k_comp_ev2 = k_comp_ev2))
         ), 
         col=c(#"blue",
           "blue", "green", "red"), 
         lty=c(#3,
           1, 1, 1), 
         lwd=5, 
         bty="n", cex=1.5) 
  
  legend("bottomleft", bty="n", cex=1.5,
         legend=c(
           #expression(paste("Parameters:")),
           substitute(paste(alpha, " = ", a), list(a=round(a, 3))),
           substitute(paste(gamma, " = ", c), list(c=round(c, 3))),
           substitute(paste(zeta, " = ", f), list(f=round(f, 3))),
           substitute(paste(eta, " = ", g), list(g=round(g, 3))),
           substitute(paste(theta, " = ", h), list(h=round(h, 3))),
           substitute(paste("p = ", p), list(p=round(p, 3)))
         ))
  print(plot1)
  print(Sys.time() - Start.Time)
}




#Simulation

#output_file <- "yourpath/plot_output_4_s_D_new.pdf"


#pdf(output_file, width=8,height=6)
plot_curves_k_comp_ev(a=-2, #U0 -> A
                      c=log(1), #A -> D
                      f=log(.75), #A -> Y
                      g=2, #U0 -> D
                      h=log(1), #U0 -> Y
                      p=.2, #P(U_0=1)
                      k_comp_ev= 1, #shape of weibull distribution for the competing event 
                      k_outcome= 1, # shape of weibull distribution for the outcome
                      cens.rate = 1/10, # rough rate of how much individuals are independently censored
                      comp.rate=25,
                      haz=1/1000, #hazard of the event of interest and the competing event, if all parameters are 1
                      n = 10^5, #sample size # time when we expect the time-dependent confounder to occur
                      max = 200,
                      ymax=.01)
#dev.off()
