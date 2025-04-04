###################################################
### Programm to derive the bias with competing confounding
### Author: Jost Viebrock
### 04/04/2025
############################
library(ggplot2)

# Constants
lambda_Y <- 1 / 1000
lambda_D <- 1 / 40
zeta <- log(.75)
gamma <- log(1)
p <- 0.2

# inv.logit function
inv_logit <- function(x) {
  return(1 / (1 + exp(-x)))
}

# Define the logistic function
inv_logit <- function(x) {
  return(1 / (1 + exp(-x)))
}

# Define the function for the bias
F_diff <- function(k, lambda_Y, lambda_D, zeta, gamma, theta, eta, alpha, p) {
  term1 <- (1 - exp(-lambda_Y * exp(zeta) * k - lambda_D * exp(gamma) * k)) * 
    (lambda_Y * exp(zeta)) / (lambda_Y * exp(zeta) + lambda_D * exp(gamma)) * 
    ((1 - p) / 2) / (p * inv_logit(alpha) + (1 - p) / 2)
  
  term2 <- (1 - exp(-lambda_Y * exp(zeta + theta) * k - lambda_D * exp(gamma + eta) * k)) * 
    (lambda_Y * exp(zeta + theta)) / (lambda_Y * exp(zeta + theta) + lambda_D * exp(gamma + eta)) * 
    (inv_logit(alpha) * p) / (p * inv_logit(alpha) + (1 - p) / 2)
  
  term3 <- (1 - exp(-lambda_Y * k - lambda_D * k)) * 
    (lambda_Y) / (lambda_Y + lambda_D) * 
    ((1 - p) / 2) / (p * (1 - inv_logit(alpha)) + (1 - p) / 2)
  
  term4 <- (1 - exp(-lambda_Y * exp(theta) * k - lambda_D * exp(eta) * k)) * 
    (lambda_Y * exp(theta)) / (lambda_Y * exp(theta) + lambda_D * exp(eta)) * 
    ((1 - inv_logit(alpha)) * p) / (p * (1 - inv_logit(alpha)) + (1 - p) / 2)
  
  result <- term1 + term2 - term3 - term4
  return(result)
}


##### The following plotting functions only differ by the header and the way how its pollted. They plot the same function.

### only for competing confounding 
generate_bias_plot_1 <- function(theta_vals, eta_vals, alpha_vals, filename) {
  library(ggplot2)
  
  # Define constants
  lambda_Y <- 1 / 1000
  lambda_D <- 1 / 40
  zeta <- log(0.75)
  gamma <- log(1)
  p <- 0.2
  
  # Generate data frame
  k_values <- seq(0, 200, length.out = 100)
  data <- expand.grid(k = k_values, index = 1:length(theta_vals))
  data$theta <- theta_vals[data$index]
  data$eta <- eta_vals[data$index]
  data$alpha <- alpha_vals[data$index]
  
  # Compute F_diff (assuming F_diff function is defined)
  data$F_diff <-  F_diff(k=data$k, lambda_Y, lambda_D, zeta, gamma, theta=data$theta, eta=data$eta, alpha=data$alpha, p) -
    F_diff(k=data$k, lambda_Y, lambda_D, zeta, gamma, theta=data$theta, eta=data$eta, 0, p)
  
  # Round values for legend
  data$theta_legend <- round(data$theta, 2)
  data$eta_legend <- round(data$eta, 2)
  
  # Erstelle ein neues Label für Facetten
  data$facet_label <- sprintf("α = %.0f, θ = %.0f", data$alpha, data$theta)
  
  # Plot
  plot <- ggplot(data, aes(x = k, y = F_diff, color = factor(eta_legend))) +
    geom_line(size = 2.5) +
    facet_wrap(~ facet_label, labeller = as_labeller(identity)) +
    theme_minimal() +
    labs(title = " ",
         x = "Time", 
         y = "Bias",
         linetype = expression(theta), 
         color = expression(eta)) +
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
    scale_y_continuous(labels = function(y) format(y, scientific = FALSE)) +
    theme(
      legend.position = "bottom",  
      axis.title = element_text(size = 20),   
      axis.text = element_text(size = 18),    
      legend.text = element_text(size = 16),  
      legend.title = element_text(size = 18),  
      strip.text = element_text(size = 22)  
    )
  
  # Save plot
  print(plot)
  #cairo_pdf(paste0("yourpath/", filename, ".pdf"), width = 8, height = 6)
  #print(plot)
  #dev.off()
}

theta_vals <-  c( log(1),  log(1),  log(1),log(1))
eta_vals <- c( 2,  -2,  5,-5)
alpha_vals <- c( -2,  -2,  -2,-2)

generate_bias_plot_1(theta_vals, eta_vals, alpha_vals, "Comp.Conf.only")


### for both confounding - different header

generate_bias_plot_2 <- function(theta_vals, eta_vals, alpha_vals, filename) {
  library(ggplot2)
  
  # Define constants
  lambda_Y <- 1 / 1000
  lambda_D <- 1 / 40
  zeta <- log(0.75)
  gamma <- log(1)
  p <- 0.2
  
  # Generate data frame
  k_values <- seq(0, 200, length.out = 100)
  data <- expand.grid(k = k_values, index = 1:length(theta_vals))
  data$theta <- theta_vals[data$index]
  data$eta <- eta_vals[data$index]
  data$alpha <- alpha_vals[data$index]
  
  # Compute F_diff (assuming F_diff function is defined)
  data$F_diff <-  F_diff(k=data$k, lambda_Y, lambda_D, zeta, gamma, theta=data$theta, eta=data$eta, alpha=data$alpha, p) -
    F_diff(k=data$k, lambda_Y, lambda_D, zeta, gamma, theta=data$theta, eta=data$eta, 0, p)
  
  # Combine theta and eta for legend (use text-based labels and expression for Greek letters)
  data$legend_label <- paste0(
    expression(θ), " = ", round(data$theta, 2), ", ", 
    expression(η), " = ", round(data$eta, 2)
  )
  
  # Plot
  plot <- ggplot(data, aes(x = k, y = F_diff, color = legend_label)) + 
    geom_line(size = 2.5) +  # All lines are solid
    facet_wrap(~ alpha, labeller = label_bquote(alpha == .(alpha))) +
    theme_minimal() +
    labs(title = " ",
         x = "Time", 
         y = "Bias",
         color = "") +
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
    scale_y_continuous(labels = function(y) format(y, scientific = FALSE)) +
    theme(
      legend.position = "bottom",  # Legende unter die Grafik verschieben
      axis.title = element_text(size = 20),   # Achsenbeschriftung größer
      axis.text = element_text(size = 18),    # Zahlen an den Achsen größer
      legend.text = element_text(size = 16),  # Legendenbeschreibung größer
      legend.title = element_text(size = 16),  # Legendentitel größer
      strip.text = element_text(size = 22) 
    )
  
  # Save plot
  print(plot)
  #cairo_pdf(paste0("yourpath/", filename, ".pdf"), width = 8, height = 6)
  #print(plot)
  #dev.off()
}
theta_vals <-  c(2,  -2, log(1.2))
eta_vals <- c(2,  -2,  log(1.5))
alpha_vals <- c(-2,  -2, -2)

generate_bias_plot_2(theta_vals, eta_vals, alpha_vals, "Both.Conf")
