library(shiny)
library(nortest)
library(stats)
library(VGAM)

generate_random_numbers_additive_congruence <- function(n, x0, x1, m) {
  m <- nchar(as.character(x0))
  result <- c(as.integer(x0), as.integer(x1))
  for (i in 3:n) {
    xi <- (result[i-1] + x0) %% 10^m
    result <- c(result, xi)
  }
  result <- result / 10^m
  return(result)
}

generate_random_numbers_linear_congruence <- function(n, seed, a, b, m) {
  m <- nchar(as.character(seed))
  result <- c(as.integer(seed))
  for (i in 2:n) {
    next_value <- (a * result[i-1] + b) %% 10^m
    result <- c(result, next_value)
  }
  result <- result / 10^m
  return(result)
}

generate_random_numbers_one_seed <- function(n, u) {
  result <- c(u)
  k <- nchar(as.character(u))
  for (i in 1:(n-1)) {
    u <- u^2
    u_str <- sprintf('%0*s', k*2, as.character(u))
    start <- ifelse(k %% 2 == 1, (k - 1) %/% 2, k %/% 2)
    u <- as.numeric(substr(u_str, start + 1, start + k))
    result <- c(result, u)
  }
  return(result / (10^k))
}

generate_random_numbers_two_seeds <- function(n, u1, u2) {
  result <- c(u1, u2)
  k <- nchar(as.character(u1))
  u <- u1
  for (i in 3:n) {
    u <- ifelse(i %% 2 == 1, u * u2, u * u1)
    u_str <- sprintf('%0*s', k*2, as.character(u))
    start <- ifelse(k %% 2 == 1, (k - 1) %/% 2, k %/% 2)
    u <- as.numeric(substr(u_str, start + 1, start + k))
    result <- c(result, u)
  }
  return(result / (10^k))
}

generate_random_numbers_seed_fix <- function(n, u, fix) {
  result <- c(u)
  k <- nchar(as.character(u))
  for (i in 2:n) {
    u <- u * fix
    u_str <- sprintf('%0*s', k*2, as.character(u))
    start <- ifelse(k %% 2 == 1, (k - 1) %/% 2, k %/% 2)
    u <- as.numeric(substr(u_str, start + 1, start + k))
    result <- c(result, u)
  }
  return(result / (10^k))
}

ranbern <- function(size, p, uniform_data) {
  y <- ifelse(uniform_data < p, 1, 0)
  hist(y, breaks = 0:2, main = "Bernoulli Distribution", xlab = "Outcome", col = "skyblue")
  return(y)
}

ranpoisson <- function(size, lambda, x) {
  y <- sapply(x, function(u) {
    n <- 0
    while (u > ppois(n, lambda)) {
      n <- n + 1
    }
    n
  })
  hist(y, main = "Poisson Distribution", xlab = "Value", col = "lightgreen")
  return(y)
}

rangeometric <- function(size, theta, x) {
  y <- sapply(x, function(u) {
    n <- 0
    while (u > pgeom(n, theta)) {
      n <- n + 1
    }
    n
  })
  hist(y, main = "Geometric Distribution", xlab = "Value", col = "salmon")
  return(y)
}

ranbino <- function(size, n, prob, x) {
  y <- sapply(x, function(u) {
    count <- 0
    while (u > pbinom(count, n, prob)) {
      count <- count + 1
    }
    count
  })
  hist(y, main = "Binomial Distribution", xlab = "Value", col = "purple")
  return(y)
}

ranunif <- function(size, a, b, x) {
  y <- floor(x * (b - a + 1)) + a
  hist(y, breaks = seq(a, b, 1), main = "Uniform Distribution", xlab = "Value", col = "gold")
  return(y)
}

rannegbinom <- function(size, r, p, x) {
  y <- sapply(x, function(u) {
    n <- 0
    while (u > pnbinom(n, r, p)) {
      n <- n + 1
    }
    n
  })
  hist(y, main = "Negative Binomial Distribution", xlab = "Value", col = "coral")
  return(y)
}


ranexp <- function(size, lambda, x) {
  y <- (-1 / lambda) * log(1 - x)
  hist(y, main = "Exponential Distribution", xlab = "Value", col = "coral")
  return(y)
}


rancauchy <- function(size, sigma, x) {
  y <- sigma * tan(pi * (x - 0.5))
  hist(y, main = "Cauchy Distribution", xlab = "Value", col = "lightpink")
  return(y)
}

rannorm <- function(size, mean, sd, x) {
  y <- qnorm(x) * sd + mean
  hist(y, main = "Normal Distribution", xlab = "Value", col = "lightblue")
  return(y)
}
ranhyper <- function(size, m, n, k,x) {
  y <- numeric(size)
  for (i in 1:size) {
    u <- x[i]
    j <- 0
    while (u > phyper(j, m, n, k)) {
      j <- j + 1
    }
    y[i] <- j
  }
  hist(y)
  return(y)
}
ranrayleigh <- function(size, sigma, x) {
  y <- sigma * sqrt(-log(1 - x))
  hist(y)
  return(y)
}
best_test <- function(data, dist_type, ...) {
  test_result <- list(p_value = NULL, test = NULL)
  
  params <- list(...)
  
  if (dist_type == "bernoulli_distribution") {
    test_result <- list(test = "Chi-Square", p_value = chisq.test(table(data))$p.value)
    
  } else if (dist_type == "poisson_distribution" && !is.null(params$lambda)) {
    test_result <- list(test = "Chi-Square", p_value = chisq.test(table(data))$p.value)
    
  } else if (dist_type == "geometric_distribution" && !is.null(params$theta)) {
    test_result <- list(test = "Chi-Square", p_value = chisq.test(table(data))$p.value)
    
  } else if (dist_type == "binomial_distribution" && !is.null(params$size) && !is.null(params$prob)) {
    test_result <- list(test = "Chi-Square", p_value = chisq.test(table(data))$p.value)
    
  } else if (dist_type == "uniform_distribution" && !is.null(params$a) && !is.null(params$b)) {
    test_result <- list(test = "Kolmogorov-Smirnov", p_value = ks.test(data, "punif", min = params$a, max = params$b)$p.value)
    
  } else if (dist_type == "negative_binomial_distribution" && !is.null(params$r) && !is.null(params$p)) {
    test_result <- list(test = "Chi-Square", p_value = chisq.test(table(data))$p.value)
    
  } else if (dist_type == "hypergeometric_distribution" && !is.null(params$m) && !is.null(params$n) && !is.null(params$k)) {
    expected_probs <- dhyper(0:max(data), params$m, params$n, params$k)
    test_result <- list(test = "Chi-Square", p_value = chisq.test(table(data), p = expected_probs, rescale.p = TRUE)$p.value)
    
  } else if (dist_type == "categorical_distribution" && !is.null(params$probs)) {
    test_result <- list(test = "Chi-Square", p_value = chisq.test(table(data), p = params$probs, rescale.p = TRUE)$p.value)
    
  } else if (dist_type == "exponential_distribution" && !is.null(params$lambda)) {
    test_result <- list(test = "Kolmogorov-Smirnov", p_value = ks.test(data, "pexp", params$lambda)$p.value)
    
  } else if (dist_type == "cauchy_distribution" && !is.null(params$sigma)) {
    test_result <- list(test = "Kolmogorov-Smirnov", p_value = ks.test(data, "pcauchy", location = 0, scale = params$sigma)$p.value)
    
  } else if (dist_type == "rayleigh_distribution" && !is.null(params$sigma)) {
    rayleigh_cdf <- function(x, sigma) 1 - exp(-x^2 / (2 * sigma^2))
    test_result <- list(test = "Kolmogorov-Smirnov", p_value = ks.test(data, rayleigh_cdf, params$sigma)$p.value)
    
  } else if (dist_type == "pareto_distribution" && !is.null(params$lambda) && !is.null(params$xm)) {
    pareto_cdf <- function(x) 1 - (params$xm / x)^params$lambda
    test_result <- list(test = "Kolmogorov-Smirnov", p_value = ks.test(data, pareto_cdf)$p.value)
    
  } else if (dist_type == "normal_distribution" && !is.null(params$mean) && !is.null(params$sd)) {
    test_result <- list(test = "Kolmogorov-Smirnov", p_value = ks.test(data, "pnorm", params$mean, params$sd)$p.value)
  }

  if (!is.null(test_result$p_value)) {
    if (test_result$p_value < 0.05) {
      return(paste("Based on p-value", test_result$p_value, "from the", test_result$test, "test, the data does not follow the", dist_type, "distribution."))
    } else {
      return(paste("Based on p-value", test_result$p_value, "from the", test_result$test, "test, the data follows the", dist_type, "distribution."))
    }
  } else {
    return("Error: Invalid distribution or missing parameters.")
  }
}


ui <- navbarPage(
  titlePanel("Random Number Simulation"),
  
  tabsetPanel(
    tabPanel("Uniform Number Generation",
             sidebarLayout(
               sidebarPanel(
                 selectInput("uniform_method", "Uniform Generation Method:",
                             choices = c("Use runif" = "runif",
                                         "Additive Congruence" = "additive_congruence",
                                         "Linear Congruence" = "linear_congruence",
                                         "Two Seeds" = "two_seeds",
                                         "One Seed" = "one_seed",
                                         "Seed and Fix" = "seed_fix"),
                             selected = "runif"),
                 numericInput("size", "Sample Size:", 100),
                 conditionalPanel(condition = "input.uniform_method == 'one_seed'",
                                  numericInput("seed", "Seed:", value = 123)),
                 conditionalPanel(condition = "input.uniform_method == 'seed_fix'",
                                  numericInput("fix", "Fix Value:", value = 2),
                                  numericInput("seed_fix_seed", "Seed:", value = 456)),
                 conditionalPanel(condition = "input.uniform_method == 'additive_congruence'",
                                  numericInput("x0", "Initial Value x0:", value = 123),
                                  numericInput("x1", "Initial Value x1:", value = 456)),
                 conditionalPanel(condition = "input.uniform_method == 'linear_congruence'",
                                  numericInput("linear_seed", "Seed:", value = 123),
                                  numericInput("slope", "Slope (a):", value = 456),
                                  numericInput("intercept", "Intercept (b):", value = 1013904223)),
                 conditionalPanel(condition = "input.uniform_method == 'two_seeds'",
                                  numericInput("seed1", "Seed 1:", value = 123),
                                  numericInput("seed2", "Seed 2:", value = 456)),
                 actionButton("generate_uniform", "Generate Uniform Numbers"),
                 actionButton("Analyze_uniform","Analyze numbers")
               ),
               mainPanel(
                 textOutput("uniform_output")
                 
               )
             )),
    
    tabPanel("Generate Numbers from Different Distributions",
             sidebarLayout(
               sidebarPanel(
                 selectInput("dist_choice", "Choose Distribution Method:",
                             choices = c(
                               "Pareto Distribution" = "pareto_distribution",
                               "Poisson Distribution" = "poisson_distribution",
                               "Geometric Distribution" = "geometric_distribution",
                               "Binomial Distribution" = "binomial_distribution",
                               "Bernoulli Distribution" = "bernoulli_distribution",
                               "Uniform Distribution" = "uniform_distribution",
                               "Negative Binomial Distribution" = "negative_binomial_distribution",
                               "Hypergeometric Distribution" = "hypergeometric_distribution",
                               "Categorical Distribution" = "categorical_distribution",
                               "Exponential Distribution" = "exponential_distribution",
                               "Cauchy Distribution" = "cauchy_distribution",
                               "Rayleigh Distribution" = "rayleigh_distribution",
                               "Normal Distribution" = "normal_distribution"
                             ),
                             selected = "exponential_distribution"),
                 
                 conditionalPanel(condition = "input.dist_choice == 'pareto_distribution'",
                                  numericInput("lambda_pareto", "Lambda:", 2),
                                  numericInput("xm_pareto", "xm:", 1)),
                 conditionalPanel(condition = "input.dist_choice == 'poisson_distribution'",
                                  numericInput("lambda", "Lambda:", 3)),
                 conditionalPanel(condition = "input.dist_choice == 'geometric_distribution'",
                                  numericInput("theta_geom", "Theta:", 0.5)),
                 conditionalPanel(condition = "input.dist_choice == 'binomial_distribution'",
                                  numericInput("theta_binom", "Trials (n):", 10),
                                  numericInput("prob_binom", "Probability (p):", 0.5)),
                 conditionalPanel(condition = "input.dist_choice == 'bernoulli_distribution'",
                                  numericInput("p", "Probability (p):", 0.5)),
                 conditionalPanel(condition = "input.dist_choice == 'uniform_distribution'",
                                  numericInput("a_unif", "Min (a):", 0),
                                  numericInput("b_unif", "Max (b):", 1)),
                 conditionalPanel(condition = "input.dist_choice == 'negative_binomial_distribution'",
                                  numericInput("r_neg_binom", "Failures (r):", 5),
                                  numericInput("p_neg_binom", "Probability (p):", 0.5)),
                 conditionalPanel(condition = "input.dist_choice == 'hypergeometric_distribution'",
                                  numericInput("m", "M:", 10),
                                  numericInput("n", "N:", 5),
                                  numericInput("k", "K:", 3)),
                 conditionalPanel(condition = "input.dist_choice == 'categorical_distribution'",
                                  textInput("probs_cat", "Probabilities (comma-separated):", "0.2,0.3,0.5")),
                 conditionalPanel(condition = "input.dist_choice == 'exponential_distribution'",
                                  numericInput("lambda_exp", "Lambda:", 1)),
                 conditionalPanel(condition = "input.dist_choice == 'cauchy_distribution'",
                                  numericInput("sigma_cauchy", "Sigma:", 1)),
                 conditionalPanel(condition = "input.dist_choice == 'rayleigh_distribution'",
                                  numericInput("sigma_rayleigh", "Sigma:", 1)),
                 conditionalPanel(condition = "input.dist_choice == 'normal_distribution'",
                                  numericInput("mean_norm", "Mean:", 0),
                                  numericInput("sd_norm", "Standard Deviation:", 1)),
                 
                 actionButton("generate_numbers", "Generate Numbers"),
                 actionButton("analyze_distribution", "Analyze")
               ),
               mainPanel(
                 plotOutput("dist_plot"),
                 uiOutput("dist_output"),
                 plotOutput("plots"),
                 textOutput("analysis_result"),
                 plotOutput("qq_plot")
               )
             ))
  )
)

server <- function(input, output, session) {

  uniform_data <- reactiveVal()
  distribution_data <- reactiveVal()
  
  observeEvent(input$generate_uniform, {
    size <- input$size
    method <- input$uniform_method
    result <- switch(method,
                     runif = runif(size, 0, 1),
                     additive_congruence = generate_random_numbers_additive_congruence(size, input$x0, input$x1, 7),
                     linear_congruence = generate_random_numbers_linear_congruence(size, input$linear_seed, input$slope, input$intercept, 7),
                     two_seeds = generate_random_numbers_two_seeds(size, input$seed1, input$seed2),
                     one_seed = generate_random_numbers_one_seed(size, input$seed),
                     seed_fix = generate_random_numbers_seed_fix(size, input$seed_fix_seed, input$fix))
    
    uniform_data(result)
    
    output$uniform_output <- renderText({
      paste("Generated Uniform Data:", paste(result, collapse = ", "))
    })
  })
  
  observeEvent(input$Analyze_uniform, {
    data <- uniform_data()
    if (is.null(data)) {
      showNotification("Please generate uniform data first.", type = "error")
      return()
    }
    ks_test <- ks.test(data, "punif")
    ad_test <- ad.test(data)
    chi_square <- sum((data - mean(data))^2 / mean(data))
    showNotification(
      paste(
        "Kolmogorov-Smirnov Test p-value:", round(ks_test$p.value, 4),
        "\nAnderson-Darling Test p-value:", round(ad_test$statistic, 4),
        "\nChi-Square Statistic:", round(chi_square, 4)
      ), 
      duration = 5, 
      type = "message"
    )
  })  
  observeEvent(input$generate_numbers, {
    data <- uniform_data()
    if (is.null(data) || length(data) == 0) {
      showNotification("Please generate uniform data first.", type = "error")
      return(NULL)
    }
    probs_cat <- as.numeric(unlist(strsplit(input$probs_cat, ",")))
    dist_function <- switch(input$dist_choice,
                            "bernoulli_distribution" = ranbern(input$size, input$p, data),
                            "poisson_distribution" = ranpoisson(input$size, input$lambda, data),
                            "geometric_distribution" = rangeometric(input$size, input$theta_geom, data),
                            "binomial_distribution" = ranbino(input$size, input$theta_binom, input$prob_binom, data),
                            "uniform_distribution" = ranunif(input$size, input$a_unif, input$b_unif, data),
                            "negative_binomial_distribution" = rannegbinom(input$size, input$r_neg_binom, input$p_neg_binom, data),
                            "exponential_distribution" = ranexp(input$size, input$lambda_exp, data),
                            "cauchy_distribution" = rancauchy(input$size, input$sigma_cauchy, data),
                            "normal_distribution" = rannorm(input$size, input$mean_norm, input$sd_norm, data),
                            "rayleigh_distribution" = ranrayleigh(input$size, input$sigma_rayleigh, data),
                            "categorical_distribution" = rancategorical(input$size, probs_cat, data),
                            "hypergeometric_distribution" = ranhyper(input$size, input$m, input$n, input$k, data),
                            "pareto_distribution" = ranpareto(input$size, input$lambda_pareto, input$xm_pareto, data))
    
    distribution_data(dist_function)
  })

  output$dist_plot <- renderPlot({
    if (!is.null(distribution_data())) {
      x <- distribution_data()
    
      par(mfrow = c(1, 2))  
      

      hist(x, col = "skyblue", main = "Histogram", xlab = "Values")
      

      boxplot(x, main = "Box Plot", col = "lightgreen", horizontal = TRUE)
    }
  })
  

  observeEvent(input$analyze_distribution, {
    data <- distribution_data()
    if (is.null(data) || length(data) == 0) {
      showNotification("Please generate distribution data first.", type = "error")
      return(NULL)
    }
    
 
    dist_type <- input$dist_choice
    analysis_result <- NULL
    

    tryCatch({
      analysis_result <- switch(dist_type,
                                "bernoulli_distribution" = best_test(data, dist_type, p = input$p),
                                "poisson_distribution" = best_test(data, dist_type, lambda = input$lambda),
                                "geometric_distribution" = best_test(data, dist_type, theta = input$theta_geom),
                                "binomial_distribution" = best_test(data, dist_type, size = input$theta_binom, prob = input$prob_binom),
                                "uniform_distribution" = best_test(data, dist_type, a = input$a_unif, b = input$b_unif),
                                "negative_binomial_distribution" = best_test(data, dist_type, r = input$r_neg_binom, p = input$p_neg_binom),
                                "hypergeometric_distribution" = best_test(data, dist_type, m = input$m, n = input$n, k = input$k),
                                "categorical_distribution" = {
                                  probs <- as.numeric(unlist(strsplit(input$probs_cat, ",")))
                                  best_test(data, dist_type, probs = probs)
                                },
                                "exponential_distribution" = best_test(data, dist_type, lambda = input$lambda_exp),
                                "cauchy_distribution" = best_test(data, dist_type, sigma = input$sigma_cauchy),
                                "rayleigh_distribution" = best_test(data, dist_type, sigma = input$sigma_rayleigh),
                                "pareto_distribution" = best_test(data, dist_type, lambda = input$lambda_pareto, xm = input$xm_pareto),
                                "normal_distribution" = best_test(data, dist_type, mean = input$mean_norm, sd = input$sd_norm),
   
                                {
                                  showNotification("Invalid distribution or missing parameters.", type = "error")
                                  NULL
                                })
      
      output$analysis_result <- renderText({
        if (is.null(analysis_result)) "Analysis could not be completed. Please check parameters." 
        else analysis_result
      })
      

      output$plots <- renderPlot({
        dist_type <- input$dist_choice
               x_vals <- seq(0, 10, length.out = 50)  
        x_vals_b <- 0:input$theta_binom
        probs_cat <- as.numeric(unlist(strsplit(input$probs_cat, ",")))
        type <- "l"
     
        theoretical_density <- switch(dist_type,
                                      "normal_distribution" = {
                                        x_vals_norm <- seq(input$mean_norm - 3 * input$sd_norm, input$mean_norm + 3 * input$sd_norm, length.out = 50)
                                        dnorm(x_vals_norm, mean = input$mean_norm, sd = input$sd_norm)
                                        },
                                      "poisson_distribution" = {
                                        type <- "h"
                                        x_vals <- 0:input$size
                                        dpois(x_vals, lambda = input$lambda)
                                        }
                                      ,
                                      "binomial_distribution" = {
                                        type <- "h"
                                        x_vals <- 0:input$theta_binom
                                        theoretical_density <- dbinom(x_vals, size = input$theta_binom, prob = input$prob_binom)
                                        },
                                      "geometric_distribution" = {
                                        x_vals <- 0:50
                                        dgeom(x_vals, input$theta_geom)}
                                      ,
                                      "exponential_distribution" = dexp(x_vals, rate = input$lambda_exp),
                                      "cauchy_distribution" = dcauchy(x_vals, scale = input$sigma_cauchy),
                                      "uniform_distribution" = {
                                        x_vals <- seq(input$a_unif,input$b_unif,lengh.out=)
                                        dunif(x_vals, min = input$a_unif, max = input$b_unif)
                                        },
                                      "negative_binomial_distribution" = {
                                        type <- "h"
                                        x_vals <- 0:50
                                        theoretical_density <- dnbinom(x_vals, size = input$r_neg_binom, prob = input$p_neg_binom) 
                                      },
                                      "bernoulli_distribution" = {
                                        type <- "h"
                                        x_vals <- c(0, 1)  
                                        theoretical_density <- dbinom(x_vals, size = 1, prob = input$p)
                                      },
                                      "hypergeometric_distribution" = {
                                          x_vals <- 0:min(input$k, input$m)
                                          theoretical_density <- dhyper(x_vals, m = input$m, n = input$n, k = input$k)
                                      },
                                      "categorical_distribution" = {
                                        num_categories <- length(probs_cat)    
                                        x_vals <- 1:num_categories               
                                        theoretical_density <- probs_cat      
                                      },
                                      "pareto_distribution" = dpareto(x_vals,shape = input$lambda_pareto,scale = input$xm_pareto),
                                      "rayleigh_distribution" = {
                                        sigma <- input$sigma_rayleigh
                                        if (sigma > 0) {
                                          d_rayleigh <- (x_vals / sigma^2) * exp(-x_vals^2 / (2 * sigma^2))  
                                          d_rayleigh[x_vals < 0] <- 0 
                                 
                                          d_rayleigh
                                        } else {
                                          rep(NA, length(x_vals)) 
                                        
                                        }
                                      },
                                      NULL)
        
   
        if (!is.null(theoretical_density)) {
          plot(x_vals, theoretical_density, type = type, col = "blue", lwd = 2,
               main = paste("Theoretical Density for", dist_type),
               xlab = "Values", ylab = "Density")
        }
        else {
          showNotification("Density plot is not available for this distribution type.", type = "error")
        }
      })
      

      output$qq_plot <- renderPlot({
        if (!is.null(data)) {
          theoretical <- switch(dist_type,
                                "normal_distribution" = qnorm(ppoints(data), mean = input$mean_norm, sd = input$sd_norm),
                                "poisson_distribution" = qpois(ppoints(data), lambda = input$lambda),
                                "binomial_distribution" = qbinom(ppoints(data), size = input$theta_binom, prob = input$prob_binom),
                                "geometric_distribution" = qgeom(ppoints(data), prob = input$theta_geom),
                                "exponential_distribution" = qexp(ppoints(data), rate = input$lambda_exp),
                                "cauchy_distribution" = qcauchy(ppoints(data), scale = input$sigma_cauchy),
                                "uniform_distribution" = qunif(ppoints(data), min = input$a_unif, max = input$b_unif),
                                "negative_binomial_distribution" = qnbinom(ppoints(data), size = input$r_neg_binom, prob = input$p_neg_binom),
                                "rayleigh_distribution" = sqrt(-2 * input$sigma_rayleigh^2 * log(1 - ppoints(data))), 
                                "pareto_distribution" = input$xm_pareto * (1 - ppoints(data))^(-1 / input$lambda_pareto),
                                "hypergeometric_distribution" = qhyper(ppoints(data), m = input$m, n = input$n, k = input$k),
                                "categorical_distribution" = as.numeric(as.factor(data)), 
                                "bernoulli_distribution" = qbinom(ppoints(data), size = 1, prob = input$p),  
                                NULL)
          
          if (!is.null(theoretical)) {
            qqplot(theoretical, data, main = paste("QQ Plot for", dist_type), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
            abline(0, 1, col = "red")
          } else {
            showNotification("QQ plot is not available for this distribution type.", type = "error")
          }
        } else {
          showNotification("Please generate data first.", type = "error")
        }
      })
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
}


shinyApp(ui = ui, server = server)

