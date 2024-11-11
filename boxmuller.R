library(shiny)
library(shinythemes)
library(ggplot2)

boxmuller <- function(n) {
  n <- ceiling(n / 2)
  u1 <- runif(n)
  u2 <- runif(n)
  z1 <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
  z2 <- sqrt(-2 * log(u1)) * sin(2 * pi * u2)
  return(c(z1, z2))
}

normal_random <- function(n, mu, sigma2) {
  standard_normals <- boxmuller(n)
  return(standard_normals * sqrt(sigma2) + mu)
}

chi_square_random <- function(n, k) {
  x <- numeric(n)
  for (i in 1:n) {
    u <- boxmuller(k)
    x[i] <- sum(u^2)
  }
  return(x)
}

t_student_random <- function(n, k) {
  x <- chi_square_random(n, k)
  z <- boxmuller(n)
  return(z / sqrt(x / k))
}

f_random <- function(n, k1, k2) {
  x1 <- chi_square_random(n, k1)
  x2 <- chi_square_random(n, k2)
  return((x1 / k1) / (x2 / k2))
}

log_normal_random <- function(n, mu, sigma2) {
  x <- normal_random(n, mu, sigma2)
  return(exp(x))
}

cauchy_random <- function(n, mu, k) {
  z1 <- boxmuller(n)
  z2 <- boxmuller(n)
  return(z1 / z2 * k + mu)
}

#------------------------------------------- UI -------------------------------------------------------
ui <- fluidPage(
  titlePanel("Random Number Generators with Box-Muller Transform Method"),
  
  theme = shinytheme("superhero"),
  
  tags$style(HTML("
       #generate,#compareToggle {
         background-color: #5981B1;
         color: white;
         border: none;
         padding: 10px 10px;
         text-align: center;
         font-size: 14px;
         cursor: pointer;
         border-radius: 15px; 
       }
       
 ")),
  sidebarLayout(
    sidebarPanel(
      selectInput("distribution", "Choose a Distribution:", 
                  choices = c("Normal", "Chi-Square", "t-Student", "F", "Log-Normal", "Cauchy")),
      numericInput("n", "Sample Size (n):", value = 100),
      
      # Conditional inputs based on distribution choice
      conditionalPanel(
        condition = "input.distribution == 'Normal'",
        numericInput("mu", "Mean (mu):", value = 0),
        numericInput("sigma2", "Variance (sigma^2):", value = 1)
      ),
      conditionalPanel(
        condition = "input.distribution == 'Chi-Square'",
        numericInput("k", "Degrees of Freedom (k):", value = 10)
      ),
      conditionalPanel(
        condition = "input.distribution == 't-Student'",
        numericInput("k", "Degrees of Freedom (k):", value = 10)
      ),
      conditionalPanel(
        condition = "input.distribution == 'F'",
        numericInput("k1", "Degrees of Freedom (k1):", value = 10),
        numericInput("k2", "Degrees of Freedom (k2):", value = 10)
      ),
      conditionalPanel(
        condition = "input.distribution == 'Log-Normal'",
        numericInput("mu", "Mean of Log (mu):", value = 0),
        numericInput("sigma2", "Variance of Log (sigma^2):", value = 1)
      ),
      conditionalPanel(
        condition = "input.distribution == 'Cauchy'",
        numericInput("mu", "Location (mu):", value = 0),
        numericInput("k", "Scale (k):", value = 1)
      ),
      
      actionButton("generate", "Generate"),
      actionButton("compareToggle", "Compare")
    ),
    
    mainPanel(
      
      conditionalPanel(
        condition = "input.generate > 0",
        plotOutput("histPlot"),
        h4("Generated Numbers"),
        tableOutput("dataTable"),
      ),
      
    )
  )
)

# ----------------------------------------------- server ----------------------------------------
server <- function(input, output, session) {
  compare_state <- reactiveVal(FALSE) # حالت پیش فرض: "Compare"
  
  observeEvent(input$compareToggle, {
    compare_state(!compare_state()) # تغییر وضعیت دکمه
    
    # به‌روزرسانی متن دکمه با توجه به وضعیت
    updateActionButton(session, "compareToggle",
                       label = if (compare_state()) "Don’t Compare" else "Compare")
  })
  
  data <-eventReactive(input$generate,{
    req(input$generate) # اجرا فقط وقتی که کاربر دکمه "Generate" را فشار دهد
    n <- input$n
    
    generated_data <- switch(input$distribution,
                             "Normal" = normal_random(n, input$mu, input$sigma2),
                             "Chi-Square" = chi_square_random(n, input$k),
                             "t-Student" = t_student_random(n, input$k),
                             "F" = f_random(n, input$k1, input$k2),
                             "Log-Normal" = log_normal_random(n, input$mu, input$sigma2),
                             "Cauchy" = cauchy_random(n, input$mu, input$k))
    
    rnorm_data <- switch(input$distribution,
                         "Normal" = rnorm(n, mean = input$mu, sd = sqrt(input$sigma2)),
                         "Chi-Square" = rchisq(n, df = input$k),
                         "t-Student" = rt(n, df = input$k),
                         "F" = rf(n, df1 = input$k1, df2 = input$k2),
                         "Log-Normal" = rlnorm(n, meanlog = input$mu, sdlog = sqrt(input$sigma2)),
                         "Cauchy" = rcauchy(n, location = input$mu, scale = input$k))
    
    combined_data <- data.frame(Value = c(generated_data, rnorm_data),
                                Method = rep(c("Box-Muller", "rnorm"), each = n))
    return(combined_data)
  })
  
  output$histPlot <- renderPlot({
    plot_data <- data()
    
    if (compare_state()) {
      ggplot(plot_data, aes(x = Value, fill = Method)) +
        geom_histogram(position = "identity", alpha = 0.4, bins = 30, color = "black") +
        scale_fill_manual(values = c("Box-Muller" = "#3498db", "rnorm" = "#e74c3c")) +
        labs(title = paste("Comparison of", input$distribution, "Generated and rnorm Distributions"),
             x = "Sample Value", y = "Frequency") +
        theme_minimal(base_size = 15) +
        theme(legend.position = "top") +
        geom_vline(aes(xintercept = mean(plot_data$Value[plot_data$Method == "Box-Muller"])), 
                   color = "#2980b9", linetype = "dashed", size = 1) +
        geom_vline(aes(xintercept = mean(plot_data$Value[plot_data$Method == "rnorm"])), 
                   color = "#c0392b", linetype = "dashed", size = 1)
    } else {
      ggplot(subset(plot_data, Method == "Box-Muller"), aes(x = Value)) +
        geom_histogram(fill = "steelblue", bins = 30, color = "black", alpha = 0.3) +
        labs(title = paste(input$distribution, "Distribution - Box-Muller Generated"),
             x = "Sample Value", y = "Frequency") +
        theme_classic()
    }
  })
  
  output$dataTable <- renderTable({
    generated_data <- data()$Value[data()$Method == "Box-Muller"]
    matrix_data <- matrix(generated_data, ncol = 17, byrow = TRUE)
    as.data.frame(matrix_data)
  }, colnames = FALSE, rownames = FALSE)
}

shinyApp(ui = ui, server = server)
