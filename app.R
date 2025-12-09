#Carregar as bibliotecas
library(shiny)
library(shinydashboard) 
library(readr) 
library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(microbenchmark)
library(rsconnect)


# 1. Distância Euclidiana (Loops)
euclidean_distance <- function(ind_i, ind_j) {
  ind_i <- as.numeric(ind_i); ind_j <- as.numeric(ind_j)
  soma <- 0
  for (k in seq_along(ind_i)) {
    if (!is.na(ind_i[k]) && !is.na(ind_j[k])) {
      soma <- soma + (ind_i[k] - ind_j[k])^2
    }
  }
  return(sqrt(soma))
}

# 2. Correlação de Pearson (Loops)
pearson_correlation <- function(ind_i, ind_j) {
  ind_i <- as.numeric(ind_i); ind_j <- as.numeric(ind_j)
  mean_i <- mean(ind_i, na.rm = TRUE); mean_j <- mean(ind_j, na.rm = TRUE)
  num <- 0; den_i <- 0; den_j <- 0
  for (k in seq_along(ind_i)) {
    if (!is.na(ind_i[k]) && !is.na(ind_j[k])) {
      num  <- num  + (ind_i[k] - mean_i) * (ind_j[k] - mean_j)
      den_i <- den_i + (ind_i[k] - mean_i)^2
      den_j <- den_j + (ind_j[k] - mean_j)^2
    }
  }
  corr <- num / sqrt(den_i * den_j)
  if (is.nan(corr) || (den_i * den_j) == 0) { return(0) } else { return(corr) }
}

# 3. Média de Diferenças Absolutas (Loops)
mean_abs_diff <- function(ind_i, ind_j) {
  ind_i <- as.numeric(ind_i); ind_j <- as.numeric(ind_j)
  soma <- 0; count <- 0
  for (k in seq_along(ind_i)) {
    if (!is.na(ind_i[k]) && !is.na(ind_j[k])) {
      soma <- soma + abs(ind_i[k] - ind_j[k])
      count <- count + 1
    }
  }
  if (count == 0) { return(NA) } else { return(soma / count) }
}

# 4. Contagem Acima do Threshold (Loops)
count_above_threshold <- function(ind_i, ind_j, threshold = 10) {
  ind_i <- as.numeric(ind_i); ind_j <- as.numeric(ind_j)
  count <- 0
  for (k in seq_along(ind_i)) {
    if (!is.na(ind_i[k]) && !is.na(ind_j[k])) {
      if (abs(ind_i[k] - ind_j[k]) > threshold) {
        count <- count + 1
      }
    }
  }
  return(count)
}

# 5. Índice Composto 
composite_index <- function(ind_i, ind_j, threshold = 10) {
  d_euclid <- euclidean_distance(ind_i, ind_j)
  corr      <- pearson_correlation(ind_i, ind_j)
  mad      <- mean_abs_diff(ind_i, ind_j)
  cat_val  <- count_above_threshold(ind_i, ind_j, threshold)
  
  ci <- (d_euclid + abs(corr) + mad) / (cat_val + 1)
  return(ci)
}


# 6. Função para Construir Matriz
build_matrix <- function(data, index_fun, threshold = 10) {
  n <- nrow(data)
  M <- matrix(0, n, n)
  for (i in 1:n) {
    ind_i <- as.numeric(data[i, ])
    for (j in 1:n) {
      ind_j <- as.numeric(data[j, ])
      M[i, j] <- index_fun(ind_i, ind_j, threshold = threshold)
    }
  }
  return(M)
}


# --- Versão Vetorizada ---

mad_matrix <- function(X) {
  n <- nrow(X)
  M <- matrix(0, n, n)
  for (i in 1:n) {
    diffs <- abs(sweep(X, 2, X[i, ], "-"))
    M[i, ] <- rowMeans(diffs, na.rm = TRUE)
  }
  return(M)
}

cat_matrix <- function(X, threshold = 10) {
  n <- nrow(X)
  M <- matrix(0, n, n)
  for (i in 1:n) {
    diffs <- abs(sweep(X, 2, X[i, ], "-"))
    M[i, ] <- rowSums(diffs > threshold, na.rm = TRUE)
  }
  return(M)
}

composite_matrix_vec <- function(X, threshold = 10) {
  mat_euclid <- as.matrix(dist(X, method = "euclidean"))
  mat_corr   <- cor(t(X), use = "pairwise.complete.obs")
  mat_mad    <- mad_matrix(X)
  mat_cat    <- cat_matrix(X, threshold)
  S_vec <- (mat_euclid + abs(mat_corr) + mat_mad) / (mat_cat + 1)
  return(S_vec)
}


# ---------------- UI ----------------

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = "Otimização de Tempo Computacional em R",
    titleWidth = 400
  ),
  
  dashboardSidebar(
    width = 300,
    
    fileInput("file_upload", "Carregar Ficheiro de Biomarcadores (.csv ou .xlsx):",
              accept = c(".csv", ".txt", ".xlsx")),
    
    tags$hr(),
    
    selectInput("metric_select", "Selecionar Métrica de Similaridade:",
                choices = c("Índice Composto" = "Composite", 
                            "Distância Euclidiana" = "Euclidean", 
                            "Correlação de Pearson" = "Pearson", 
                            "Diferenças Absolutas Médias" = "MAD")),
    
    tags$hr(),
    
    sliderInput("threshold_val", "Definir Threshold para Contagem (usado no Composite):",
                min = 0.01, max = 50, value = 10.0, step = 0.1),
    
    tags$hr(),
    
    actionButton("run_analysis", "Executar Análise e Comparação", icon = icon("ruler-combined"))
  ),
  
  dashboardBody(
    fluidRow(
      box(
        title = "Heatmap da Matriz de Similaridade (Versão Vectorizada)", status = "primary", solidHeader = TRUE,
        width = 6, height = "500px",
        plotOutput("similarity_heatmap", height = "450px")
      ),
      
      box(
        title = "Comparação de Tempo de Execução do Índice Composto", status = "warning", solidHeader = TRUE,
        width = 6, height = "500px",
        plotOutput("timing_comparison", height = "200px"),
        tags$hr(),
        h4("Tempos Médios Reais (Microbenchmark):"),
        tableOutput("timing_table")
      )
    )
  )
)


# ---------------- SERVER ----------------

server <- function(input, output, session) {
  
  data_input <- reactive({
    req(input$file_upload)
    path <- input$file_upload$datapath
    
    if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
      df <- readxl::read_excel(path) %>% as.data.frame()
    } else {
      df <- read.csv2(path)
    }
    
    dados_matrix <- df %>%
      mutate(across(everything(), ~as.numeric(as.character(.)))) %>%
      as.matrix()
    
    if (ncol(dados_matrix) > 1 && all(is.na(dados_matrix[, 1]))) {
      dados_matrix <- dados_matrix[, -1]
    }
    return(dados_matrix)
  })
  
  
  analysis_results <- eventReactive(input$run_analysis, {
    dados <- data_input()
    threshold <- input$threshold_val
    
    times <- microbenchmark(
      Loops_Base = build_matrix(dados, composite_index, threshold = threshold),
      Vectorized_Matrix = composite_matrix_vec(dados, threshold = threshold),
      Modular_Loops = build_matrix(dados, composite_index, threshold = threshold),
      times = 10
    )
    
    timing_summary <- summary(times) %>%
      select(expr, mean, median, min, max) %>%
      mutate(expr = as.character(expr), mean = mean / 1e9, median = median / 1e9) %>%
      rename(Versão = expr, `Tempo Médio (s)` = mean, `Tempo Mediana (s)` = median)
    
    S_matrix <- composite_matrix_vec(dados, threshold = threshold)
    
    mat_euclid <- as.matrix(dist(dados, method = "euclidean"))
    mat_corr   <- cor(t(dados), use = "pairwise.complete.obs")
    mat_mad    <- mad_matrix(dados)
    
    return(list(
      similarity_matrix = S_matrix,
      timing_results = timing_summary,
      mat_euclid = mat_euclid,
      mat_corr = mat_corr,
      mat_mad = mat_mad
    ))
  }, ignoreNULL = FALSE)
  
  
  # Heatmap
  output$similarity_heatmap <- renderPlot({
    results <- analysis_results()
    
    if (input$metric_select == "Composite") { M <- results$similarity_matrix
    } else if (input$metric_select == "Euclidean") { M <- results$mat_euclid
    } else if (input$metric_select == "Pearson") { M <- results$mat_corr
    } else if (input$metric_select == "MAD") { M <- results$mat_mad }
    
    matriz_df <- M %>%
      as.data.frame() %>%
      mutate(Var1 = 1:nrow(M)) %>%
      tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "Valor") %>%
      mutate(Var2 = as.integer(gsub("V", "", Var2)))
    
    ggplot(matriz_df, aes(x = Var1, y = Var2, fill = Valor)) +
      geom_tile() +
      scale_fill_gradientn(colors = c("white", "skyblue", "darkblue")) +
      labs(title = paste("Heatmap da Matriz de", input$metric_select),
           x = "Indivíduo i", y = "Indivíduo j", fill = input$metric_select) +
      theme_minimal()
  })
  
  
  # Gráfico do Tempo
  output$timing_comparison <- renderPlot({
    results <- analysis_results()
    req(results$timing_results)
    timing <- results$timing_results
    
    ggplot(timing, aes(x = Versão, y = `Tempo Médio (s)`, fill = Versão)) +
      geom_bar(stat = "identity") +
      labs(title = "Comparação de Tempos de Execução (Tempo Médio)",
           y = "Tempo (segundos)", x = "") +
      theme_classic()
  })
  
  # Tabela
  output$timing_table <- renderTable({
    results <- analysis_results()
    req(results$timing_results)
    results$timing_results
  }, digits = 5, striped = TRUE, bordered = TRUE)
  
}

shinyApp(ui, server)




