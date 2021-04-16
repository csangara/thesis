library(shiny)

ui <- fluidPage(
  selectInput("plot_type", "Select plot type",
              list("Performance (average)" = "perf_avg",
                   "Performance (distribution)" = "perf_dist",
                   "UMAP of predictions" = "umap_pred" )),
  uiOutput("plot_type_args"),
  textOutput("result")
)

server <- function(input, output) {
  output$plot_type_args <- renderUI({
    if (input$plot_type == "perf_avg"){
      tagList(
      selectInput("x_axis", "x-axis",
                  list("Method" = "method", "Dataset type" = "dataset_type",
                       "Dataset" = "dataset")),
      selectInput("y_axis", "metric",
                  list("Correlation" = "corr", "RMSE" = "RMSE",
                       "Accuracy" = "accuracy", "Sensitivity" = "sensitivity",
                       "Specificity" = "specificity", "Precision" = "precision",
                       "F1 score" = "F1")),
      actionButton("plot_now", "Plot")
      )
    }
  })
  
  data <- eventReactive(input$plot_now, {
    paste("You chose", input$x_axis, "and", input$y_axis)
  })
  
  output$result <- renderText({
    data()
  })

}

shinyApp(ui = ui, server = server)