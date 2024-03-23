## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(colourpicker) # you might need to install this

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BF591 Assignment 7"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fileupload", "Load differential expression results", accept = ".csv"),
      radioButtons("x_variable", "Choose the column for X-axis",
                   choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                   selected = "log2FoldChange"),
      radioButtons("y_variable", "Choose the column for Y-axis",
                   choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                   selected = "padj"),
      colourInput("color1", "Base Point Color", value = "#FFCF56"),#22577A
      colourInput("color2", "Highlight Point Color", value = "#22577A"),
      sliderInput("magnitude_slider", "Select the magnitude of the p adjusted coloring:",
                  min = -300, max = -1, value = -150),
      actionButton("plot_button", "Plot")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("volcano")),
        tabPanel("Table", tableOutput("table"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #' load_Data
    #'
    #' @details Okay this one is a little weird but bear with me here. This is 
    #' still a "function", but it will take no arguments. The `reactive({})` bit 
    #' says "if any of my inputs (as in, input$...) are changed, run me again". 
    #' This is useful when a user clicks a new button or loads a new file. In 
    #' our case, look for the uploaded file's datapath argument and load it with 
    #' read.csv. Return this data frame in the normal return() style.
    load_data <- reactive({
      req(input$fileupload) 
      data <- read.csv(input$fileupload$datapath)
      return(data)
        return()
    })
    
    #' Volcano plot
    #'
    #' @param dataf The loaded data frame.
    #' @param x_name The column name to plot on the x-axis
    #' @param y_name The column name to plot on the y-axis
    #' @param slider A negative integer value representing the magnitude of
    #' p-adjusted values to color. Most of our data will be between -1 and -300.
    #' @param color1 One of the colors for the points.
    #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
    #'
    #' @return A ggplot object of a volcano plot
    #' @details I bet you're tired of these plots by now. Me too, don't worry.
    #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
    #' Write a normal volcano plot using geom_point, and integrate all the above 
    #' values into it as shown in the example app. The testing script will treat 
    #' this as a normal function.
    #' 
    #' !!sym() may be required to access column names in ggplot aes().
    #'
    #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
    volcano_plot <-
        function(dataf, x_name, y_name, slider, color1, color2) {
          plot <- ggplot(dataf, aes_string(x = x_name, y = sprintf("-log10(%s) + 1e-10", y_name), color = "factor(padj < 10^slider)")) +
            geom_point(size = 1, alpha = 0.7) +
            scale_color_manual(values = c(color2, color1), name = sprintf("padj < 1 × 10^%d", slider)) +
            labs(x = x_name, y = sprintf("-log10(%s)", y_name)) +
            theme_minimal() +
            theme(legend.position = "bottom", legend.direction = "horizontal")  
          
          return(plot)
        }
    
    #' Draw and filter table
    #'
    #' @param dataf Data frame loaded by load_data()
    #' @param slider Negative number, typically from the slider input.
    #'
    #' @return Data frame filtered to p-adjusted values that are less than 
    #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
    #' displayed.
    #' @details Same as above, this function is a standard R function. Tests will 
    #' evaluate it normally. Not only does this function filter the data frame to 
    #' rows that are above the slider magnitude, it should also change the format 
    #' of the p-value columns to display more digits. This is so that it looks 
    #' better when displayed on the web page. I would suggest the function 
    #' `formatC()`
    #'
    #' @examples draw_table(deseq_df, -210)
    #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
      filtered_data <- dataf %>%
        filter(padj < 10^slider) %>%  
        mutate(across(c(pvalue, padj), ~formatC(., format = "e", digits = 3)))  #Format p-value columns
      
      return(filtered_data)
    }
    
    #' These outputs aren't really functions, so they don't get a full skeleton, 
    #' but use the renderPlot() and renderTabel() functions to return() a plot 
    #' or table object, and those will be displayed in your application.
    output$volcano <- renderPlot({
      req(input$plot_button)
      data <- load_data()
      plot <- volcano_plot(data, input$x_variable, input$y_variable,
                           input$magnitude_slider, input$color1, input$color2)
      return(plot)
    }) # replace this NULL

    
    # Same here, just return the table as you want to see it in the web page
    output$table <- renderTable({
      data <- load_data()
      filtered_table <- draw_table(data, input$magnitude_slider)
      return(filtered_table)
    }) # replace this NULL
}

# Run the application
shinyApp(ui = ui, server = server)
