library(shiny)
library(Seurat)
library(dplyr)

data <- readRDS('astrocytes.rds')

# Define UI for application that draws plots
ui <- navbarPage("Graphing Interface",
                 tabPanel(title = 'DimPlot',
                          sidebarPanel(
                              #type of group to plot
                              selectInput(inputId = 'Dim_group', label = 'Group by: ', selected = "ident",
                                          c('Cluster Group' = 'ident', 'Age' = 'youngvsold', 'Sample ID' = 'orig.ident')),
                              #options for type of reduction data to display
                              selectInput(inputId = 'Dim_reduction', label = 'Reduction Type:', selected = 'umap',
                                          c('UMAP' = 'umap','PCA' = 'pca')),
                              
                              radioButtons(inputId = 'Dim_split', label = 'Split by:', selected = 'none', inline = TRUE,
                                                 choices = c('None' = 'none','Age' = 'youngvsold', 'Sample' = 'orig.ident')),
                              #confirmation to plot
                              actionButton("dimplot", "Plot!")
                          ),
                          mainPanel(
                              plotOutput("DimPlot")
                          )
                 ),
                 tabPanel(title = 'VlnPlot',
                          sidebarPanel(
                              #gene text input
                              textInput(inputId = 'Vln_feature', label = 'Gene: (First letter must be capitalised)', 
                                        value = "", width = NULL, placeholder = "E.g. Aqp4"),
                              #options for type of data to display
                              selectInput(inputId = 'Vln_assay', label = 'Data Type:', selected = 'RNA',
                                          c('Raw' = 'RNA','Processed' = 'SCT')),
                              #confirmation to plot
                              actionButton("vlnplot", "Plot!")
                          ),
                          mainPanel(
                              plotOutput("VlnPlot")
                          )
                 ),
                 tabPanel(title = 'FeaturePlot',
                          sidebarPanel(
                              #gene text input
                              textInput(inputId = 'Feature_feature', label = 'Gene: (First letter must be capitalised)', 
                                        value = "", width = NULL, placeholder = "E.g. Aqp4"),
                              #confirmation to plot
                              actionButton("featureplot", "Plot!")
                          ),
                          mainPanel(
                              plotOutput("FeaturePlot")
                          )
                 ),
                 tabPanel(title = 'DotPlot',
                          sidebarPanel(
                              #gene text input
                              textInput(inputId = 'Dot_feature', label = 'Gene: (First letter must be capitalised)', 
                                        value = "", width = NULL, placeholder = "E.g. Aqp4"),
                              #options for type of data to display
                              selectInput(inputId = 'Dot_assay', label = 'Data Type:', selected = 'RNA',
                                          c('Raw' = 'RNA','Processed' = 'SCT')),
                              #confirmation to plot
                              actionButton("dotplot", "Plot!")
                          ),
                          mainPanel(
                              plotOutput("DotPlot")
                          )
                 ),
                 tabPanel(title = 'DoHeatmap',
                          sidebarPanel(
                              #gene text input
                              textInput(inputId = 'Heat_feature', label = 'Gene: (First letter must be capitalised)', 
                                        value = "", width = NULL, placeholder = "E.g. Aqp4"),
                              #confirmation to plot
                              actionButton("heatmap", "Plot!")
                          ),
                          mainPanel(
                              plotOutput("DoHeatmap")
                          )
                 )
                
)

# Define server logic required to draw plots
server <- function(input, output) {
    
    output$DimPlot <- renderPlot({
        #Start with empty plot
        if (input$dimplot == 0)
            return()
        
        #freeze the input till confirmation is given to plot
        input$dimplot
        group = isolate(input$Dim_group)
        reduction = isolate(input$Dim_reduction)
        split_by = isolate(input$Dim_split)
        #plot the data
        if (split_by == 'none') {
            DimPlot(data, label = T, pt.size = 0.4, group.by = group, split.by = NULL, 
                    reduction = reduction, ncol = 2)
        } else {
            DimPlot(data, label = T, pt.size = 0.4, group.by = group, split.by = split_by, 
                    reduction = reduction, ncol = 2)
        }
        
    })
    
    output$VlnPlot <- renderPlot({
        #Start with empty plot
        if (input$vlnplot == 0)
            return()
        
        #freeze the input till confirmation is given to plot
        input$vlnplot
        features = isolate(input$Vln_feature)
        assay = isolate(input$Vln_assay)
        
        #plot the data
        VlnPlot(data, features = strsplit(features, ' ')[[1]], assay = assay, ncol = 1)
        })
    
    output$FeaturePlot <- renderPlot({
        #Start with empty plot
        if (input$featureplot == 0)
            return()
        
        #freeze the input till confirmation is given to plot
        input$featureplot
        features = isolate(input$Feature_feature)
        
        #plot the data
        FeaturePlot(data, features = strsplit(features, ' ')[[1]], ncol = 1)
    })
    
    output$DotPlot <- renderPlot({
        #Start with empty plot
        if (input$dotplot == 0)
            return()
        
        #freeze the input till confirmation is given to plot
        input$dotplot
        features = isolate(input$Dot_feature)
        assay = isolate(input$Dot_assay)
        
        #plot the data
        DotPlot(data, features = strsplit(features, ' ')[[1]], assay = assay)
    })
    
    output$DoHeatmap <- renderPlot({
        #Start with empty plot
        if (input$heatmap == 0)
            return()
        
        #freeze the input till confirmation is given to plot
        input$heatmap
        features = isolate(input$Heat_feature)
        #plot the data
        DoHeatmap(data, features = strsplit(features, ' ')[[1]], assay = 'SCT')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
