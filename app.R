library(shiny)
library(shinysky)
library(glue)
library(ggplot2)
library(DESeq2)
library(shinyWidgets)

df <- readRDS("rld_tidy.rds")
all_gene_symbols <- unique(df$gene)

num_col = 5
group_labels <- c("Female 20 DPH", "Female 50 DPH", "Male 20 DPH", "Male 50 DPH")
paired_palette <- c("#E69F00", "#D55E00", "#009E73", "#0072B2")


ui <- fluidPage(
  titlePanel("Developmental RA Gene Expression"),
  ## Autocomplete format:
  textInput.typeahead(id = "gene_symbol",
                      placeholder = "Enter gene symbol",
                      local = data.frame(name = c(all_gene_symbols)),
                      valueKey = "name",
                      tokens = c(1:length(all_gene_symbols)),
                      template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-name'>{{name}}</p>")
                      ),
  actionBttn(
    inputId = "add",
    label = "Add",
    style = "simple",
    color = "primary",
    size = "sm"
  ),
  actionBttn(
    inputId = "clear",
    label = "Clear plot",
    style = "bordered",
    color = "primary",
    size = "sm"
  ),
  downloadBttn(
    outputId = "download",
    label = "Download plot",
    style = "bordered",
    color = "primary",
    size = "sm"
  ),
  br(),
  br(),
  plotOutput(outputId = "abundance"),
  
)

server <- function(input, output, session) {
  # create reactive values object
  myValues <- reactiveValues()
  
  # when "Add" button is pressed, append to genelist
  # condition checks that gene is not already in list
  # calculate width for plot based on number of genes
  observeEvent(input$add, {
    if (!(isolate(input$gene_symbol) %in% isolate(myValues$genelist))) {
      myValues$genelist <-
        c(isolate(myValues$genelist),
          isolate(input$gene_symbol))
    }
    if (length(isolate(myValues$genelist <= num_col))) {
      myValues$width = (2 * length(isolate(myValues$genelist)))
    } else {
      myValues$width = 2 * num_col
    }
  })
  
  # build plot whenever "Add" button is clicked
  plot_reactive <- eventReactive(input$add, {
    data <- df[df$gene %in% isolate(myValues$genelist), ]
    p <- ggplot(data, aes(x=group, y=log2)) +
      geom_point(aes(color=group), position=position_jitter(w = 0.1,h = 0)) +  
      labs(color="") +
      facet_wrap(~ gene, ncol=num_col, scales="free_y") + 
      xlab("") +
      ylab("log2 abundance") +
      scale_x_discrete(labels=group_labels) +
      scale_color_manual(values = paired_palette, labels=group_labels) +
      theme_bw() +
      theme(strip.background = element_rect(fill="white"),
            axis.text.x = element_text(angle=45, vjust=1, hjust=1),
            axis.text = element_text(size=14),
            axis.title = element_text(size=14),
            strip.text.x = element_text(size = 20, colour = 'white'),
            strip.background.x = element_rect(fill="black"),
            legend.position="none")
  })
  
  # display plot
  output$abundance <- renderPlot({
      print(plot_reactive())
  })
  
  
  # download plot
  output$download <- downloadHandler(
    filename = "transcript_abundance.png",
    content  = function(file) {
      ggsave(
        file,
        plot = plot_reactive(),
        height = ceiling(length(isolate(
          myValues$genelist
        )) / num_col) * 4,
         width = isolate(myValues$width)
      )
    }
  )
  
  # restart session (clear all input data)
  observeEvent(input$clear, {
    session$reload()
  })

  }

shinyApp(ui = ui, server = server)
