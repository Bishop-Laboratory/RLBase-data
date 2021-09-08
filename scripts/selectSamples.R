#' App used for selecting samples for the ML model to train upon
#' Usage: `Rscript selectSamples.R configtsv outfile host port`

library(shiny)
library(tidyverse)
library(DT)
library(ggprism)

# Magic
RLFSRDAFILE <- "rlbase-data/misc/rlfsRes.rda"
MANIFEST_FINAL <- "rlbase-data/rlbase_manifest_final.tsv"

# Load all previously-selected samples 
# This ensures selections are saved even if session ends
load('misc-data/todiscard.rda')

if (! interactive()) {
  # Get args
  args <- commandArgs(trailingOnly = TRUE)
  
  # Get output dir
  configtsv <- args[1]
  host <- args[2]
  port <- as.numeric(args[3])
} else {
  configtsv <- "rlbase-data/rlpipes-out/config.tsv"
  host <- "0.0.0.0"
  port <- 4848
}


# Get RLFS data
load(RLFSRDAFILE)

# Get the manifest
manifest <- read_tsv(configtsv) 

# Wrangle for display
toDT <- manifest %>% 
  select(experiment, condType, mode, condition) %>%
  filter(experiment %in% names(rlfsRes),
         ! duplicated(experiment)) %>%
  arrange(condType, experiment) %>%
  unique() %>%
  mutate(
    # Check boxes are already checked if sample is in the todiscard.rda vector
    discard = paste0('<input type="checkbox" name="', experiment, '" value="', TRUE, '" ', ifelse(experiment %in% discard, "checked", ""), ' />')
  )
  

# UI panels
RLFS_panel <- function() {
  tagList(
    br(),
    fluidRow(
      column(
        width = 6,
        htmlOutput(outputId = "RLFSOutHTML")
      ),
      column(
        width = 6,
        plotOutput('zScorePlot')
      )
    ),
    fluidRow(
      column(
        width = 6,
        plotOutput('pValPlot')
      ),
      column(
        width = 6,
        plotOutput('FFTPlot')
      )
    )
  )
}

# UI Page
ui <- fluidPage(
  
  # App title ----
  titlePanel("Select samples"),
  p("Select the samples which do not match their 'condType' to 'discard'. Selections are saved automatically in 'misc-data/todiscard.rda'"),
  fluidRow(
    column(
      width = 6,
      fluidRow(
        column(
          width = 12,
          DTOutput("rlbaseSamples")
        )
      ),
      fluidRow(
        column(
          width = 12,
          verbatimTextOutput("sel")
        )
      )
    ),
    column(
      width = 6,
      RLFS_panel()
    )
  )
  
)


# Server function
server <- function(input, output, session) {
  
  observeEvent(input$save, {
    print("press")
  })
  
  output$rlbaseSamples <- renderDT({
    toDT %>% column_to_rownames("experiment")
  }, selection = 'single', escape=FALSE, server = FALSE, 
  options = list(dom = 't', paging = FALSE, ordering = FALSE, scrollY="600px"),
    callback = JS("table.rows().every(function(i, tab, row) {
            var $this = $(this.node());
            $this.attr('id', this.data()[0]);
            $this.addClass('shiny-input-checkboxgroup');
          });
          Shiny.unbindAll(table.table().node());
          Shiny.bindAll(table.table().node());")
  )
  
  output$sel <- renderPrint({
     discard <- sapply(toDT$experiment, function(i) input[[i]])
     discard <- names(which(sapply(discard, function(x) ! is.null(x))))
     save(discard, file = "misc-data/todiscard.rda")
     str(discard)
  })
  
  
  
  
  current_samp <- reactive({
    # Get selected row from datatable
    selectedRow <- ifelse(is.null(input$rlbaseSamples_rows_selected), 
                          1, 
                          input$rlbaseSamples_rows_selected)
    
    # Get current sample ID
    current_samp <- toDT %>%
      filter(dplyr::row_number() == selectedRow) %>%
      pull(experiment)
    
    current_samp
  })
  
  
  output$zScorePlot <- renderPlot({
    
    # Get LZ
    lz <- rlfsRes %>%
      pluck(current_samp(),
            "Z-scores", 
            "regioneR::numOverlaps")
    
    # Plot
    data.frame("zscore" = lz$shifted.z.scores,
               "shift" = lz$shifts) %>%
      ggplot(aes(y = zscore, x = shift)) +
      geom_vline(color = "firebrick", xintercept = 0, linetype = "dashed") +
      geom_line(size = 1) +
      ggtitle('ZScore around RLFS', subtitle = current_samp()) +
      ylab("Peak Enrichment (Z-Score)") +
      scale_y_continuous(expand = c(0,0)) +
      xlab("Distance to RLFS (bp)") +
      theme_prism(base_size = 15)
  }) %>%
    bindCache(current_samp())
  
  output$RLFSOutHTML <- renderUI({
    
    pval <- rlfsRes %>%
      pluck(current_samp(),
            "perTestResults", 
            "regioneR::numOverlaps",
            "pval") 
    
    # Get the data for this sample
    tagList(
      p(
        style = paste0("color: ", ifelse(pval < .05, "blue", "red")),
        round(pval, 5)
      )
    )
  })
  
  output$FFTPlot <- renderPlot({
    
    # Get LZ
    lz <- rlfsRes %>%
      pluck(current_samp(),
            "Z-scores", 
            "regioneR::numOverlaps")
    
    # Plot
    data.frame("fftval" = Re(fft(lz$shifted.z.scores)),
               "freq" = seq(lz$shifts)) %>%
      ggplot(aes(y = fftval, x = freq)) +
      geom_hline(color = "firebrick", yintercept = 0, linetype = "dashed") +
      geom_line(size = 1) +
      ggtitle('Fourier Transform of ZScore around RLFS', subtitle = current_samp()) +
      ylab("FT-ZScore (Real Component)") +
      xlab("Relative Frequency") +
      theme_prism(base_size = 15)
  }) %>%
    bindCache(current_samp())
  
  
  output$pValPlot <- renderPlot({
    # Get LZ
    pt <- rlfsRes %>%
      pluck(current_samp(),
            "perTestResults", 
            "regioneR::numOverlaps")
    
    # Plot
    regioneR:::plot.permTestResults(pt)
    
  }) %>%
    bindCache(current_samp())
}

app <- shinyApp(ui, server)
runApp(app, port = port, host=host)
