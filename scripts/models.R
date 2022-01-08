#' App used for selecting samples for the ML model to train upon
#' Usage: `Rscript selectSamples.R configtsv outfile host port`

library(shiny)
library(tidyverse)
library(DT)
library(ggprism)
library(shinyWidgets)

# Magic
RLFSRDAFILE <- "rlbase-data/misc/rlfsRes.rda"
MANIFEST_FINAL <- "rlbase-data/rlpipes-out/config.tsv"
RLFSDIR <- "rlbase-data/rlpipes-out/rlfs_rda/"
TODISCARD <- 'misc-data/model/todiscard.rda'
PREVDONE <- "misc-data/model/prevdone.rda"


# Load all previously-selected samples 
# This ensures selections are saved even if session ends
load(TODISCARD)
load(PREVDONE)

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
manifest <- read_tsv(configtsv, show_col_types = FALSE) 

# Clean cols
if ("label.x" %in% colnames(manifest)) {
  to_rep <- which(manifest$label.x != manifest$label.y)
  manifest$label <- manifest$label.x
  manifest$label[to_rep] <- "NEG"
  manifest$condition[to_rep] <- "ActD"
  manifest <- unique(dplyr::select(manifest, -label.x, -label.y))
}

# Check for previous
if (file.exists(file.path(RLFSDIR, "manifest.csv"))) {
  prev <- read_csv(file.path(RLFSDIR, "manifest.csv"))
} else {
  prev <- tibble(id = "")
}

# Wrangle for display
toDT <- manifest %>% 
  dplyr::select(experiment, label) %>%
  dplyr::filter(experiment %in% names(rlfsRes),
         ! duplicated(experiment)) %>%
  unique() %>%
  mutate(
    # Check boxes are already checked if sample is in the todiscard.rda vector
    discard = paste0('<input type="checkbox" name="', experiment, 
                     '" value="', TRUE, '" ',
                     ifelse(experiment %in% discard, "checked", ""), ' />'),
    previously_evaluated = ifelse(! (experiment %in% prev$id | experiment %in% {{ discard }}),
                                 "<p style='color:orange'>No</p>",
                                 "<p style='color:blue'>Yes</p>")
  ) %>%
  arrange(desc(previously_evaluated), desc(label), experiment) %>%
  dplyr::rename(label = label)
  

# UI panels
RLFS_panel <- function() {
  tagList(
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    fluidRow(
      column(
        width = 10,
        plotOutput('zScorePlot', height = "500px")
      )
    )
  )
}

# UI Page
ui <- fluidPage(
  
  # App title ----
  titlePanel("Select samples to exclude from ML training"),
  fluidRow(
    column(
      width = 8,
      hr(),
      span(strong("Purpose: "), "The ML model used in RLSuite will need to adapt ",
           "to new studies, datasets, and techniques to remain useful. When new data ",
           "is added to RLBase, the model should be rebuilt so that it continually improves over time.",
           " However, the model cannot be supplied blindly with training data because",
           " many samples do not match their author-assigned labels",
           ". For instance, RNaseH1 treatment (a 'NEG' label)",
           " often fails to remove R-loops effectively and the sample can strongly resemble a 'POS' sample (e.g., ERX3974965	- ",
           strong(a(href="https://rlbase-data.s3.amazonaws.com/misc/example_failed_RNH.png", target="_blank", "link")),
           "). Therefore, human intervention is needed to ", strong("select the samples which clearly do not match their label"), 
           " so they may be exluded from the training data."),
      hr(),
       span(strong("Usage: "), " Select the samples which clearly do not match their",
            " 'label' to 'discard' (i.e., the ML model will not train on them). ",
            strong("Selections are saved automatically."),
            paste0(" New samples which have not yet",
                   " been evaluated are indicated in the 'previously_evaluated' column.",
                   " When ready, please click 'Build model' to generate a new model and HTML report.")), 
      hr(),
      h4("Label guide"),
      br(),
      span(strong("POS:"), " A positive sample for R-loop mapping (e.g., typical S9.6 DRIP-Seq)."),
      br(),
      span(strong("NEG:"), " A sample which does not map R-loops faithfully (e.g., S9.6 DRIP-Seq + RNaseH1 treatment or Input control)."),
      br(),
      br(),
      span("See examples of POS and NEG ", strong(a(href="https://rlbase-data.s3.amazonaws.com/misc/examples_for_select_samples.png", target="_blank", "here."))),
      hr(),
      span(paste0("*NOTE: Some samples will be particularly difficult to judge because ",
                  "they share many characteristics with the ideal examples but do not ",
                  "look quite correct. An example is ERX2277510 and ERX2277511.",
                  " In these cases, it is better to not exclude them so that the ",
                  "model can train on these more difficult-to-judge samples.")),
      hr()
    )
  ),
  fluidRow(
    column(
      width = 1,
      shiny::actionButton(inputId = "done", label = "Build model", class="btn btn-danger")
    ),
    column(
      width = 2,
      textOutput("donetime")
    )
  ),
  fluidRow(
    column(
      width = 1,
      a(h4("Report", class = "btn btn-info action-button" , 
           style = "fontweight:600"), target = "_blank",
        href = "https://rlbase-data.s3.amazonaws.com/misc/model/FFT-classifier.html")
    ),
    column(
      width = 1,
      a(h4("Feature Prep Model", class = "btn btn-info action-button" , 
           style = "fontweight:600"), target = "_blank",
        href = "https://rlbase-data.s3.amazonaws.com/misc/model/prepFeatures.rda")
    ),
    column(
      width = 1,
      a(h4("Classifier", class = "btn btn-info action-button" , 
           style = "fontweight:600"), target = "_blank",
        href = "https://rlbase-data.s3.amazonaws.com/misc/model/fftModel.rda")
    )
  ),
  fluidRow(
    column(
      width = 6,
      fluidRow(
        column(
          width = 12,
          hr(),
          h4("Discards"),
          verbatimTextOutput("sel"),
          hr(),
        )
      ),
      fluidRow(
        column(
          width = 12,
          DTOutput("rlbaseSamples"),
          hr(),
          br()
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
  
  dtdat <- reactiveVal(value = toDT)
  
  
  observeEvent(input$save, {
    print("press")
  })
  
  output$rlbaseSamples <- renderDT({
    dtdat() %>% column_to_rownames("experiment")
  }, selection = list(mode = 'single', selected = 1), escape=FALSE, server = FALSE, 
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
     message(discard)
     message(length(discard))
     save(discard, file = TODISCARD)
     str(discard)
  })
  
  donetime <- reactiveVal(prevdone)
  
  output$donetime <- reactive({
    donetime()
  })
  
  observeEvent(input$done, {
    # Create a Progress object
    progress <- shiny::Progress$new()
    
    progress$set(message = "Writing manifest.", value = .1)
    on.exit(progress$close())
    
    Sys.sleep(2)
    # Get samples to discard
    load(TODISCARD)
    
    # Get the manifest
    manifest <- read_tsv(MANIFEST_FINAL, show_col_types = FALSE)
    if ("label.x" %in% colnames(manifest)) {
      to_rep <- which(manifest$label.x != manifest$label.y)
      manifest$label <- manifest$label.x
      manifest$label[to_rep] <- "NEG"
      manifest$condition[to_rep] <- "ActD"
      manifest <- unique(dplyr::select(manifest, -label.x, -label.y))
    }
    manifestModel <- manifest %>% 
      dplyr::filter(experiment %in% names(rlfsRes),
             ! duplicated(experiment),
             ! experiment %in% discard) %>%
      arrange(label, experiment) %>%
      mutate(
        group = label,
        filename = paste0(experiment, "_", genome, ".rlfs.rda")
      ) %>%
      dplyr::select(id=experiment, group, filename) %>% unique()
    
    write_csv(manifestModel, file = file.path(RLFSDIR, "manifest.csv"))
    
    progress$set(message = "Building model.", detail = "takes ~1-2 minutes", value = .5)
    rmarkdown::render(input = "scripts/FFT-clasiffier.Rmd", output_file = "../misc-data/model/FFT-classifier.html")
    
    progress$set(message = "Uploading to AWS.", detail="HTML", value = .9)
    aws.s3::put_object(
      file = "misc-data/model/FFT-classifier.html", 
      object = "misc/model/FFT-classifier.html", 
      bucket = "s3://rlbase-data/"
    )
    progress$set(message = "Uploading to AWS.", detail="Models", value = .9)
    aws.s3::put_object(
      file = "misc-data/model/fftModel.rda", 
      object = "misc/model/fftModel.rda", 
      bucket = "s3://rlbase-data/"
    )
    aws.s3::put_object(
      file = "misc-data/model/prepFeatures.rda", 
      object = "misc/model/prepFeatures.rda", 
      bucket = "s3://rlbase-data/"
    )
    Sys.sleep(1)
    
    progress$close()
    prevdone <- paste0("Last run at: ", Sys.time())
    save(prevdone, file = PREVDONE)
    donetime(prevdone)
    newdtdat <- dtdat() %>%
      mutate(
        previously_evaluated = "<p style='color:blue'>Yes</p>"
      )
    dtdat(newdtdat)
    
    shinyWidgets::sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Links Updated",
      type = "success"
    )
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
    
    # Pvql
    pval <- rlfsRes %>%
      pluck(current_samp(),
            "perTestResults", 
            "regioneR::numOverlaps",
            "pval") 
    
    # Plot
    data.frame("zscore" = lz$shifted.z.scores,
               "shift" = lz$shifts) %>%
      ggplot(aes(y = zscore, x = shift)) +
      geom_vline(color = "firebrick", xintercept = 0, linetype = "dashed") +
      geom_line(size = 1) +
      ggtitle('ZScore around RLFS', subtitle = current_samp()) +
      ylab("Peak Enrichment (Z-Score)") +
      labs(caption = paste0("pval (permTes) ~ ", round(pval, 5))) +
      scale_y_continuous(expand = c(0,0)) +
      xlab("Distance to RLFS (bp)") +
      theme_prism(base_size = 15) +
      theme(plot.caption = element_text(
        face = "bold",
        colour = ifelse(pval < .05, "blue", "red"), size = 20
      )) 
  }) %>%
    bindCache(current_samp())
}

app <- shinyApp(ui, server)
runApp(app, port = port, host=host)
