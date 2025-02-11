library(shiny)
library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)
library(DT)
library(tidyr)
library(plotly)

addResourcePath("/assets", file.path(getwd(), "www"))

# Load Data for both apps
combined_data <- read.csv("./data_folder/kir_1uM_final.csv", row.names = 2)
mutant_wild_kinases <- combined_data %>%
  select(-CAS, -SMILES, -cSMILES, -`Dose..µM.`)
mutant_wild_kinases[is.na(mutant_wild_kinases)] <- 100
column_choices <- sub("\\.", "(", colnames(mutant_wild_kinases))  # Replace the first dot with '(' due to R code problem where ( and ) becomes .
column_choices <- gsub("\\.", ")", column_choices) # Replace the last dot with ')' due to R code problem where ( ) becomes .

colnames(mutant_wild_kinases) <- column_choices
gini_scores <- read.csv("./data_folder/gini_coefficients.csv", row.names = 1)


all_kinases <- read.csv('./data_folder/all_kinases_all_lineages_zeroes.csv')
new_data <- read.csv("./data_folder/kir_1uM_final.csv")
column_choices2 <- sub("\\.", "(", colnames(new_data))  # Replace the first dot with '(' due to R code problem where ( and ) becomes .
column_choices2 <- gsub("\\.", ")", column_choices2) # Replace the last dot with ')' due to R code problem where ( ) becomes .

colnames(new_data) <- column_choices2
new_wt_data <- new_data[, 1:(5+409)]
new_wt_data[, 6:ncol(new_wt_data)] <- 100 - new_wt_data[, 6:ncol(new_wt_data)]
wt_kinases <- colnames(new_wt_data)[6:ncol(new_wt_data)]
mutant_kinases <- colnames(mutant_wild_kinases)[410:ncol(mutant_wild_kinases)]

# Function to calculate KISS score (from cancer_app2.R)
calculate_KISS <- function(row, kinases) {
  kinases <- kinases[kinases %in% wt_kinases]
  if (length(kinases) == 0) return(NA)
  
  if (all(is.na(row[kinases]))) return(NA)
  
  if (length(row[kinases]) > 1) {
    on_target <- exp(mean(log(row[kinases]), na.rm = TRUE))
  } else {
    on_target <- row[kinases]
  }
  
  off_target_cols <- row[!names(row) %in% kinases]
  IP_1 <- mean(off_target_cols, na.rm = TRUE)
  sy2 <- var(off_target_cols, na.rm = TRUE)
  sz2 <- var(row, na.rm = TRUE)
  slo2 <- (max(off_target_cols, na.rm = TRUE)^2) / (2 * sum(!is.na(off_target_cols)))
  
  IP_2 <- (sy2 - slo2) / (sz2 - slo2)
  IP <- (IP_1 + 100 * IP_2^5) / 2
  KISS <- on_target - IP
  return(KISS)
}

# Function for looking up kinases (from cancer_app2.R)
lookup <- function(kinase) {
  # If no kinase input, return NULL
  if (is.null(kinase) || kinase == "") return(NULL)
  
  # Find the first partial match for the kinase
  matched_kinase <- wt_kinases[grep(tolower(kinase), tolower(wt_kinases))][1]
  
  # If no match found, return NULL
  if (is.na(matched_kinase)) {
    return(NULL)
  }
  
  # Create a dataframe with relevant columns (CAS, Compound, Dose, Kinase Inhibition)
  df <- data.frame(
    CAS = new_wt_data$CAS,
    Compound = new_wt_data$Compound,
    Dose = new_wt_data$`Dose()µM)`,  # Assuming correct format here
    Kinase_Inhibition = new_wt_data[[matched_kinase]],
    stringsAsFactors = FALSE
  )
  
  # If all values in Kinase_Inhibition are NA, return NULL
  if (all(is.na(df$Kinase_Inhibition))) return(NULL)
  
  # Calculate KISS score
  df$KISS <- apply(new_wt_data[, 6:ncol(new_wt_data)], 1, function(row) {
    score <- calculate_KISS(row, kinases = matched_kinase)
    if (is.nan(score)) return(NA)
    return(score)
  })
  
  # Filter rows where Kinase_Inhibition and KISS are not NA
  df <- df[!is.na(df$Kinase_Inhibition) & !is.na(df$KISS), ]
  
  # If the filtered dataframe is empty, return NULL
  if (nrow(df) == 0) return(NULL)
  
  return(df)
}

# Function to generate radar plot
radar_plot_for_kinases <- function(wild_df, kinase, input_kinase) {
  num_drugs <- nrow(wild_df)
  drug_names <- rownames(wild_df)
  
  angles <- seq(0, 2 * pi, length.out = num_drugs + 1)
  
  inhibition_levels <- c(25, 50, 75, 100)
  
  values <- 100 - wild_df[[kinase]]
  values <- c(values, values[1])
  
  par(mar = rep(1, 4)) 
  plot.new() 
  plot.window(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1)
  
  # Draw circular grid for inhibition levels
  for (level in inhibition_levels) {
    lines(cos(angles) * level / 100, sin(angles) * level / 100, lty = 2, col = 'darkgrey')
    text(0, level / 100, paste0('≤ ', level, '% Inhibition'), col = 'darkgrey', cex = 1.2)
  }
  
  # Draw radial lines for inhibition values with varying alpha and length
  for (i in 1:num_drugs) {
    # Adjust inhibition values to thresholds
    if (values[i] <= 25) {
      x_pos <- cos(angles[i]) * 0.25
      y_pos <- sin(angles[i]) * 0.25
      segments(0, 0, x_pos, y_pos, col = rgb(30, 144, 255, maxColorValue = 255), lwd = 0.5)  
    } else if (values[i] <= 50) {
      x_pos <- cos(angles[i]) * 0.5 
      y_pos <- sin(angles[i]) * 0.5
      segments(0, 0, x_pos, y_pos, col = rgb(30, 144, 255, maxColorValue = 255), lwd = 1) 
    } else if (values[i] <= 75) {
      x_pos <- cos(angles[i]) * 0.75 
      y_pos <- sin(angles[i]) * 0.75
      segments(0, 0, x_pos, y_pos, col = rgb(30, 144, 255, maxColorValue = 255), lwd = 1.5)  
    } else {
      x_pos <- cos(angles[i]) * 1 
      y_pos <- sin(angles[i]) * 1
      segments(0, 0, x_pos, y_pos, col = rgb(30, 144, 255, maxColorValue = 255), lwd = 4) 
    }
  }
  
  # Add kinase name in the center
  text(0, 0, labels = input_kinase, cex = 2.0, col = "black", font = 2)
  
  # Add drug labels for values with inhibition < 10
  for (i in 1:num_drugs) {
    if (wild_df[[kinase]][i] < 10) {
      label_pos_x <- cos(angles[i]) * (values[i] / 100 + 0.15)
      label_pos_y <- sin(angles[i]) * (values[i] / 100 + 0.15)
      text(label_pos_x, label_pos_y, labels = drug_names[i], cex = 1.0, col = "black")
    }
  }
}

# UI for merged app
ui <- fluidPage(
  # Custom style to remove background and box shadow
  tags$style(HTML("
    .container-fluid {
      background-color: transparent !important;
      box-shadow: none !important;
    }
  ")),
  
  tags$div(
    style = "display: flex; justify-content: space-between; align-items: center; width: 100%;",
    
    # Title (centered)
    div(
      titlePanel(div(HTML('<b>KIRHub: Kinase Inhibitor Repurposing Hub</b>'), style = "font-size: 52px; color: steelblue; text-align: center; width: 100%;"))
    ),
    
    # Logo with link, positioned on the top-right and with space at the bottom
    div(
      tags$a(href = "https://research.fredhutch.org/gujral/en/lab-members.html", target = "_blank",
             tags$img(src = "assets/logo.png", height = "80px", style = "position: absolute; right: 20px; top: 10px; margin-bottom: 20px;"))  # Added margin-bottom
    )
  ),
  
  tags$br(),
  
  tags$style(HTML("
    .well {
      background-color: transparent !important;
      box-shadow: none !important;
    }
  ")),
  
  sidebarLayout(
    sidebarPanel(
      h4(div(HTML('<b>Select an App to begin Kinase Analysis</b>'), style = "font-size: 35px; text-align: center; color: black;")),
      
      # Use divs with CSS styling to equally space out the buttons
      div(style = "display: flex; justify-content: center; gap: 15px;",  # Adjusted spacing
          actionButton("wild_button", HTML("Identifying Drug Candidates <br> to Inhibit Wild-Type Kinase(s)"), 
                       class = "btn btn-success", 
                       style = "width: 620px; font-size: 35px; padding: 12px 18px; text-align: center; 
                          text-decoration: underline; background-color: #5cb85c; color: white; 
                          border: 2px solid #5cb85c; border-radius: 5px; white-space: normal;"),
          
          actionButton("mutation_button", HTML("Identifying Drug Candidates <br> to Inhibit Mutant Kinase(s)"), 
                       class = "btn btn-primary", 
                       style = "width: 620px; font-size: 35px; padding: 12px 18px; text-align: center; 
                          text-decoration: underline; background-color: #337ab7; color: white; 
                          border: 2px solid #337ab7; border-radius: 5px; white-space: normal;"),
          
          actionButton("lineage_button", HTML("Exploring Kinase Importance <br> Across Cancer Lineages"), 
                       class = "btn btn-info", 
                       style = "width: 620px; font-size: 35px; padding: 12px 18px; text-align: center; 
                          text-decoration: underline; background-color: #5bc0de; color: white; 
                          border: 2px solid #5bc0de; border-radius: 5px; white-space: normal;")
      ),
      width = 2.3
    ),
  
    mainPanel(
      # Introduction Section (Placed BELOW the buttons)
      tags$div(
        style = "text-align: center; max-width: 1500px; margin-left: 200px; padding: 10px;",
        tags$p(HTML("
          <b>KIRHub</b> is an interactive web-based platform designed for researchers and clinicians 
          to explore the relationship between clinically-approved kinase inhibitors and their effects on 
          both wild-type and oncogenic kinase variants. The application is fully browser-based, 
          requiring no data uploads, ensuring ease of use and accessibility. 
          Click on the logo on the top right to reach out regarding any questions or inquiries.<br><br>
        "), style = "font-size: 18px; color: black;")
      ),
      
      uiOutput("app_ui"),
      width = 15
    )
  )
)

# Server logic for merged app
server <- function(input, output, session) {
  
  # Reactive values to keep track of which app is active
  active_app <- reactiveVal("")
  
  # Observe when the Mutation Combinations button is clicked
  observeEvent(input$mutation_button, {
    active_app("mutation")
  })
  
  # Observe when the Wild Kinase button is clicked
  observeEvent(input$wild_button, {
    active_app("wild")
  })
  
  # Observe when the Cancer Lineages button is clicked
  observeEvent(input$lineage_button, {
    active_app("lineage")
  })
  
  # UI output for the selected app
  output$app_ui <- renderUI({
    tags$head(
      tags$style(HTML("
      /* Set uniform font properties across the UI */
      * {
        font-family: Arial, sans-serif;
        font-size: 16px;
        color: #333;
      }
      /* Override specific text styles if needed */
      h4, h3, .btn, .form-control, .dataTables_wrapper {
        font-size: 20px;
      }
    "))
    )
    
    if (active_app() == "mutation") {
      
      ### UI FROM mutations_app1.R (Mutation Combinations)
      
      fluidPage(
        titlePanel(div(HTML('<b>Identifying Drug(s) to Inhibit Mutant Kinase(s)</b>'), style = "font-size: 30px; text-align: center; color: steelblue;")),
        
        sidebarLayout(
          sidebarPanel(
            h4(div(HTML('<b>Step 1: Select Kinase / Mutation</b>'), style = "font-size: 20px; text-align: left; color:steelblue")),
            selectizeInput(
              "first_mutation",
              label = NULL,
              choices = c("", mutant_kinases),
              multiple = FALSE,
              selected = "",
              options = list(placeholder = 'Select kinase / mutation of interest...')
            ),
            
            h4(div(HTML('<b>Step 2: Select Second Kinase / Mutation in Combination (Optional)</b>'), style = "font-size: 20px; text-align: left; color:steelblue")),
            selectizeInput(
              "second_mutation",
              label = NULL,
              choices = c("", mutant_kinases),
              multiple = FALSE,
              selected = "None",
              options = list(placeholder = 'Select kinase mutation in combination with above (optional)...')
            ),
            
            tags$div(
              style = "margin-top: 20px; padding: 10px; border: 1px solid #ddd; background-color: white;",
              tags$h4(div(HTML('<b>User Instructions:</b>'), style = "text-align: left; color:steelblue; font-size: 20px;")),
              tags$ol(
                tags$li(div(HTML("<b> Select a specific mutant kinase </b> for which you would like to identify an FDA-approved inhibitor (e.g., BRAF V600E, FGFR2 fusions, ALK fusions), then search for approved drugs that effectively target this mutation.")), style = "font-size: 15px;"),
                HTML("<br>"),
                tags$li(div(HTML("<b> (Optional) </b> Select an additional mutation to find FDA-approved drugs that inhibit the combination of mutations. If none, leave this blank.")), style = "font-size: 15px;"),
                tags$hr()
              )
            ),
            width = 3
          ),
          
          mainPanel(
            tags$div(
              h4(uiOutput("selected_mutations_text")),
              style = "font-size: 16px; color: steelblue"
            ),
            
            dataTableOutput("drug_table"),
            
            uiOutput("conditional_download_and_description"),
            
            tags$hr(),  # Horizontal line separator
            tags$br(), tags$br(),
            
            # Radar plot UI - Shows 1 or 2 radar plots based on selected kinases
            conditionalPanel(
              condition = "input.first_mutation != ''",
              plotOutput("radarPlot1", height = "600px"),
              tags$br(),
              downloadButton("downloadRadarPlot1", "Download Radar Plot 1 (.png)"),
              tags$hr()
            ),
            
            conditionalPanel(
              condition = "input.second_mutation != '' && input.second_mutation != 'None'",
              plotOutput("radarPlot2", height = "600px"),
              tags$br(),
              downloadButton("downloadRadarPlot2", "Download Radar Plot 2 (.png)"),
              tags$hr()
            ),
            
            width = 9
          )
        )
      )
      
      ### END OF UI FROM mutations_app1.R
      
    } else if (active_app() == "wild") {
      
      ### UI FROM Wild Kinases Only App with Paralogs
      
      fluidPage(
        titlePanel(div(HTML('<b>Identifying Drug(s) to Inhibit Wild-Type Kinases</b>'), style = "font-size: 30px; text-align: center; color: steelblue;")),
        
        sidebarLayout(
          sidebarPanel(
            h4(div(HTML('<b>Step 1: Select Kinase </b>'), style = "font-size: 20px; text-align: left; color:steelblue")),
            selectizeInput(
              "first_mutation",
              label = NULL,
              choices = c("", wt_kinases),
              multiple = FALSE,
              selected = "",
              options = list(placeholder = 'Select kinase of interest...')
            ),
            
            h4(div(HTML('<b>Step 2: Select Second Kinase / Paralog in Combination (Optional)</b>'), style = "font-size: 20px; text-align: left; color:steelblue")),
            selectizeInput(
              "second_mutation",
              label = NULL,
              choices = c("", wt_kinases),
              multiple = FALSE,
              selected = "None",
              options = list(placeholder = 'Select kinase in combination with above (optional)...')
            ),
            
            tags$div(
              style = "margin-top: 20px; padding: 10px; border: 1px solid #ddd; background-color: white;",
              tags$h4(div(HTML('<b>User Instructions:</b>'), style = "text-align: left; color:steelblue; font-size: 20px;")),
              tags$ol(
                tags$li(div(HTML("<b> Select a specific kinase </b> for which you would like to identify an FDA-approved inhibitor (e.g., ABL1).")), style = "font-size: 15px;"),
                HTML("<br>"),
                tags$li(div(HTML("<b> (Optional) </b> Select an additional kinase to find FDA-approved drugs that inhibit the combination kinases. If none, leave this blank.")), style = "font-size: 15px;"),
                tags$hr()
              )
            ),
            width = 3
          ),
          
          mainPanel(
            tags$div(
              h4(uiOutput("selected_wild_text")),
              style = "font-size: 16px; color: steelblue"
            ),
            
            dataTableOutput("drug_table"),
            
            uiOutput("conditional_download_and_des"),
            
            tags$hr(),  # Horizontal line separator
            tags$br(), tags$br(),
            
            # Radar plot UI - Shows 1 or 2 radar plots based on selected kinases
            conditionalPanel(
              condition = "input.first_mutation != ''",
              plotOutput("radarPlot1", height = "600px"),
              tags$br(),
              downloadButton("downloadRadarPlot1", "Download Radar Plot 1 (.png)"),
              tags$hr()
            ),
            
            conditionalPanel(
              condition = "input.second_mutation != '' && input.second_mutation != 'None'",
              plotOutput("radarPlot2", height = "600px"),
              tags$br(),
              downloadButton("downloadRadarPlot2", "Download Radar Plot 2 (.png)"),
              tags$hr()
            ),
            
            # KISS Score Scatter Plots - Shows 1 or 2 plots based on selected kinases
            conditionalPanel(
              condition = "input.first_mutation != ''",
              plotOutput("kissPlot1", height = "600px"),
              tags$br(),
              downloadButton("downloadKissPlot1", "Download KISS Plot 1 (.png)"),
              tags$hr()
            ),
            
            conditionalPanel(
              condition = "input.second_mutation != '' && input.second_mutation != 'None'",
              plotOutput("kissPlot2", height = "600px"),
              tags$br(),
              downloadButton("downloadKissPlot2", "Download KISS Plot 2 (.png)"),
              tags$hr()
            ),
            
            width = 9
          )
        )
      )
      
      ### END of UI from wild kinases app.R
      
    } else if (active_app() == "lineage") {
      
      ### UI FROM cancer_app2.R (Cancer Lineages)
      
      fluidPage(
        titlePanel(div(HTML('<b>Exploring Kinase Importance Across Cancer Lineages</b>'), style = "font-size: 30px; text-align: center; color: steelblue;")),
        
        sidebarLayout(
          sidebarPanel(
            h4(div(HTML('<b>Step 1: Select a Kinase</b>'), style = "font-size: 20px; text-align: left; color: steelblue;")),
            selectizeInput(
              "kinase",
              label = NULL,
              choices = c("", colnames(all_kinases)[4:ncol(all_kinases)]),
              multiple = FALSE,
              selected = NULL,
              options = list(placeholder = 'Select kinase of interest...')
            ),
            tags$div(
              style = "margin-top: 20px; padding: 10px; border: 1px solid #ddd; background-color: white;",
              tags$h4(div(HTML('<b>User Instructions:</b>'), style = "text-align: left; color: steelblue; font-size: 20px;")),
              tags$ol(
                tags$li(div(HTML("<b> Select the kinase </b> for which you would like to see the important cancer lineages and sub-lineages. Additionally the user can obtain a plot of the KISS Score vs. Inhibition for the kinase to see which FDA Approved drugs best target this kinase"), style = "font-size: 15px;")),
                tags$br(),
                tags$li(div(HTML("<b>Note: 'Important Kinase' is defined </b> as  a kinase that has a gene effect score of less than -0.5 in a particular cancer lineage and its sub-lineages. These gene effect scores data is obtained from the <a href='https://depmap.org/portal/' target='_blank'>DepMap Portal website</a>."), style = "font-size: 15px;")),
                tags$br(),
                tags$hr()
              )
            ),
            width = 3
          ),
          
          mainPanel(
            tags$div(
              h4(uiOutput("selected_kinase_text")),
              style = "font-size: 16px; color: steelblue"
            ),
            DTOutput("lineageTable"),
            tags$br(),
            conditionalPanel(
              condition = "input.kinase != ''",
              downloadButton("download_table", "Download Table of Results (.csv)")
            ),
            tags$hr(),
            
            conditionalPanel(
              condition = "input.kinase != ''",
              tags$div(
                style = "margin-top: 20px;",
                h4(HTML("<b>Column Descriptions</b>"), style = "color: steelblue; font-size: 20px;"),
                tags$ul(
                  tags$li(HTML("<b>Lineage 1:</b> The primary cancer lineage the kinase is critical in."), style = "font-size: 15px;"),
                  tags$li(HTML("<b>Lineage 2:</b> The sub lineage the kinase is critical in."), style = "font-size: 15px;"),
                  tags$li(HTML("<b>Lineage 3:</b> The most specific sub lineage the kinase is critical in."), style = "font-size: 15px;"),
                  tags$li(HTML("<b>Important Counts:</b> Number of cancer cell lines of the particular lineage, the kinase is critical in."), style = "font-size: 15px;"),
                  tags$li(HTML("<b>Total Counts:</b> Total number of cancer cell lines where the kinase is observed in."), style = "font-size: 15px;"),
                  tags$li(HTML("<b>Percentage of Important Counts:</b> Number of cancer cell lines where the kinase is critical / Total number of cancer cell lines where kinase is observed (%)"), style = "font-size: 15px;")
                )
              )
            ),
            tags$hr(),
            
            # # Radar plot UI (added)      
            # plotOutput("radarPlot", height = "600px"),  # Set radar plot height to 600px
            # tags$br(),
            # conditionalPanel(
            #   condition = "input.kinase != ''",
            #   downloadButton("downloadRadarPlot", "Download Radar Plot (.png)")
            # ),
            # tags$hr(),
            # 
            # 
            # # New KISS Score Plot Panel
            # plotOutput("kissPlot", height = "600px"),
            # tags$br(),
            # conditionalPanel(
            #   condition = "input.kinase != '' && output.kissPlot !== null",
            #   downloadButton("download_kiss_plot", "Download KISS Plot (.png)")
            # ),
            # tags$hr(),
            
            tags$br(), tags$br(),
            plotOutput("barPlot", height = "600px"),
            tags$br(),
            conditionalPanel(
              condition = "input.kinase != ''",
              downloadButton("download_plot", "Download Plot (.png)")
            ),
            conditionalPanel(
              condition = "input.kinase == ''"
            ),
            width = 9
          )
        )
      )
      
      ### END OF UI FROM cancer_app2.R
      
    } # else {
    #   tagList(
    #     tags$br(), tags$br(),
    #     h4(div("Please select an app to begin.", style = "text-align: center; font-weight: bold; font-size: 18px;"))
    #   )
    # }
  })
  
  # Server logic for Mutation Combinations app (from mutations_app1.R)
  observeEvent(active_app(), {
    if (active_app() == "mutation") {
      
      ### SERVER LOGIC FROM mutations_app1.R
      
      # Reset the selectizeInput to start empty
      observe({
        if (is.null(input$first_mutation)) {
          updateSelectizeInput(session, "first_mutation", selected = "")
        }
      })
      
      observeEvent(input$first_mutation, {
        # Get the current first mutation selection
        selected_first_mutation <- input$first_mutation
        
        # Update the choices for the second mutation, excluding the selected first mutation
        updateSelectizeInput(session, "second_mutation",
                             choices = c("None", colnames(mutant_wild_kinases[mutant_kinases])[colnames(mutant_wild_kinases[mutant_kinases]) != selected_first_mutation]),
                             selected = "None")
      })
      
      # Reactive data to handle the filtered table
      filtered_data <- reactive({
        selected_first_mutation <- input$first_mutation
        selected_second_mutation <- input$second_mutation
        
        # threshold <- 80  # Previously used inhibition threshold
        
        if (!is.null(selected_first_mutation) && selected_first_mutation != "") {
          first_mutation_columns <- selected_first_mutation
          
          # Get inhibition values for all drugs targeting the selected kinase
          drug_inhibition_values <- mutant_wild_kinases[, first_mutation_columns, drop = FALSE]
          
          # Convert to a data frame and retain row names (drug names)
          drug_inhibition_df <- data.frame(
            Kinase_Inhibitor = rownames(drug_inhibition_values),
            Inhibition = 100 - as.numeric(drug_inhibition_values[, first_mutation_columns]),
            stringsAsFactors = FALSE
          )
          
          # Remove NAs and select the top 10 drugs based on inhibition (sorted in descending order)
          drug_inhibition_df <- drug_inhibition_df[!is.na(drug_inhibition_df$Inhibition), ]
          drug_inhibition_df <- drug_inhibition_df[order(-drug_inhibition_df$Inhibition), ]
          top_10_drugs <- head(drug_inhibition_df, 10)
          
          inhibiting_drugs_set <- top_10_drugs$Kinase_Inhibitor
          
          # --- Old threshold-based logic (commented out) ---
          # inhibiting_drugs_set <- rownames(mutant_wild_kinases)[apply(mutant_wild_kinases[, first_mutation_columns, drop = FALSE], 1, function(row) any(as.numeric(row) < threshold))]
          
          # If a second mutation is selected, find the common inhibitors for both
          if (!is.null(selected_second_mutation) && selected_second_mutation != "" && selected_second_mutation != "None") {
            second_mutation_columns <- selected_second_mutation
            
            # Get inhibition values for second mutation
            second_drug_inhibition_values <- mutant_wild_kinases[, second_mutation_columns, drop = FALSE]
            
            second_drug_inhibition_df <- data.frame(
              Kinase_Inhibitor = rownames(second_drug_inhibition_values),
              Inhibition = 100 - as.numeric(second_drug_inhibition_values[, second_mutation_columns]),
              stringsAsFactors = FALSE
            )
            
            second_drug_inhibition_df <- second_drug_inhibition_df[!is.na(second_drug_inhibition_df$Inhibition), ]
            second_drug_inhibition_df <- second_drug_inhibition_df[order(-second_drug_inhibition_df$Inhibition), ]
            top_10_second_drugs <- head(second_drug_inhibition_df, 10)
            
            second_inhibiting_drugs <- top_10_second_drugs$Kinase_Inhibitor
            
            # Find drugs common to both kinase inhibition lists
            inhibiting_drugs_set <- intersect(inhibiting_drugs_set, second_inhibiting_drugs)
          }
          
          unique_inhibiting_drugs <- unique(inhibiting_drugs_set)
          gini_selected_drugs <- gini_scores[unique_inhibiting_drugs, "Gini", drop = FALSE]
          
          first_inhibition_col <- paste0("Inhibition (", selected_first_mutation, ") (%)")
          second_inhibition_col <- paste0("Inhibition (", selected_second_mutation, ") (%)")
          
          if (!is.null(selected_second_mutation) && selected_second_mutation != "" && selected_second_mutation != "None") {
            output_data <- data.frame(
              "Kinase Inhibitor" = unique_inhibiting_drugs,
              "Dose (um)" = rep("1.0", length(unique_inhibiting_drugs)),
              first_inhibition_col = 100 - as.numeric(mutant_wild_kinases[unique_inhibiting_drugs, selected_first_mutation]),
              second_inhibition_col = 100 - as.numeric(mutant_wild_kinases[unique_inhibiting_drugs, selected_second_mutation]),
              "Gini Score" = gini_selected_drugs$Gini,
              row.names = NULL
            )
            
            colnames(output_data)[1] <- "Kinase Inhibitor"
            colnames(output_data)[2] <- "Dose (μm)"
            colnames(output_data)[3] <- first_inhibition_col
            colnames(output_data)[4] <- second_inhibition_col
            colnames(output_data)[5] <- "Gini Score"
          } else {
            output_data <- data.frame(
              "Kinase Inhibitor" = unique_inhibiting_drugs,
              "Dose (um)" = rep("1.0", length(unique_inhibiting_drugs)),
              first_inhibition_col = 100 - as.numeric(mutant_wild_kinases[unique_inhibiting_drugs, selected_first_mutation]),
              "Gini Score" = gini_selected_drugs$Gini,
              row.names = NULL
            )
            colnames(output_data)[1] <- "Kinase Inhibitor"
            colnames(output_data)[2] <- "Dose (μm)"
            colnames(output_data)[3] <- first_inhibition_col
            colnames(output_data)[4] <- "Gini Score"
          }
          
          output_data <- output_data %>%
            arrange(desc(!!sym(first_inhibition_col)))
          
          return(output_data)
        } else {
          return(NULL)
        }
      })
      
      # Render the table
      output$drug_table <- renderDataTable({
        filtered_data()
      })
      
      # Conditional UI for download button and column descriptions
      output$conditional_download_and_description <- renderUI({
        if (!is.null(filtered_data()) && nrow(filtered_data()) > 0) {
          tagList(
            downloadButton("download_csv", "Download Table of Results (.csv)"),
            tags$br(),
            tags$hr(),
            tags$div(
              style = "margin-top: 20px;",
              h4(HTML("<b>Column Descriptions</b>"), style = "color:steelblue; font-size: 20px;"),
              tags$ul(
                tags$li(HTML("<b>Kinase Inhibitor:</b> The drug or compound used for inhibition."), style = "font-size: 15px;"),
                tags$li(HTML("<b>Dose (μm):</b> The dose (in micromolar) of the compound used in the study."), style = "font-size: 15px;"),
                tags$li(HTML("<b>Inhibition (kinase) (%):</b> The percentage inhibition of the selected kinase(s) by the drug."), style = "font-size: 15px;"),
                tags$li(HTML("<b>Gini Score:</b> A measure of selectivity for the drug. Higher values indicate greater selectivity."), style = "font-size: 15px;")
              )
            )
          )
        }
      })
      
      # CSV Download
      output$download_csv <- downloadHandler(
        filename = function() {
          all_mutations <- if (!is.null(input$second_mutation) && input$second_mutation != "None" && input$second_mutation != "") {
            paste(input$first_mutation, input$second_mutation, sep = "_")
          } else {
            input$first_mutation
          }
          paste(all_mutations, "_", "results_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(filtered_data(), file, row.names = FALSE)
        }
      )
      
      output$selected_mutations_text <- renderUI({
        if (!is.null(input$first_mutation) && input$first_mutation != "") {
          all_mutations <- if (!is.null(input$second_mutation) && input$second_mutation != "None" && input$second_mutation != "") {
            c(input$first_mutation, input$second_mutation)
          } else {
            input$first_mutation
          }
          
          HTML(paste('<b style="font-size: 16px;">Kinase Inhibitors for Selected Mutation(s):</b> ', 
                     paste(all_mutations, collapse = " and "), 
                     '<span style="font-size: 16px;"></span>'))
        } else {
          tagList(
            HTML('<span style="font-size: 16px;">Select kinase(s) of interest</span>'),
            tags$hr()
          )
        }
      })
      
      # Radar Plot Rendering (1 or 2 plots)
      output$radarPlot1 <- renderPlot({
        if (is.null(input$first_mutation) || input$first_mutation == "") return(NULL)
        
        title_text <- paste("Radar Plot for", input$first_mutation)  # Set title dynamically
        
        # Render Radar Plot with Title
        par(mar = c(3, 3, 3, 3))  # Adjust margins for title
        plot.new()
        title(main = title_text, col.main = "navy", font.main = 2, cex.main = 1.8)  # Set title properties
        
        radar_plot_for_kinases(mutant_wild_kinases, input$first_mutation, input$first_mutation)
      })
      
      output$radarPlot2 <- renderPlot({
        if (is.null(input$second_mutation) || input$second_mutation == "None") return(NULL)
        
        title_text <- paste("Radar Plot for", input$second_mutation)  # Set title dynamically
        
        # Render Radar Plot with Title
        par(mar = c(3, 3, 3, 3))  # Adjust margins for title
        plot.new()
        title(main = title_text, col.main = "navy", font.main = 2, cex.main = 1.8)  # Set title properties
        
        radar_plot_for_kinases(mutant_wild_kinases, input$second_mutation, input$second_mutation)
      })
      
      # Download Radar Plots
      output$downloadRadarPlot1 <- downloadHandler(
        filename = function() { paste(input$first_mutation, "_radar_plot.png", sep = "") },
        content = function(file) {
          png(file, width = 800, height = 800, res = 150)
          radar_plot_for_kinases(mutant_wild_kinases, input$first_mutation, input$first_mutation)
          dev.off()
        }
      )
      
      output$downloadRadarPlot2 <- downloadHandler(
        filename = function() { paste(input$second_mutation, "_radar_plot.png", sep = "") },
        content = function(file) {
          if (is.null(input$second_mutation) || input$second_mutation == "None") return(NULL)
          png(file, width = 800, height = 800, res = 150)
          radar_plot_for_kinases(mutant_wild_kinases, input$second_mutation, input$second_mutation)
          dev.off()
        }
      )
      
      ### END OF SERVER LOGIC FROM mutations_app1.R
      
    }
  })
  
  # Server logic for Wild Kinases app (from wild_kinases.R)
  observeEvent(active_app(), {
    if (active_app() == "wild") {
      
      ### SERVER LOGIC FROM wild_kinases.R
      
      # Reset the selectizeInput to start empty
      observe({
        if (is.null(input$first_mutation)) {
          updateSelectizeInput(session, "first_mutation", selected = "")
        }
      })
      
      observeEvent(input$first_mutation, {
        selected_mutation <- input$first_mutation
        
        # Update the choices for the second mutation, excluding the selected first mutation
        updateSelectizeInput(session, "second_mutation",
                             choices = c("None", wt_kinases[wt_kinases != selected_mutation]),  # Exclude first mutation
                             selected = "None")
      })
      
      # Reactive data to handle the filtered table
      filtered_data <- reactive({
        selected_first_mutation <- input$first_mutation
        selected_second_mutation <- input$second_mutation
        
        # threshold <- 80  # Previously used inhibition threshold
        
        if (!is.null(selected_first_mutation) && selected_first_mutation != "") {
          first_mutation_columns <- selected_first_mutation
          
          # Get inhibition values for all drugs targeting the selected kinase
          drug_inhibition_values <- mutant_wild_kinases[, first_mutation_columns, drop = FALSE]
          
          # Convert to a data frame and retain row names (drug names)
          drug_inhibition_df <- data.frame(
            Kinase_Inhibitor = rownames(drug_inhibition_values),
            Inhibition = 100 - as.numeric(drug_inhibition_values[, first_mutation_columns]),
            stringsAsFactors = FALSE
          )
          
          # Remove NAs and select the top 10 drugs based on inhibition (sorted in descending order)
          drug_inhibition_df <- drug_inhibition_df[!is.na(drug_inhibition_df$Inhibition), ]
          drug_inhibition_df <- drug_inhibition_df[order(-drug_inhibition_df$Inhibition), ]
          top_10_drugs <- head(drug_inhibition_df, 10)
          
          inhibiting_drugs_set <- top_10_drugs$Kinase_Inhibitor
          
          # --- Old threshold-based logic (commented out) ---
          # inhibiting_drugs_set <- rownames(mutant_wild_kinases)[apply(mutant_wild_kinases[, first_mutation_columns, drop = FALSE], 1, function(row) any(as.numeric(row) < threshold))]
          
          # If a second kinase is selected, find the common inhibitors for both
          if (!is.null(selected_second_mutation) && selected_second_mutation != "" && selected_second_mutation != "None") {
            second_mutation_columns <- selected_second_mutation
            
            # Get inhibition values for second kinase
            second_drug_inhibition_values <- mutant_wild_kinases[, second_mutation_columns, drop = FALSE]
            
            second_drug_inhibition_df <- data.frame(
              Kinase_Inhibitor = rownames(second_drug_inhibition_values),
              Inhibition = 100 - as.numeric(second_drug_inhibition_values[, second_mutation_columns]),
              stringsAsFactors = FALSE
            )
            
            second_drug_inhibition_df <- second_drug_inhibition_df[!is.na(second_drug_inhibition_df$Inhibition), ]
            second_drug_inhibition_df <- second_drug_inhibition_df[order(-second_drug_inhibition_df$Inhibition), ]
            top_10_second_drugs <- head(second_drug_inhibition_df, 10)
            
            second_inhibiting_drugs <- top_10_second_drugs$Kinase_Inhibitor
            
            # Find drugs common to both kinase inhibition lists
            inhibiting_drugs_set <- intersect(inhibiting_drugs_set, second_inhibiting_drugs)
          }
          
          unique_inhibiting_drugs <- unique(inhibiting_drugs_set)
          gini_selected_drugs <- gini_scores[unique_inhibiting_drugs, "Gini", drop = FALSE]
          
          first_inhibition_col <- paste0("Inhibition (", selected_first_mutation, ") (%)")
          second_inhibition_col <- paste0("Inhibition (", selected_second_mutation, ") (%)")
          
          if (!is.null(selected_second_mutation) && selected_second_mutation != "" && selected_second_mutation != "None") {
            output_data <- data.frame(
              "Kinase Inhibitor" = unique_inhibiting_drugs,
              "Dose (um)" = rep("1.0", length(unique_inhibiting_drugs)),
              first_inhibition_col = 100 - as.numeric(mutant_wild_kinases[unique_inhibiting_drugs, selected_first_mutation]),
              second_inhibition_col = 100 - as.numeric(mutant_wild_kinases[unique_inhibiting_drugs, selected_second_mutation]),
              "Gini Score" = gini_selected_drugs$Gini,
              row.names = NULL
            )
            
            colnames(output_data)[1] <- "Kinase Inhibitor"
            colnames(output_data)[2] <- "Dose (μm)"
            colnames(output_data)[3] <- first_inhibition_col
            colnames(output_data)[4] <- second_inhibition_col
            colnames(output_data)[5] <- "Gini Score"
          } else {
            output_data <- data.frame(
              "Kinase Inhibitor" = unique_inhibiting_drugs,
              "Dose (um)" = rep("1.0", length(unique_inhibiting_drugs)),
              first_inhibition_col = 100 - as.numeric(mutant_wild_kinases[unique_inhibiting_drugs, selected_first_mutation]),
              "Gini Score" = gini_selected_drugs$Gini,
              row.names = NULL
            )
            colnames(output_data)[1] <- "Kinase Inhibitor"
            colnames(output_data)[2] <- "Dose (μm)"
            colnames(output_data)[3] <- first_inhibition_col
            colnames(output_data)[4] <- "Gini Score"
          }
          
          output_data <- output_data %>%
            arrange(desc(!!sym(first_inhibition_col)))
          
          return(output_data)
        } else {
          return(NULL)
        }
      })
      
      # Render the table
      output$drug_table <- renderDataTable({
        filtered_data()
      })
      
      # Conditional UI for download button and column descriptions
      output$conditional_download_and_des <- renderUI({
        if (!is.null(filtered_data()) && nrow(filtered_data()) > 0) {
          tagList(
            downloadButton("download_csv", "Download Table of Results (.csv)"),
            tags$br(),
            tags$hr(),
            tags$div(
              style = "margin-top: 20px;",
              h4(HTML("<b>Column Descriptions</b>"), style = "color:steelblue; font-size: 20px;"),
              tags$ul(
                tags$li(HTML("<b>Kinase Inhibitor:</b> The drug or compound used for inhibition."), style = "font-size: 15px;"),
                tags$li(HTML("<b>Dose (μm):</b> The dose (in micromolar) of the compound used in the study."), style = "font-size: 15px;"),
                tags$li(HTML("<b>Inhibition (kinase) (%):</b> The percentage inhibition of the selected kinase(s) by the drug."), style = "font-size: 15px;"),
                tags$li(HTML("<b>Gini Score:</b> A measure of selectivity for the drug. Higher values indicate greater selectivity."), style = "font-size: 15px;")
              )
            )
          )
        }
      })
      
      # CSV Download
      output$download_csv <- downloadHandler(
        filename = function() {
          all_mutations <- if (!is.null(input$second_mutation) && input$second_mutation != "None" && input$second_mutation != "") {
            paste(input$first_mutation, input$second_mutation, sep = "_")
          } else {
            input$first_mutation
          }
          paste(all_mutations, "_", "results_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(filtered_data(), file, row.names = FALSE)
        }
      )
      
      output$selected_wild_text <- renderUI({
        if (!is.null(input$first_mutation) && input$first_mutation != "") {
          all_mutations <- if (!is.null(input$second_mutation) && input$second_mutation != "None" && input$second_mutation != "") {
            c(input$first_mutation, input$second_mutation)
          } else {
            input$first_mutation
          }
          
          HTML(paste('<b style="font-size: 16px;">Kinase Inhibitors for Selected Kinase(s):</b> ', 
                     paste(all_mutations, collapse = " and "), 
                     '<span style="font-size: 16px;"></span>'))
        } else {
          tagList(
            HTML('<span style="font-size: 16px;">Select kinase(s) of interest</span>'),
            tags$hr()
          )
        }
      })
      
      # Radar Plot Rendering (1 or 2 plots)
      output$radarPlot1 <- renderPlot({
        if (is.null(input$first_mutation) || input$first_mutation == "") return(NULL)
        
        title_text <- paste("Radar Plot for", input$first_mutation)  # Set title dynamically
        
        # Render Radar Plot with Title
        par(mar = c(3, 3, 3, 3))  # Adjust margins for title
        plot.new()
        title(main = title_text, col.main = "navy", font.main = 2, cex.main = 1.8)  # Set title properties
        
        radar_plot_for_kinases(mutant_wild_kinases, input$first_mutation, input$first_mutation)
      })
      
      output$radarPlot2 <- renderPlot({
        if (is.null(input$second_mutation) || input$second_mutation == "None") return(NULL)
        
        title_text <- paste("Radar Plot for", input$second_mutation)  # Set title dynamically
        
        # Render Radar Plot with Title
        par(mar = c(3, 3, 3, 3))  # Adjust margins for title
        plot.new()
        title(main = title_text, col.main = "navy", font.main = 2, cex.main = 1.8)  # Set title properties
        
        radar_plot_for_kinases(mutant_wild_kinases, input$second_mutation, input$second_mutation)
      })
      
      # Download Radar Plots
      output$downloadRadarPlot1 <- downloadHandler(
        filename = function() { paste(input$first_mutation, "_radar_plot.png", sep = "") },
        content = function(file) {
          png(file, width = 800, height = 800, res = 150)
          radar_plot_for_kinases(mutant_wild_kinases, input$first_mutation, input$first_mutation)
          dev.off()
        }
      )
      
      output$downloadRadarPlot2 <- downloadHandler(
        filename = function() { paste(input$second_mutation, "_radar_plot.png", sep = "") },
        content = function(file) {
          if (is.null(input$second_mutation) || input$second_mutation == "None") return(NULL)
          png(file, width = 800, height = 800, res = 150)
          radar_plot_for_kinases(mutant_wild_kinases, input$second_mutation, input$second_mutation)
          dev.off()
        }
      )
      
      # KISS Score Plot for First Kinase
      output$kissPlot1 <- renderPlot({
        if (is.null(input$first_mutation) || input$first_mutation == "") return(NULL)
        
        data1 <- lookup(input$first_mutation)
        if (!is.null(data1) && nrow(data1) > 0) {
          ggplot(data1, aes(x = Kinase_Inhibition, y = KISS)) +
            geom_point(color = "navy", alpha = 1.0, size = 3) +
            labs(
              title = paste("KISS Score vs. Inhibition % for", input$first_mutation),
              x = "Kinase Inhibition (%)",
              y = "KISS Score"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14)
            ) +
            geom_text(
              data = data1 %>% top_n(3, wt = KISS),
              aes(label = Compound),
              vjust = 0.3, hjust = -0.1, fontface = "bold", size = 6
            ) +
            scale_x_reverse()
        }
      })
      
      # KISS Score Plot for Second Kinase
      output$kissPlot2 <- renderPlot({
        if (is.null(input$second_mutation) || input$second_mutation == "None") return(NULL)
        
        data2 <- lookup(input$second_mutation)
        if (!is.null(data2) && nrow(data2) > 0) {
          ggplot(data2, aes(x = Kinase_Inhibition, y = KISS)) +
            geom_point(color = "red", alpha = 1.0, size = 3) +
            labs(
              title = paste("KISS Score vs. Inhibition % for", input$second_mutation),
              x = "Kinase Inhibition (%)",
              y = "KISS Score"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14)
            ) +
            geom_text(
              data = data2 %>% top_n(3, wt = KISS),
              aes(label = Compound),
              vjust = 0.3, hjust = -0.1, fontface = "bold", size = 6
            ) +
            scale_x_reverse()
        }
      })
      
      # Download KISS Plot 1
      output$downloadKissPlot1 <- downloadHandler(
        filename = function() { paste(input$first_mutation, "_KISS_plot.png", sep = "") },
        content = function(file) {
          data1 <- lookup(input$first_mutation)
          if (is.null(data1) || nrow(data1) == 0) return(NULL)
          
          p <- ggplot(data1, aes(x = Kinase_Inhibition, y = KISS)) +
            geom_point(color = "navy", alpha = 1.0, size = 3) +
            labs(
              title = paste("KISS Score vs. Inhibition % for", input$first_mutation),
              x = "Kinase Inhibition (%)",
              y = "KISS Score"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14)
            ) +
            geom_text(
              data = data1 %>% top_n(3, wt = KISS),
              aes(label = Compound),
              vjust = 0.3, hjust = -0.1, fontface = "bold", size = 6
            ) +
            scale_x_reverse()
          
          ggsave(file, plot = p, device = "png", width = 8, height = 6, dpi = 150)
        }
      )
      
      # Download KISS Plot 2
      output$downloadKissPlot2 <- downloadHandler(
        filename = function() { paste(input$second_mutation, "_KISS_plot.png", sep = "") },
        content = function(file) {
          if (is.null(input$second_mutation) || input$second_mutation == "None") return(NULL)
          
          data2 <- lookup(input$second_mutation)
          if (is.null(data2) || nrow(data2) == 0) return(NULL)
          
          p <- ggplot(data2, aes(x = Kinase_Inhibition, y = KISS)) +
            geom_point(color = "red", alpha = 1.0, size = 3) +
            labs(
              title = paste("KISS Score vs. Inhibition % for", input$second_mutation),
              x = "Kinase Inhibition (%)",
              y = "KISS Score"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14)
            ) +
            geom_text(
              data = data2 %>% top_n(3, wt = KISS),
              aes(label = Compound),
              vjust = 0.3, hjust = -0.1, fontface = "bold", size = 6
            ) +
            scale_x_reverse()
          
          ggsave(file, plot = p, device = "png", width = 8, height = 6, dpi = 150)
        }
      )
      
      ### END OF SERVER LOGIC FROM wild_kinases.R
      
    }
  })
  
  # Server logic for Cancer Lineages app (from cancer_app2.R)
  observeEvent(active_app(), {
    if (active_app() == "lineage") {
      
      ### SERVER LOGIC FROM cancer_app2.R
      
      filtered_data <- reactive({
        if (input$kinase == "") return(NULL)
        
        total_counts_original <- all_kinases %>%
          group_by(lineage_1, lineage_2, lineage_3) %>%
          summarise(total_count = n()) %>%
          ungroup()
        
        filtered_kinase_data <- all_kinases %>%
          filter(all_kinases[[input$kinase]] != 0)
        
        lineage_summary_filtered <- filtered_kinase_data %>%
          group_by(lineage_1, lineage_2, lineage_3) %>%
          summarise(count = n()) %>%
          ungroup()
        
        # Join the filtered summary with the total counts from the original data
        lineage_summary <- lineage_summary_filtered %>%
          left_join(total_counts_original, by = c("lineage_1", "lineage_2", "lineage_3")) %>%
          mutate(percentage = (count / total_count) * 100) %>%
          arrange(desc(percentage))
        
        return(lineage_summary)
      })
      
      output$lineageTable <- renderDT({
        data <- filtered_data()
        if (is.null(data)) return(NULL)
        
        column_names <- c("Lineage 1", "Lineage 2", "Lineage 3", "Important Counts", "Total Counts", "Percentage of Important Counts (%)")
        
        datatable(data, colnames = column_names, options = list(pageLength = 15), rownames = FALSE) %>%
          formatRound('percentage', 2) # percentage column to 2 decimal places
      })
      
      # output$radarPlot <- renderPlot({
      #   if (is.null(input$kinase) || input$kinase == "") return(NULL)
      #   
      #   # Find the first partial match for the kinase
      #   matched_kinase <- colnames(mutant_wild_kinases)[grep(tolower(input$kinase), tolower(colnames(mutant_wild_kinases)))][1]
      #   
      #   # If no partial match is found, display the custom message
      #   if (is.na(matched_kinase)) {
      #     plot.new()
      #     text(0.5, 0.5, paste("No FDA Drug Inhibition data available for", input$kinase, "kinase."), cex = 2, col = "red", font = 2)
      #     return(NULL)
      #   }
      #   
      #   # If a partial match is found, proceed with rendering the radar plot
      #   radar_plot_for_kinases(mutant_wild_kinases, matched_kinase, input$kinase)
      # })
      # 
      # output$downloadRadarPlot <- downloadHandler(
      #   filename = function() {
      #     paste(input$kinase, "_radar_plot.png", sep = "")
      #   },
      #   content = function(file) {
      #     tryCatch({
      #       # Set up PNG device with proper dimensions and resolution
      #       png(file, width = 800, height = 800, res = 150)
      #       
      #       # Find the matched kinase for the download handler
      #       matched_kinase <- colnames(mutant_wild_kinases)[grep(tolower(input$kinase), tolower(colnames(mutant_wild_kinases)))][1]
      #       
      #       # Check if matched_kinase is not NA before plotting
      #       if (!is.na(matched_kinase)) {
      #         radar_plot_for_kinases(mutant_wild_kinases, matched_kinase, input$kinase)
      #       } else {
      #         # If no match found, create a blank plot with a message
      #         plot.new()
      #         text(0.5, 0.5, paste("No FDA Drug Inhibition data available for", input$kinase, "kinase."), cex = 2, col = "red", font = 2)
      #       }
      #       
      #       dev.off()  # Close the PNG device to finalize the file
      #     }, error = function(e) {
      #       message("Error in generating the radar plot: ", e$message)
      #     })
      #   }
      # )
      # 
      # # New Plot for KISS Score vs Kinase Inhibition
      # kiss_data <- reactive({
      #   if (input$kinase == "") return(NULL)
      #   lookup(input$kinase)
      # })
      # 
      # output$kissPlot <- renderPlot({
      #   if (is.null(input$kinase) || input$kinase == "") return(NULL)
      #   
      #   data <- kiss_data()
      #   if (is.null(data) || nrow(data) == 0) {
      #     # If no data is found, display the custom message
      #     plot.new()
      #     text(0.5, 0.5, paste("No FDA drug Inhibition data available for", input$kinase, "kinase."), cex = 2, col = "red", font = 2)
      #     return(NULL)
      #   }
      #   
      #   ggplot(data, aes(x = Kinase_Inhibition, y = KISS)) +
      #     geom_point(color = "navy", alpha = 1.0, size = 3) +
      #     labs(
      #       title = paste("KISS Score vs. Inhibition % for", input$kinase),
      #       x = "Kinase Inhibition (%)",
      #       y = "KISS Score"
      #     ) +
      #     theme_minimal() +
      #     theme(
      #       plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      #       axis.title = element_text(size = 16),
      #       axis.text = element_text(size = 14)
      #     ) +
      #     geom_text(
      #       data = data %>% top_n(3, wt = KISS),
      #       aes(label = Compound),
      #       vjust = 0.3, hjust = -0.1, fontface = "bold", size = 6
      #     ) +
      #     scale_x_reverse()
      # })
      # 
      # output$download_kiss_plot <- downloadHandler(
      #   filename = function() {
      #     paste(input$kinase, "_KISS_plot_", Sys.Date(), ".png", sep = "")
      #   },
      #   content = function(file) {
      #     data <- kiss_data()
      #     if (is.null(data)) return(NULL)
      #     
      #     p <- ggplot(data, aes(x = Kinase_Inhibition, y = KISS)) +
      #       geom_point(color = "navy", alpha = 1.0) +
      #       labs(
      #         title = paste("KISS Score vs Inhibition % for", input$kinase),
      #         x = "Kinase Inhibition (%)",
      #         y = "KISS Score"
      #       ) +
      #       theme_minimal() +
      #       theme(
      #         plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      #         axis.title = element_text(size = 16),
      #         axis.text = element_text(size = 14)
      #       ) +
      #       geom_text(
      #         data = data %>% top_n(3, wt = KISS),
      #         aes(label = Compound),
      #         vjust = -1, hjust = 0.5, fontface = "bold", size = 4
      #       ) +
      #       scale_x_reverse()
      #     
      #     ggsave(file, plot = p, device = "png")
      #   }
      # )
      
      # Existing Bar Plot for Kinase Occurrences
      output$barPlot <- renderPlot({
        data <- filtered_data()
        if (is.null(data)) return(NULL)
        
        lineage_counts <- data %>%
          group_by(lineage_1) %>%
          summarise(total_count = sum(count))
        
        ggplot(lineage_counts, aes(x = reorder(lineage_1, total_count), y = total_count)) +
          geom_bar(stat = "identity", fill = "orangered") +
          coord_flip() +
          labs(
            title = paste("Number of Occurrences of", input$kinase, "in Lineage(s) 1"),
            x = "Lineage(s) 1",
            y = "Number of Important Kinases"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Increase plot title size, make it bold, and center it
            axis.title = element_text(size = 16),  # Increase size of axis labels
            axis.text = element_text(size = 14)  # Increase size of axis ticks
          )
      })
      
      # Display selected kinase
      output$selected_kinase_text <- renderUI({
        if (!is.null(input$kinase) && input$kinase != "") {
          HTML(paste('<b>Selected Kinase:</b> ', input$kinase, '<span style="font-size: 22px;"></span>'))
        } else {
          tagList(
            HTML('<span style="font-size: 16px;">Select kinase of interest</span>'),
            tags$hr()
          )
        }
      })
      
      # Download table as CSV
      output$download_table <- downloadHandler(
        filename = function() {
          paste(input$kinase, "_lineage_summary_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          data <- filtered_data()
          write.csv(data, file, row.names = FALSE)
        }
      )
      
      # Download plot as PNG
      output$download_plot <- downloadHandler(
        filename = function() {
          paste(input$kinase, "_bar_plot_", Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
          data <- filtered_data()
          if (is.null(data)) return(NULL)
          
          # Create a bar plot for download
          lineage_counts <- data %>%
            group_by(lineage_1) %>%
            summarise(total_count = sum(count))
          
          p <- ggplot(lineage_counts, aes(x = reorder(lineage_1, total_count), y = total_count)) +
            geom_bar(stat = "identity", fill = "orangered") +
            coord_flip() +
            labs(
              title = paste("Number of Occurrences of", input$kinase, "in Lineage(s) 1"),
              x = "Lineage(s) 1",
              y = "Total Occurrences"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Center plot title and make it bold
              axis.title = element_text(size = 16),  # Increase axis label size
              axis.text = element_text(size = 14)  # Increase axis tick size
            )
          
          # Save plot to file
          ggsave(file, plot = p, device = "png")
        }
      )
      
      ### END OF SERVER LOGIC FROM cancer_app2.R
      
    }
  })
}


options <- list()
if (!interactive()) {
  options$port <- 3838
  options$launch.browser <- FALSE
  options$host <- "0.0.0.0"
} else {
  options$port <- 786
}

shinyApp(ui, server, options = options)