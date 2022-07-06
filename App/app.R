# Shiny App Code#
#################

#### SETUP #### 
# source model code
source("contact_tracing_v5.R")

# libraries
library(shinythemes)
library(shinyjs)
library(plotly)
library(RColorBrewer)
library(tidyverse)

# housekeeping
# radio button values
choiceNames= c("Fraction of high contact group cases detected in community", 
               "Fraction of contacts successfully traced",
               "Isolation and quarantine efficacy (first generation)", 
               "Isolation and quarantine efficacy (subsequent generations)")
choiceValues = c("HiMSM_prob.det", "contact_trace_prob", "adh", "adh2")

#### UI #### 
ui <- fluidPage(
  # THEME
  theme=shinytheme("simplex"),
  
  # JS - for reset function
  useShinyjs(),
  
  # TITLE
  titlePanel("COVID-19 Contact Tracing Model"),
  
  # SIDEBAR
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        
        # MAIN PARAMETERS
        tabPanel("Main", fluid=TRUE,
                 
                 h4("Community testing"),
                 h5("These sliders describe detection of infection in the community (among individuals who are not traced contacts)."),
                 sliderInput("HiMSM_prob.det", "Fraction of cases detected (high contact)", 
                             min=0, max=1, value=.5, step = 0.1),                               
                 sliderInput("LoMSM_prob.det", "Fraction of cases detected (low contact)", 
                             min=0, max=1, value=.1, step = 0.1),  
                 
                 h4("Contact tracing"),
                 sliderInput("contact_trace_prob", "Fraction of contacts successfully traced", min=0, max=1, value=0.25, step = 0.05), 
                 
                 h5("Quarantine and treatment efficacy is the fraction of transmission prevented in traced contacts."),
                 sliderInput("adh", "First generation of traced contacts", 
                             min=0.5, max=.95, value=.75, step = 0.05),
                 sliderInput("adh2", "Second generation of traced contacts", 
                             min=0.5, max=.95, value=.75, step = 0.05), 
                 
                 h4("Epidemiology"),
                 sliderInput("SAR", "Secondary Attach Rate", 
                             min=0.1, max=0.5, value=0.2, step = 0.05),  
                 sliderInput("HiMSM_contacts","Contacts, High Contact group",
                             min=3, max=30, value=5, step =1),    
                 sliderInput("LoMSM_contacts","Contacts, Low Contact group", 
                             min=0, max=5, value=1.5, step =0.5),  
                 sliderInput("LoMSM_prob","Fraction of Populion in Low Contact group", 
                             min=0, max=1, value=0.6, step =0.1), 
                 
                 radioButtons("xaxis", "X-axis variable:", 
                              choiceNames= choiceNames, 
                              choiceValues = choiceValues,
                              selected = NULL, inline = FALSE, width = NULL),
                 
                 # save and download
                 actionButton("restore_all", "Restore original inputs")#,
                 #downloadButton(outputId = "download_Inputs", 
                 #               label = 'Download inputs',
                 #               class= "mybutton"),
                 #downloadButton(outputId = "download_Data", 
                 #               label = 'Download estimates',
                 #               class= "mybutton")
                 
        ),
        
        # ADVANCED PARAMETERS
        tabPanel("Advanced", fluid = TRUE,
                 
                # h4("Detection"),
                # 
                # # program
                # sliderInput("rel_trans", "Relative risk of transmission among detected cases (vs. undetected)", min=0, max=1, value=0.5, step = 0.01), 
                # sliderInput("test_uptake", "Fraction of eligible contacts tested", min=0, max=1, value=0.9, step = 0.01),
                # 
                # # epidemiology
                # h4("Pre-symptomatic"),
                # sliderInput("P_RR", "Relative risk of transmission (vs. symptomatic)", min=0, max=1.5, value=1, step = 0.01),
                # sliderInput("P_dur", "Duration of pre-symptomatic transmission", min=0, max=3, value=1.5, step = 0.01),  
                # 
                # h4("Symptomatic"),
                # sliderInput("S_dur", "Duration of symptomatic transmission", min=0, max=10, value=4, step = 0.01),   
                # 
                # h4("Asymptomatic"),
                # sliderInput("A_RR", "Relative risk of transmission (vs. symptomatic)", min=0, max=1.5, value=.7, step = 0.01),                               
                # sliderInput("A_dur", "Duration of asymptomatic transmission", min=0, max=10, value=5.5, step = 0.01),   
                # 
                # sliderInput("A_prob", "Fraction of cases that are asymptomatic", min=0, max=1, value=0.4, step = 0.01),
                # 
                # h4("Strategy"),
                # HTML("<strong>'Contact tracing only'</strong> compares contact tracing vs no contact tracing at the same level of testing.  <strong>'Testing scale-up + contact tracing</strong>'
                #  compares contact tracing to no contact tracing with a user-selected baseline level of testing.  This is sensitive to the 
                #  to assumptions about both the baseline detection fraction and the behavioral benefits of testing (RR of transmission among detected cases)."),
                # 
                # # comparator             
                # radioButtons("comparator", "", choices = c("Contact tracing only", "Testing scale-up + contact tracing"),
                #              inline = FALSE, width = NULL, selected = "Contact tracing only",
                #              choiceValues = NULL),
                # conditionalPanel(
                #   "input.comparator == 'Testing scale-up + contact tracing'",
                #   sliderInput("baseline_S_prob.det", "Baseline fraction of symptomatic cases detected", min = 0, max = 1, value = .2, step = 0.01),
                #   sliderInput("baseline_A_prob.det", "Baseline fraction of asymptomatic cases detected", min = 0, max = 1, value = 0, step = 0.01)
                #   
                # ),
                 
      )
        
      ),
      width=3),
    
    # RESULTS          
    mainPanel(tabsetPanel(type = "tabs",
                          
                          # model output                
                          tabPanel("Model results", 
                                   h4(""),
                                   htmlOutput("txt"),
                                   
                                   # output plots
                                   plotlyOutput("plot1", height = "450", width = "1300"),
                                   plotlyOutput("plot2", height = "450", width = "1300"),
                                   htmlOutput("txt2"),
                                   
                                   includeMarkdown("content/start.md"),
                          ),
                          
                          # upload and view inputs
                          tabPanel("Upload inputs" #,
                          #         h4("Input data"),
                          #         fileInput("file", "Choose CSV File",
                          #                   multiple = FALSE,
                          #                   accept = c("text/csv",
                          #                              "text/comma-separated-values,text/plain",
                          #                              ".csv")),
                          #         tableOutput("tbl2")
                          ),
                          
                          # documentation
                          tabPanel("Documentation"#,
                          #         includeMarkdown("content/model.md")
                          )
    ))),
  hr()#,
  #includeMarkdown("content/footer.md")
)

#### SERVER #### 
server <- function(input, output, session) {
  
### INPUTS ###
  
  # LOAD INPUTS
#  observeEvent(input$file, {
#    
#    inFile = input$file
#    if(!is.null(inFile)){
#      # Load inputs
#      uploaded_inputs <- read.csv(inFile$datapath)
#      # Update each input
#      for(i in 1:nrow(uploaded_inputs)){
#        if(!grepl("testing|comparator|xaxis", uploaded_inputs$inputId[i])) { 
#          updateSliderInput(session,
#                            inputId = uploaded_inputs$inputId[i],
#                            value = uploaded_inputs$value[i])
#        }else{
#          updateRadioButtons(session,inputId = uploaded_inputs$inputId[i],
#                             selected = uploaded_inputs$value[i])
#        }
#      }
#    }
#    
#  })
#  
#  # MAKE INPUT TABLE
#  inputs = reactive({
#    inFile = input$file
#    if (is.null(inFile)) {
#      return(NULL)
#    } else{ return(read.csv(inFile$datapath))}
#  })
#  output$tbl2 = renderTable({ inputs() })

### MODEL ###
  # RESTORE BUTTON
  observeEvent(input$restore_all, {
    vars = c(SAR, # secondary attack rate
             HiMSM_contacts, # avg. (daily) contacts, High Contact group 
             LoMSM_contacts, # avg. (daily) contacts, Low Contact group 
             duration, # duration of infectiousness
             HiMSM_prob.det, # detection probability, High Contact
             LoMSM_prob.det, # detection probability, Low Contact
             LoMSM_prob, # fraction of infections that are Low Contact. Moot?
             contact_trace_prob, # probability of contact tracing 
             comparator, # character string, assigns a scenario
             baseline_HiMSM_prob.det, # prob detection w/o tracing (Hi)
             baseline_LoMSM_prob.det, # prob detection w/o tracing (Lo)
             test_uptake, # probability tested (if traced? without traced?)
             adh, # adherence / % redux in transmission
             adh2, # adherence / % redux in transmission
             vax, # fraction of MSM population with vax
             rel_trans, # relative number of secondary infections (detected 
             # compared to undetected)
             xaxis)
    
    for(i in 1:length(vars)) reset(vars[i])
  })
  
  # RUN CONTACT TRACING MODEL
  out <- reactive({
    
    get_R(SAR                      = input$SAR, 
          HiMSM_contacts           = input$HiMSM_contacts, 
          LoMSM_contacts           = input$LoMSM_contacts, 
          duration                 = 21,
          HiMSM_prob.det           = input$HiMSM_prob.det,
          LoMSM_prob.det           = input$LoMSM_prob.det, 
          LoMSM_prob               = input$LoMSM_prob,
          contact_trace_prob       = input$contact_trace_prob, 
          comparator               = "Contact tracing only",
          baseline_HiMSM_prob.det  = 0.1, 
          baseline_LoMSM_prob.det  = 0.05, 
          test_uptake              = 0.5, 
          adh                      = input$adh, 
          adh2                     = input$adh2,
          vax                      = input$vax,
          rel_trans                = 0.5, 
          xaxis                    = input$xaxis)
  
    })
  
  # MAKE TABLE
  output$tbl = renderTable({ out() })
  
  # MAKE PLOTS 
  plots = reactive({ 
    make_plots(out(), xaxis = choiceNames[which(choiceValues==input$xaxis)])
    })
  
  # MAKE PLOT TOP ROW
  output$plot1 <- renderPlotly({
    plots()[[1]]
  })
  
  # MAKE PLOT BOTTOM ROW
  output$plot2 <- renderPlotly({
    subplot(style(plots()[[2]], showlegend = FALSE),
            plots()[[3]])
  })
  
  # MAKE TEXT
  output$txt = renderText({
    paste("<font size='4'><strong>What can contact tracing achieve?</strong> 
    This model shows how contact tracing impacts the effective reproduction 
    number, R(t).")
    })
  
  output$txt2 = renderText({
    paste("<font size = '4'><strong>Model approach</strong>: if the epidemic 
          started with the basic reproductive number <strong>R<sub>0</sub>=", 
          input$R0, "</strong>, and physical distancing measures achieved 
          <strong>R(t)=", input$Rt, "</strong>, what happens when we add contact 
          tracing?You can adjust parameters with the sliders on the left.  For 
          more details, see documentation tab.</font>  ", sep = "")
  })


  
  ### OUTPUTS ###
  # DOWNLOAD INPUTS
  #output$download_Inputs <- downloadHandler(
  #  filename = function() {
  #    paste("inputs_", Sys.Date(), ".csv", sep="")
  #  },
  #  content = function(file) {
  #    # Define inputs to save
  #    inputs_to_save <- c("P_RR", "P_dur", "S_dur",
  #                                "A_RR", "A_dur",
  #                                "S_prob.det", "A_prob.det", "A_prob", "contact_trace_prob", 
  #                                 "R0", "Rt", "comparator",
  #                                "baseline_S_prob.det", "baseline_A_prob.det", "test_uptake", "adh", "adh2", "rel_trans", "xaxis")
  #
  #    # Declare inputs
  #    inputs <- NULL
  #    # Append all inputs before saving to folder
  #    for(input.i in inputs_to_save){
  #      inputs <- append(inputs, input[[input.i]])
  #    }
  #    # Inputs data.frame
  #    inputs_data_frame <- data.frame(inputId = inputs_to_save, value = inputs)
  #    # Save Inputs
  #    write.csv(inputs_data_frame, file, row.names = FALSE)
  #  }  
  #)
  #
  ## DOWNLOAD ESTIMATES
  #output$download_Data <- downloadHandler(
  #  filename = function() {
  #    paste("data_", Sys.Date(), ".csv", sep="")
  #  },
  #  content = function(file) {
  #    write.csv(out() %>% select(-point), file)
  #  }  
  #)
#
}
#
##### APP CALL ####
shinyApp(ui = ui, server = server)
