#
# This is the user-interface definition of the PnB Designer web application
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# install.packages("shiny")
# install.packages("purrr")
# install.packages("plyr")
# install.packages("dplyr")
# install.packages("DT")
# install.packages("shinythemes")
# install.packages("shinyBS")
# install.packages("shinyjs")
# install.packages("V8"))
# install.packages("shinybusy")
# install.packages("shinydashboard")

library(shiny)
library("purrr")
library(plyr)
library(dplyr)
library(DT)
library(shinythemes)
library(shinyBS)
library(shinyjs)
library(V8)
library(shinybusy)
library(shinydashboard)

jscode <- "shinyjs.swal = function(params) { swal.apply(this, params); }"

shinyUI( 
  
  fluidPage(
    
    tags$head(
      includeScript("https://cdnjs.cloudflare.com/ajax/libs/sweetalert/1.0.1/sweetalert.min.js"),
      includeCSS("https://cdnjs.cloudflare.com/ajax/libs/sweetalert/1.0.1/sweetalert.min.css")
    ),
    shinyjs::useShinyjs(),
    shinyjs::extendShinyjs(text = jscode),

    theme = shinytheme("yeti"),
                    
                    headerPanel(
                      title = div(img(src='GEML5.png'), "PnB Designer" ),
                      windowTitle = "PnB Designer"
                               ),

                    sidebarPanel(
                      selectInput("Editing", "Please select your editing strategy:",
                        choices = c("Prime editing", "Base editing")),
                      selectInput("Genomes", "Please select the genome of the species, you are working with:",
                        choices = c("","Human (hg38)","Mouse (mm10)","Zebrafish (GRCz11)","Rice (MSU7)","Thale Cress (TAIR9)","Common Grape (IGGP12Xv2)")),
                      selectInput("Mode", "Please select the running mode:",
                        choices = c("","Single Sample Run", "Multi Sample Run")),
                      
                      conditionalPanel(condition="input.Mode=='Single Sample Run'",
                        
                        conditionalPanel(condition="input.Editing == 'Prime editing'",
                          selectInput("SequenceInput", "Please select how you want to insert the editing location:",
                            choices = c("Genomic coordinates", "Sequence input")),
                          textInput("Variant", "Please name your variant:", "HEK3_1"),
                              
                            conditionalPanel(condition="input.SequenceInput=='Sequence input'",
                                             fluidRow(
                                             box(width = 12,
                                             splitLayout(
                                             cellWidths = c("41%","18%","41%"),
                                             HTML(
                                             '<div class="form-group shiny-input-container">
                                              <label for="UpstreamSequence">Upstream Sequence > 75 nt:</label>
                                              <input id="UpstreamSequence" type="text" dir="rtl" class="form-control" value="ATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAA"/>
                                              </div>'
                                                  ),
                                             textInput("Edit2","Edit:","delTCC"),
                                             textInput("DownstreamSequence", "Downstream Sequence > 75 nt:","TTGGGGCCCAGACTGAGCACGTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCT")
                                                        )
                                                )
                                                    )
                                           ),  
                            
                            conditionalPanel(condition="input.SequenceInput=='Genomic coordinates'",
                                             fluidRow(
                                             box(width = 12,
                                             tags$label("Please select the genomic location you want to edit (e.g. chr16 : 107422356):"),
                                             splitLayout(
                                             cellWidths = c("17%","7%","75%"),
                                             textInput("Chromosome", "", "chr9"),
                                             disabled(textInput("placeholder","",":")),
                                             textInput("Position", "",107422356)
                                                        )
                                                )
                                                    ),
                                             selectInput("Gene Orientation", "Please select the orientation of your target gene:",
                                              choices = c("+", "-")),
                                             textInput("Edit", HTML("Please select the edit you want to install (e.g. G>A, insA, delT)"), "delT"),
                                             textInput("Mutation", HTML("Please select the mutation you want to correct (e.g. C>T, insTT, delTCT)"), "")
                                             
                                            )
                                        ),
                        
                      conditionalPanel(condition = "input.Editing == 'Base editing'",
                        textInput("Variant2", "Please name your variant:", "NM_000517.6(HBA2):c.99G>A"),
                        selectInput("SequenceInput2", "Please select how you want to insert the editing location:",
                          choices = c("Genomic coordinates", "Sequence input")),
                            
                           conditionalPanel(condition="input.SequenceInput2=='Sequence input'",
                                            fluidRow(
                                            box(width = 12,
                                            splitLayout(
                                            cellWidths = c("41%","18%","41%"),
                                            HTML(
                                            '<div class="form-group shiny-input-container">
                                             <label for="UpstreamSequence">Upstream Sequence > 25 nt:</label>
                                             <input id="UpstreamSequence2" type="text" dir="rtl" class="form-control" value="ATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAA"/>
                                             </div>'
                                                 ),
                                            textInput("Edit3","Edit","T>C"),
                                            textInput("DownstreamSequence2", "Downstream Sequence > 25 nt:","CCTTGGGGCCCAGACTGAGCACGTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCT")
                                                        )
                                               )
                                                   )
                                          ),
                        
                           conditionalPanel(condition="input.SequenceInput2=='Genomic coordinates'",
                                            fluidRow(
                                            box(width = 12,
                                            tags$label("Please select the genomic location you want to edit (e.g. chr16 : 107422356):"),
                                            splitLayout(
                                            cellWidths = c("17%","7%","75%"),
                                            textInput("Chromosome2", "", "chr16"),
                                            disabled(textInput("placeholder","",":")),
                                            textInput("Location","",173128)
                                                        )
                                                )
                                                    ),
                                            selectInput("Gene Orientation2", "Please select the orientation of your target gene:",
                                              choices = c("+", "-")),
                                            selectInput("SNP", "Please select your single nucleotide change you want to correct:",
                                              choices = c("","G>A","C>T")))),
                    conditionalPanel(condition = "input.Editing == 'Prime editing'",
                      numericInput("PBS", "Please select the PBS length:", 13, min = 1, max = 17, step = 1), 
                      numericInput("RT", "Please select the RTT length:", 13, min = 1, max = 80 , step = 1)
                                                          
                                     )
                            ),
                    
            conditionalPanel(condition="input.Mode=='Multi Sample Run'",
                             fileInput("file1", "Choose CSV File (max. 1000 samples)",
                                       accept = c(
                                       "text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")),
                                       tags$hr(),
                             conditionalPanel(condition = "input.Editing == 'Prime editing'",
                              checkboxInput("expert", "Show all possible Oligos", FALSE ),
                                             ),
                           ),
  
                           tags$head(tags$script(src = "message-handler.js")),
                           useShinyjs(),
                           bsButton("search", "Search", type = "action"),
                           bsButton("reset", "Reset", type = "action"),
                           tags$head(tags$script(src = "message-handler.js")),
                           fluidRow(
                            column(1, offset = 9,
                            tags$a(href='Walk-through.pdf', target='blank', 'Instructions', download = 'Walk-through.pdf'),
                                  ),
                                   ),
          conditionalPanel(condition="input.Mode=='Multi Sample Run' && input.Editing == 'Base editing'",
                           hr(),
                           tags$a(href='Base_template.csv', target='blank', 'Base Editing template file', download = 'Base_template.csv'),
                          ),
          
          conditionalPanel(condition="input.Mode=='Multi Sample Run'&& input.Editing == 'Prime editing'",
                         
                          tags$a(href='Prime_template.csv', target='blank', 'Prime Editing template file', download = 'Prime_template.csv')
        
                          ),
          hr(),
          absolutePanel(
            left = "10%",
            sidebarLayout("Questions regarding PnB Designer?" ,tags$a(href="mailto:pngdesigner@gmail.com", target='blank', 'Contact us', link = 'pngdesigner@gmail.com'),)
          ),
          br()
            ),
     
           mainPanel(
             
             
             
                     add_busy_spinner(spin = "atom" , position = "top-right", margins = c(300, 500)),
                     tags$head(HTML("<script type='text/javascript' src='path/to/sweetAlert2.js'></script>")),
                     textOutput("myText"),
                     br(),
                     tags$head(tags$style("#myText{color: black;
                               font-size: 20px;
                               font-style: italic;
                               }"
                                         )
                              ), 
          
                     DT::dataTableOutput("mytable"),
                     textOutput("safeError"),
                     fluidRow(
                     br(),
                     conditionalPanel(condition = "input.Mode == 'Single Sample Run'",
                                      div(dataTableOutput("mytable2"), style = "font-size:40%"),
                                      br(),
                                     ),
                     column(1, offset = 9,
                     uiOutput("download")
                           ),
                            ),
                     br()
                     

                  ),

    #absolutePanel(
    #  right = "1%",
    #  top = "2%",
    #  textOutput("count")
    #)

    
          )
    )
  
