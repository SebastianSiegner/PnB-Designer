#
# This is the user-interface definition of the PnB Designer web application
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

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
#install.packages("shinyWidgets")

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
library(shinyWidgets)


#library(rsconnect)
#setwd("C:/Users/wien1/Desktop/R-wd/ShinyB&PEditing0.99 - With Reviews///")
#options(repos = BiocManager::repositories())

#rsconnect::setAccountInfo(name='cornlab',
#                      token='CFFC0187DDBF6EBD4449F9E6AE0106E7',
#                         secret='9hX7TZQHRCs4IuVVOqfwImGPXCdnsmffA4JTJ+nC')

#rsconnect::deployApp(appName="ShinyBPEditing0999")

jscode <- "shinyjs.swal = function(params) { swal.apply(this, params); }"

shinyUI( 
 
  fluidPage(
    
    tags$style("
              #content'{
   width: 600px;
   margin-left: auto;
   margin-right: auto;
}'"),

    tags$head(
      includeScript("https://cdnjs.cloudflare.com/ajax/libs/sweetalert/1.0.1/sweetalert.min.js"),
      includeCSS("https://cdnjs.cloudflare.com/ajax/libs/sweetalert/1.0.1/sweetalert.min.css")
    ),
    
    shinyjs::useShinyjs(),
    #shinyjs::extendShinyjs(text = jscode),


    theme = shinytheme("yeti"),
                    
                    headerPanel(
                      title = div(img(src='GEML5.png'), "PnB Designer" ),
                      windowTitle = "PnB Designer"
                               ),

    sidebarLayout(
    
                    sidebarPanel(
                      tags$style(type='text/css',".selectize-input {font-size: 14px; line-height: 28px;}"),
                      selectInput("Editing",label = div(style = "font-size:14px", "Please select your editing strategy:"), 
                        choices = c("Prime editing", "Base editing")),
                      selectInput("Genomes", label = div(style = "font-size:14px", "Please select the genome of the species, you are working with:"),
                        choices = c("","Human (hg38)","Mouse (mm10)","Zebrafish (GRCz11)","Rice (MSU7)","Thale Cress (TAIR9)","Common Grape (IGGP12Xv2)","None of the above")),
                      selectInput("Mode",label = div(style = "font-size:14px", "Please select the running mode:"),
                        choices = c("","Single Sample Run", "Multi Sample Run")),
                      
                      conditionalPanel(condition="input.Mode=='Single Sample Run'",
                        
                        conditionalPanel(condition="input.Editing == 'Prime editing'",
                          selectInput("SequenceInput", label = div(style = "font-size:14px", "Please select how you want to insert the editing location:"),
                            choices = c("Genomic coordinates", "Sequence input")),
                          textInput("Variant", label = div(style = "font-size:14px", "Please name your variant:"), "HEK3_1"),
                            conditionalPanel(condition="input.SequenceInput=='Sequence input'",
                                             fluidRow(
                                             box(width = 12,
                                             splitLayout(
                                             cellWidths = c("41%","18%","41%"),
                                             HTML(
                                             '<div class="form-group shiny-input-container">
                                              <label for= "UpstreamSequence"> Upstream Seq. > 150 nt:</label>
                                              <input id="UpstreamSequence" type="text" dir="rtl" class="form-control" value="GTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAA"/>
                                              </div>'
                                                  ),
                                             textInput("Edit2",label = div(style = "font-size:14px","Edit:"),"delTCC"),
                                             textInput("DownstreamSequence", "Downstream Seq. > 150 nt:","TTGGGGCCCAGACTGAGCACGTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCA")
                                                        )
                                                )
                                                    )
                                           ),  
                            
                            conditionalPanel(condition="input.SequenceInput=='Genomic coordinates'",
                                             fluidRow(
                                             box(width = 12,
                                             tags$label(div(style = "font-size:14px","Please select the genomic location you want to edit:")),
                                             splitLayout(
                                             cellWidths = c("17%","7%","75%"),
                                             textInput("Chromosome", "", "chr9"),
                                             disabled(textInput("placeholder","",":")),
                                             textInput("Position", "",107422356)
                                                        )
                                                )
                                                    ),
                                             selectInput("Gene Orientation",div(style = "font-size:14px", "Please select the orientation of your target gene:"),
                                              choices = c("+", "-")),
                                             switchInput(inputId = "SwitchButton", value = TRUE, onLabel = "Install an edit",
                                                         offLabel = "Correct mutation", ),
                                                         
                                                        conditionalPanel(condition = "input.SwitchButton == 0",
                                                        textInput("Mutation", label = div(style = "font-size:14px","Please select the mutation you want to correct (e.g. C>T, insTT):"), "")),
                                             
                                                        conditionalPanel(condition = "input.SwitchButton == 1",
                                                        textInput("Edit", div(style = "font-size:14px","Please select an edit you want to install (e.g. G>A, insA, delT):"), "delT"))
                                             
                                            )
                                        ),
                        
                      conditionalPanel(condition = "input.Editing == 'Base editing'",
                        textInput("Variant2", label = div(style = "font-size:14px","Please name your variant:"), "NM_000517.6(HBA2):c.99G>A"),
                        selectInput("SequenceInput2", label = div(style = "font-size:14px","Please select how you want to insert the editing location:"),
                          choices = c("Genomic coordinates", "Sequence input")),
                            
                           conditionalPanel(condition="input.SequenceInput2=='Sequence input'",
                                            fluidRow(
                                            box(width = 12,
                                            splitLayout(
                                            cellWidths = c("41%","18%","41%"),
                                            HTML(
                                            '<div class="form-group shiny-input-container">
                                             <label for="UpstreamSequence">Upstream Sequence > 25 nt:</label>
                                             <input id="UpstreamSequence2" type="text" dir="rtl" class="form-control" value="AACCCCACCCCTCACTCTGCTTCTCCCCGCAGGAT"/>
                                             </div>'
                                                 ),
                                            textInput("Edit3","Edit","A>G"),
                                            textInput("DownstreamSequence2", "Downstream Sequence > 25 nt:","TTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTC")
                                                        )
                                               )
                                                   )
                                          ),
                        
                           conditionalPanel(condition="input.SequenceInput2=='Genomic coordinates'",
                                            fluidRow(
                                            box(width = 12,
                                            tags$label(div(style = "font-size:14px","Please select the genomic location you want to edit:")),
                                            splitLayout(
                                            cellWidths = c("17%","7%","75%"),
                                            textInput("Chromosome2", "","chr16"),
                                            disabled(textInput("placeholder","",":")),
                                            textInput("Location","",173128)
                                                        )
                                                )
                                                    ),
                                            selectInput("Gene Orientation2", label = div(style = "font-size:14px","Please select the orientation of your target gene:"),
                                              choices = c("+", "-")),
                                            selectInput("SNP", label = div(style = "font-size:14px","Please select your single nucleotide change you want to correct:"),
                                              choices = c("","G>A","C>T","T>C","A>G")))),
                    conditionalPanel(condition = "input.Editing == 'Prime editing'",
                      numericInput("PBS", label = div(style = "font-size:14px","Please select the PBS length:"), 13, min = 1, max = 17, step = 1), 
                      numericInput("RT", label = div(style = "font-size:14px","Please select the RTT length:"), 13, min = 1, max = 80 , step = 1),
                      #checkboxInput("PE3", "Please tick, if you want to use the PE3 system with an additional nicking guide", FALSE),                  
                                     )
                            ),
                    
            conditionalPanel(condition="input.Mode=='Multi Sample Run'",
                             fileInput("file1", label = div(style = "font-size:11px","Choose CSV File (max. 1000 samples)"),
                                       accept = c(
                                       "text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")),
                                       tags$hr(),
                             conditionalPanel(condition = "input.Editing == 'Prime editing'",
                              checkboxInput("expert", "Show all possible pegRNAs", FALSE ),
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
            sidebarLayout("Questions regarding PnB Designer?" ,tags$a(href="mailto:pnbdesigner@gmail.com", target='blank', 'Contact us', link = 'pngdesigner@gmail.com'),)
          ),
          br()
          
                    ,width = 4),
          
           mainPanel(

             tags$head(HTML("<script type='text/javascript' src='path/to/sweetAlert2.js'></script>")),

                     add_busy_spinner(spin = "atom" , position = "top-right", margins = c(300, 500)),
                      
                     
                     textOutput("myText"),
                     br(),
                     tags$head(tags$style("#myText{color: black;
                               font-size: 20px;
                               font-style: italic;
                               }"
                                         )
                              ), 
                     textOutput("safeError"),
                     br(),
             HTML("<div style =' overflow: auto; scrollbar-width: none; zoom: 80%;' >"),
             DT::dataTableOutput("mytable"),
             HTML("</div>"),
                          
             
         conditionalPanel(condition="input.Mode=='Single Sample Run'&& input.Editing == 'Prime editing'",
                          
                          br(),
                          
                          HTML("<div style =' overflow: auto; scrollbar-width: none; zoom: 80%;' >"),
                          DT::dataTableOutput("nickingguides"),
                          HTML("</div>"),
                          
                          br(),
                          
                          HTML("<div style =' zoom: 80%;' >"),
                          column(4, offset = 2, uiOutput("download")),
                          HTML("</div>"),
                          
                          HTML("<div style =' zoom: 80%;' >"),
                          column(1, offset = 1, uiOutput("download2")),
                          HTML("</div>"),
         ),
         
         conditionalPanel(condition="input.Mode=='Multi Sample Run'& input.Editing == 'Prime editing'",
                          
          br(),
          HTML("<div style =' zoom: 80%;' >"),
          column(4, offset = 1, uiOutput("download3")),
          HTML("</div>"),
          
          HTML("<div style ='  zoom: 80%;' >"),
          column(1, offset = 2, uiOutput("download4")),
          HTML("</div>"),
                          
         ),
         
         conditionalPanel(condition="input.Editing == 'Base editing'",
                          
                          br(),
                          
                          HTML("<div style =' zoom: 80%;' >"),
                          column(1, offset = 2, uiOutput("download5")),
                          HTML("</div>"),
                          
                          br(),
                          
         ),
         
           br()
         
         ,width = 8),
  ),

          )
    
)
  
