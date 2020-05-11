

#### R.server script for PnB Designer
## DataTableOutput is produced for displaying the resulting oligos



# install.packages("devtools")
# install.packages("crayon")
# install.packages("shiny")
# devtools::install_github("timelyportfolio/sweetalertR")
# install.packages("BSgenome")

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

#BiocManager::install("BSgenome.Osativa.MSU.MSU7")

#BiocManager::install("BSgenome.Drerio.UCSC.danRer11")

#BiocManager::install("BSgenome.Vvinifera.URGI.IGGP12Xv2")

#BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")


library(devtools)
library(crayon)
library(shiny)
library(sweetalertR)
library(shinyjs)
library(V8)
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Osativa.MSU.MSU7")
library("BSgenome.Drerio.UCSC.danRer11")
library("BSgenome.Vvinifera.URGI.IGGP12Xv2")
library("BSgenome.Athaliana.TAIR.TAIR9")
library("BSgenome")


# Defines the possible input file size (At the moment set to 8 MB)

options(shiny.maxRequestSize=8*1024^2) 


# Define server logic required to produce Output Table for Oligos:

shinyServer(function(input, output, session) {

# Reset Button implementation in the server logic. 
# Important parameters are set to NULL before start of another search round.
  
  observeEvent(input$reset,{
    
    output$mytable <- NULL
    output$myText <- NULL
    output$README <- NULL
    output$mytable2 <- NULL
    c <- 1
    d <- 1
    x <- 0
    output$mytable <- NULL
    
    ## Download Button is hidden after pressing the reset button
    
    updateCheckboxInput(session, "Download", value = FALSE)
    updateButton(session, "search", disabled = FALSE)
    shinyjs::hideElement(id = "download")
    
  })

## ObserveEvent function to start a search round throught the genome after search button is pressed:
  
  observeEvent(input$search, {
    
   
    ## Target genome is chosen by the selected genome in the User Interface
    
    if((as.character(input$Genomes) == "Human (hg38)")==TRUE)
    {
      Target_Genome <- BSgenome.Hsapiens.UCSC.hg38
      
    }
    
    else if((as.character(input$Genomes) == "Mouse (mm10)")==TRUE)
    {
      Target_Genome <- BSgenome.Mmusculus.UCSC.mm10
      
    }
    else if((as.character(input$Genomes) == "Zebrafish (GRCz11)")==TRUE)
    {
      Target_Genome <- BSgenome.Drerio.UCSC.danRer11
      
    }
    else if((as.character(input$Genomes) == "Rice (MSU7)") ==TRUE)
    {
      
      Target_Genome <- BSgenome.Osativa.MSU.MSU7
      
    }
    else if((as.character(input$Genomes) == "Common Grape (IGGP12Xv2)")==TRUE)
    {
      Target_Genome <- BSgenome.Vvinifera.URGI.IGGP12Xv2
      
    }
    else if((as.character(input$Genomes) == "Thale Cress (TAIR9)")==TRUE)
    {
      Target_Genome <- BSgenome.Athaliana.TAIR.TAIR9
      
    }

#-----------------------------------------------------------------------------------------------------------------#    
#-----------------------------------------------------------------------------------------------------------------#  
    
    ## Start of a Single Sample Run, if this option is selected in the User Interface:

      if(input$Mode == "Single Sample Run")
      {

        
    ## Correction of the different Chromosome notation in the Rice genome (Chr6 instead of chr6 in all other genomes)
  
          if(input$Genomes != "Rice (MSU7)")
          {
            
            Chromosome <- input$Chromosome
            Chromosome2 <- input$Chromosome2
            
            
          }
          else if(input$Genomes == "Rice (MSU7)")
          {
            
            Chromosome <- paste("chr",substring(input$Chromosome, first = 4),sep = "")
            
            Chromosome2 <- paste("chr",substring(input$Chromosome2, first = 4),sep = "")
            
          }
        
     
#-----------------------------------------------------------------------------------------------------------------#
        
    # Base editing script:
    # If Base editing is selected in the Editing strategy this script will be run: 

      if(input$Editing == "Base editing")
        
      {
        
        
    # DownloadButton is implemented in the UI, producing a .csv ouput file:
        
        output$downloadData <- downloadHandler(
          filename = function() {
            paste('Base editing guides-', Sys.Date(), '.csv', sep='')
          },
          content = function(con) {
            write.csv(data.frame("Variant" = input$Variant, "Protospacer" = guides[,2], "EditPos." = guides[,3], "PAM" = guides[,4], "Base Editor" = guides[,5]), con, sep = ",", row.names=FALSE)
          }
        )
        
    # a, a variable for later quality control is defined: 
        
        a <- 1
        
    # create 'guides', which the later created guides will be stored to:
  
        guides <- NULL
  
    # create 'Download', which the later created protospacer sequence without html tags will be stored to:
        
        Download <- NULL
        
#-------------------------------------------------#
        
        ### Get sequence for Plus Strand using the Bioconductor package
        ### Get Sequence for Minus Strand by using reverseComplement function (bioconductor)

        if(input$SequenceInput2 == "Genomic coordinates")
        {
        
        Seq_plus <- getSeq(Target_Genome, GRanges(Chromosome2,ranges = IRanges(start = as.numeric(input$Location)-25, width = 51),))
        
        Seq_minus <- reverseComplement(Seq_plus)
        
        # Convert the DnaStringSets into character strings:
        
        Seq_plus_char <- as.character(Seq_plus)
        
        Seq_minus_char <- as.character(Seq_minus)
        
        print(Seq_minus_char)
        
        }
        if (input$SequenceInput2 == "Sequence input")
        {
          
          RIGHT = function(x,n){
            substring(x,nchar(x)-n+1)
          }
          
          
          Seq_plus <- DNAStringSet(paste(RIGHT(input$UpstreamSequence,25),substr(input$Edit3,1,1),substring(input$DownstreamSequence,0,25), sep = ""))
          
          Seq_minus <- reverseComplement(Seq_plus)
          
          # Convert the DnaStringSets into character strings:
          
          Seq_plus_char <- as.character(Seq_plus)
          
          Seq_minus_char <- as.character(Seq_minus)
          
        }
        
       # Create additional character strings, in which the mutation found in the patient will be included:
       # '-NGG' is for ABEmax with a NGG PAM; '-NGA' for NG-ABEmax/VRQR-ABEmax with a NGA PAM; '-NGCG' for NG-ABEmax;
       # '-NNGRRT' is for Sa-ABEmax; '-NNHRRT' is for SaKKH-ABEmax
        
        Seq_minus_char_NGG <- Seq_minus_char
        
        Seq_minus_char_NGA <- Seq_minus_char
        
        Seq_minus_char_NGCG <- Seq_minus_char
        
        Seq_minus_char_NNGRRT <- Seq_minus_char
        
        Seq_minus_char_NNHRRT <- Seq_minus_char
        
        Seq_minus_char_NGC <- Seq_minus_char
        
        Seq_minus_char_NYN <- Seq_minus_char
        
        
        Seq_plus_char_NGG <- Seq_plus_char
        
        Seq_plus_char_NGA <- Seq_plus_char
        
        Seq_plus_char_NGCG <- Seq_plus_char
        
        Seq_plus_char_NNGRRT <- Seq_plus_char
        
        Seq_plus_char_NNHRRT <- Seq_plus_char
        
        Seq_plus_char_NGC <- Seq_plus_char
        
        Seq_plus_char_NYN <- Seq_plus_char
        
        
        # Guide search For ABEmax (NGG PAM) depending on input'-SNP' and '-Gene Orientation':
        # For Base editing guides on the Minus Strand:

        if(input$SNP == "G>A" && input$`Gene Orientation2` == "-" || input$SNP == "C>T" && input$`Gene Orientation2` == "+" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "T>C" ){
          
        # Test wildtype sequence if input'-SNP' was assigned right and if so includes patient mutation into character string:
        # If not, see code line 216
          
         if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
         {
            substr(Seq_minus_char_NGG, start = 26, stop = 26) <- "A"
          
        # patient character string is now searched for suitable PAM sequences:

          # Perfect edible window: For ABEmax (Editing Position: 6-8)
          # Substring for PAM is selected based on the Editing window:
          
          PerfectUpstreamPAM <- as.character(substr(Seq_minus_char_NGG, 39,43))  

          # 'minus_perfecteditable' object is created to store found PAM candidates:
          
          minus_perfecteditable <- NULL
          
          # The PerfectUpstreamPAM substring is searched for all different PAM sequences
          # Matches are stored in the 'minus_perfecteditable' object
          
          for(i in 1:length(PerfectUpstreamPAM))
          {
            if (grepl( "AGG", PerfectUpstreamPAM[i]) == TRUE || grepl( "GGG", PerfectUpstreamPAM[i]) == TRUE || grepl( "CGG", PerfectUpstreamPAM[i]) == TRUE || grepl( "TGG", PerfectUpstreamPAM[i]) == TRUE ){    
              #print(c(i,PerfectUpstreamPAM[i]))
              tmp1 <- c(i,PerfectUpstreamPAM[i])
              minus_perfecteditable <- rbind(minus_perfecteditable,tmp1)
              
              
            }
          }
          
          # Quality check to only run code, if suitable PAM sequences are found:
          
          if(is.null(minus_perfecteditable) != TRUE) {
            
            minus_perfecteditabledf <- data.frame(minus_perfecteditable)
            
            minus_perfecteditabledf$X2 <- as.character(minus_perfecteditabledf$X2)
            
            minus_perfecteditabledf$X1 <- as.numeric(nrow(minus_perfecteditable))

            
            # Guide sequences are designed based on the position of the PAM sequence in the 'minus_perfecteditable' object:
            
            for (i in 1:nrow(minus_perfecteditabledf))
              for (j in 1:(nchar(minus_perfecteditabledf[1,2])-2))
              {
                {
                  if (grepl( "AGG", substr(minus_perfecteditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(minus_perfecteditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(minus_perfecteditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(minus_perfecteditabledf[i,2], j,j+2)) == TRUE ){
                    
                    # Guide vector is created with the guide number, Protospacer sequence, Edit position, PAM sequence, Base editor
                    # Resulting vectors are combined to create 'guides'
                    
                    Guide <- c(minus_perfecteditabledf[i,1],substr(Seq_minus_char_NGG[minus_perfecteditabledf[i,1]],18+j,37+j),26-(17+j),substr(minus_perfecteditabledf[i,2], j,j+2),"ABEmax/ABE8e")
                    
                    tmp2 <- NULL
                    
                    Download <- rbind(Download,Guide[2])
                    
                    for (l in 1:(nchar(minus_perfecteditabledf[1,2])-2))
                      
                    {
                      if(grepl("A", substr(Guide[2],5+l,5+l)) == TRUE)
                      {
                        
                        tmp <- 5+l
                        tmp2 <- rbind(tmp2,tmp)
                        print("yeha")
                        
                      }else{
                        
                        print("bla")
                      }
                    }
                    
                    
                    tmp2df <- data.frame(tmp2)
                    tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                    tmp3 <- data.frame(tmp3)
                    tmp4 <- data.frame(tmp3)
                    
                    
                    nrow(tmp3)
                    
                    t <- 0
                    #Guide[6] <- NULL
                    
                    for (b in 1:nrow(tmp3)) {
                      
                      if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                      {
                        Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                        t <- t+43
                        tmp4 <- tmp4[] + 43
                        
                      }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                      {
                        
                        Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                        tmp4 <- tmp4[] + 42
                        
                        
                      }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                      {
                        
                        
                        Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                        tmp4 <- tmp4[] + 43
                        
                      }
                      
                      
                      
                    }
                    
                    guides <- rbind(guides,Guide)
                    
                    
                    
                  }
                }
              }
            
          }
          
          
         }
        
          # If check of wildtype sequence with input'-SNP' does not match, a gets assigned value of NULL
          
          else if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) != TRUE){
           
           a <- NULL
           
         }
          
          
        }
        
        # For Base editing guides on the Plus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "+" || input$SNP == "C>T" && input$`Gene Orientation2` == "-" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "A>G"){
          
          if (("G" == substr(Seq_plus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_plus_char_NGG, start = 26, stop = 26) <- "A"
            
          
          
          # For Base editing guides on the Plus Strand:
            
          # Perfect edible window: For ABEmax (Editing Position: 6-8)
          # Substring for PAM is selected based on the Editing window:
          
          PerfectDownstreamPAM <- as.character(substr(Seq_plus_char_NGG, 39,43))  

          
          ## Test if there is a suitable Downstream PAM sequence
          
          plus_perfecteditable <- NULL
          
          for(i in 1:length(PerfectDownstreamPAM))
          {
            if (grepl( "AGG", PerfectDownstreamPAM[i]) == TRUE || grepl( "GGG", PerfectDownstreamPAM[i]) == TRUE || grepl( "CGG", PerfectDownstreamPAM[i]) == TRUE || grepl( "TGG", PerfectDownstreamPAM[i]) == TRUE ){    
              #print(c(i,PerfectDownstreamPAM[i]))
              tmp <- c(i,PerfectDownstreamPAM[i])
              plus_perfecteditable <- rbind(plus_perfecteditable,tmp)
              
              
            }
          }
          
          if(is.null(plus_perfecteditable) != TRUE) {
            
            
            plus_perfecteditabledf <- data.frame(plus_perfecteditable)
            
            plus_perfecteditabledf$X2 <- as.character(plus_perfecteditabledf$X2)
            
            plus_perfecteditabledf$X1 <- as.numeric(nrow(plus_perfecteditable))
            
            
            
            
            for (i in 1:nrow(plus_perfecteditabledf))
              for (j in 1:(nchar(plus_perfecteditabledf[1,2])-2))
              {
                {
                  if (grepl( "AGG", substr(plus_perfecteditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(plus_perfecteditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(plus_perfecteditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(plus_perfecteditabledf[i,2], j,j+2)) == TRUE ){
                    
                    #print(plus_perfecteditabledf[i,2])
                    #print(FANCA_char[perfect_editable_sense[i,1]])
                    
                    Guide <- c(plus_perfecteditabledf[i,1],substr(Seq_plus_char_NGG[plus_perfecteditabledf[i,1]],18+j,37+j),26-(17+j),substr(plus_perfecteditabledf[i,2], j,j+2),"ABEmax/ABE8e")
                  
                    tmp2 <- NULL
                    
                    Download <- rbind(Download,Guide[2])
                    
                    for (l in 1:(nchar(plus_perfecteditabledf[1,2])-2))
                      
                    {
                      if(grepl("A", substr(Guide[2],5+l,5+l)) == TRUE)
                      {
                        
                        tmp <- 5+l
                        tmp2 <- rbind(tmp2,tmp)
                        print("yeha")
                        
                      }else{
                        
                        print("bla")
                      }
                    }
      
                    
                    tmp2df <- data.frame(tmp2)
                    tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                    tmp3 <- data.frame(tmp3)
                    tmp4 <- data.frame(tmp3)
                    
                    
                    nrow(tmp3)
                    
                    t <- 0
                    #Guide[6] <- NULL
                    
                    for (b in 1:nrow(tmp3)) {
                      
                      if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                      {
                        Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                        t <- t+43
                        tmp4 <- tmp4[] + 43
                        
                      }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                      {
                        
                        Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                        tmp4 <- tmp4[] + 42
                        
                        
                      }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                      {
                        
                        
                        Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                        tmp4 <- tmp4[] + 43
                        
                      }
                      
                      
                      
                    }
                    
                    
                    guides <- rbind(guides,Guide)
                    
                    
                    
                  }
                }
              }
            
            
          }
          }else{
            
            a <- NULL
            
          }
          
        }
        
#---------------------------------------------#        
        
        # Guide search For NG-ABEmax, VRQR-ABEmax (NGA PAM) depending on input'-SNP' and '-Gene Orientation':
        # For Base editing guides on the Minus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "-" || input$SNP == "C>T" && input$`Gene Orientation2` == "+" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "T>C" ){
          
          
          if (input$SequenceInput2 == "Genomic coordinates" && ("G" == substr(Seq_minus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_minus_char_NGA, start = 26, stop = 26) <- "A"
            
          
          
          ## Perfect edible: NG-ABEmax window 4-6 (PAM: NGA)
          
          PerfectUpstreamPAM_NGA <- substr(Seq_minus_char_NGA, 41,45)  
          

          ## Target Plus Strand C>T or Minus Strand G>A
          
          minus_perfecteditable_NGA <- NULL;
          
          for(i in 1:length(PerfectUpstreamPAM_NGA))
          {
            if (grepl( "AGA", PerfectUpstreamPAM_NGA[i]) == TRUE || grepl( "GGA", PerfectUpstreamPAM_NGA[i]) == TRUE || grepl( "CGA", PerfectUpstreamPAM_NGA[i]) == TRUE || grepl( "TGA", PerfectUpstreamPAM_NGA[i]) == TRUE ){    
              #print(c(i,PerfectUpstreamPAM_NGA[i]))
              tmp <- c(i,PerfectUpstreamPAM_NGA[i])
              minus_perfecteditable_NGA <- rbind(minus_perfecteditable_NGA,tmp)
              
              
            }
          }
          
          if(is.null(minus_perfecteditable_NGA) != TRUE) {
            
            minus_perfecteditabledf_NGA <- data.frame(minus_perfecteditable_NGA)
            
            minus_perfecteditabledf_NGA$X2 <- as.character(minus_perfecteditabledf_NGA$X2)
            
            minus_perfecteditabledf_NGA$X1 <- as.numeric(nrow(minus_perfecteditable_NGA))
            
            
            for (i in 1:nrow(minus_perfecteditabledf_NGA))
              for (j in 1:(nchar(minus_perfecteditabledf_NGA[1,2])-2))
              {
                {
                  if (grepl( "AGA", substr(minus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE || grepl( "GGA", substr(minus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE || grepl( "CGA", substr(minus_perfecteditabledf_NGA[i,2],j,j+2)) == TRUE || grepl( "TGA", substr(minus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE ){
                    
                    #print(minus_perfecteditabledf_NGA[i,2])
                    #print(FANCA_char[perfect_editable_sense[i,1]])
                    Guide <- c(minus_perfecteditabledf_NGA[i,1],substr(Seq_minus_char_NGA[minus_perfecteditabledf_NGA[i,1]],20+j,39+j),26-(19+j),substr(minus_perfecteditabledf_NGA[i,2], j,j+2),"NG-ABEmax/ABE8e, VRQR-ABEmax")
                    
                    
                    tmp2 <- NULL
                    
                    Download <- rbind(Download,Guide[2])
                    
                    for (l in 1:(nchar(minus_perfecteditabledf_NGA[1,2])-2))
                      
                    ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 3
                      
                    {
                      if(grepl("A", substr(Guide[2],3+l,3+l)) == TRUE)
                      {
                        
                        tmp <- 3+l
                        tmp2 <- rbind(tmp2,tmp)
                        #print("yeha")
                        
                      }else{
                        
                        #print("bla")
                      }
                    }
                    
                    
                    tmp2df <- data.frame(tmp2)
                    tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                    tmp3 <- data.frame(tmp3)
                    tmp4 <- data.frame(tmp3)
                    
                    t <- 0
                    #Guide[6] <- NULL
                    
                    for (b in 1:nrow(tmp3)) {
                      
                      if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                      {
                        Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                        t <- t+43
                        tmp4 <- tmp4[] + 43
                        
                      }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                      {
                        
                        Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                        tmp4 <- tmp4[] + 42
                        
                        
                      }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                      {
                        
                        
                        Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                        tmp4 <- tmp4[] + 43
                        
                      }
                      
                      
                      
                    }
                    
                    
                    guides <- rbind(guides,Guide)
                    
                    
                    
                    
                  }
                }
              }
            
          }
          }else{
            
            a <- NULL
            
          }
          
        }
        
        # For Base editing guides on the Plus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "+" || input$SNP == "C>T" && input$`Gene Orientation2` == "-" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "A>G"){
          
          
          if (("G" == substr(Seq_plus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_plus_char_NGA, start = 26, stop = 26) <- "A"
            
          
          
          ## Perfect edible: NG-ABEmax window 4-6 (PAM: NGA)
          
          PerfectDownstreamPAM_NGA <- substr(Seq_plus_char_NGA, 41,45) 

          ## Target Plus Strand G>A or Minus Strand C>T (with NG-ABEmax)
          
          
          plus_perfecteditable_NGA <- NULL;
          
          for(i in 1:length(PerfectDownstreamPAM_NGA))
          {
            if (grepl( "AGA", PerfectDownstreamPAM_NGA[i]) == TRUE || grepl( "GGA", PerfectDownstreamPAM_NGA[i]) == TRUE || grepl( "CGA", PerfectDownstreamPAM_NGA[i]) == TRUE || grepl( "TGA", PerfectDownstreamPAM_NGA[i]) == TRUE ){    
              #print(c(i,PerfectDownstreamPAM_NGA[i]))
              tmp <- c(i,PerfectDownstreamPAM_NGA[i])
              plus_perfecteditable_NGA <- rbind(plus_perfecteditable_NGA,tmp)
              
              
            }
          }
          
          
          if(is.null(plus_perfecteditable_NGA) != TRUE) {
            
            plus_perfecteditabledf_NGA <- data.frame(plus_perfecteditable_NGA)
            
            plus_perfecteditabledf_NGA$X2 <- as.character(plus_perfecteditabledf_NGA$X2)
            
            plus_perfecteditabledf_NGA$X1 <- as.numeric(nrow(plus_perfecteditable_NGA))
            
        
            for (i in 1:nrow(plus_perfecteditabledf_NGA))
              for (j in 1:(nchar(plus_perfecteditabledf_NGA[1,2])-2))
              {
                {
                  if (grepl( "AGA", substr(plus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE || grepl( "GGA", substr(plus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE || grepl( "CGA", substr(plus_perfecteditabledf_NGA[i,2],j,j+2)) == TRUE || grepl( "TGA", substr(plus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE ){
                    
                    #print(plus_perfecteditabledf_NGA[i,2])
                    #print(FANCA_char[perfect_editable_sense[i,1]])
                    Guide <- c(plus_perfecteditabledf_NGA[i,1],substr(Seq_plus_char_NGA[plus_perfecteditabledf_NGA[i,1]],20+j,39+j),26-(19+j),substr(plus_perfecteditabledf_NGA[i,2], j,j+2),"NG-ABEmax/ABE8e, VRQR-ABEmax")
                    tmp2 <- NULL
                    
                    Download <- rbind(Download,Guide[2])
                    
                    for (l in 1:(nchar(plus_perfecteditabledf_NGA[1,2])-2))
                      
                      ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 3
                      
                    {
                      if(grepl("A", substr(Guide[2],3+l,3+l)) == TRUE)
                      {
                        
                        tmp <- 3+l
                        tmp2 <- rbind(tmp2,tmp)
                        print("yeha")
                        
                      }else{
                        
                        print("bla")
                      }
                    }
                    
                    
                    tmp2df <- data.frame(tmp2)
                    tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                    tmp3 <- data.frame(tmp3)
                    tmp4 <- data.frame(tmp3)
                    
                    
                    nrow(tmp3)
                    
                    t <- 0
                    #Guide[6] <- NULL
                    
                    for (b in 1:nrow(tmp3)) {
                      
                      if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                      {
                        Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                        t <- t+43
                        tmp4 <- tmp4[] + 43
                        
                      }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                      {
                        
                        Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                        tmp4 <- tmp4[] + 42
                        
                        
                      }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                      {
                        
                        
                        Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                        tmp4 <- tmp4[] + 43
                        
                      }
                      
                      
                      
                    }
                    
                    
                    guides <- rbind(guides,Guide)
                    
                    
                    
                    
                  }
                }
              }
            
          }
          }else{
            
            a <- NULL
            
          }
          
          
        }

#---------------------------------------------#
        
        # Guide search For NG-ABEmax (NGCG PAM) depending on input'-SNP' and '-Gene Orientation':
        # For Base editing guides on the Minus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "-" || input$SNP == "C>T" && input$`Gene Orientation2` == "+" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "T>C"){
          
          # Test wildtype sequence if input'-SNP' was assigned right and if so includes patient mutation into character string:
          # If not, see code line 216
          
          if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_minus_char_NGCG, start = 26, stop = 26) <- "A"
            
            # patient character string is now searched for suitable PAM sequences:
            
            # Perfect edible window: For NG-ABEmax NGCG PAM (Editing Position: 5-7)
            # Substring for PAM is selected based on the Editing window:
            
            PerfectUpstreamPAM_NGCG <- as.character(substr(Seq_minus_char_NGCG, 40,45))  
            
            # 'minus_perfecteditable' object is created to store found PAM candidates:
            
            minus_perfecteditable_NGCG <- NULL
            
            # The PerfectUpstreamPAM_NGCG substring is searched for all different PAM sequences
            # Matches are stored in the 'minus_perfecteditable_NGCG' object
            
            for(i in 1:length(PerfectUpstreamPAM_NGCG))
            {
              if (grepl( "AGCG", PerfectUpstreamPAM_NGCG[i]) == TRUE || grepl( "GGCG", PerfectUpstreamPAM_NGCG[i]) == TRUE || grepl( "CGCG", PerfectUpstreamPAM_NGCG[i]) == TRUE || grepl( "TGCG", PerfectUpstreamPAM_NGCG[i]) == TRUE ){    
                tmp1 <- c(i,PerfectUpstreamPAM_NGCG[i])
                minus_perfecteditable_NGCG <- rbind(minus_perfecteditable_NGCG,tmp1)
                
                
              }
            }
            
            # Quality check to only run code, if suitable PAM sequences are found:
            
            if(is.null(minus_perfecteditable_NGCG) != TRUE) {
              
              minus_perfecteditable_NGCGdf <- data.frame(minus_perfecteditable_NGCG)
              
              minus_perfecteditable_NGCGdf$X2 <- as.character(minus_perfecteditable_NGCGdf$X2)
              
              minus_perfecteditable_NGCGdf$X1 <- as.numeric(nrow(minus_perfecteditable_NGCG))
              
              
              # Guide sequences are designed based on the position of the PAM sequence in the 'minus_perfecteditable_NGCG' object:
              
              for (i in 1:nrow(minus_perfecteditable_NGCGdf))
                for (j in 1:(nchar(minus_perfecteditable_NGCGdf[1,2])-3))
                {
                  {
                    if (grepl( "AGCG", substr(minus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE || grepl( "GGCG", substr(minus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE || grepl( "CGCG", substr(minus_perfecteditable_NGCGdf[i,2],j,j+3)) == TRUE || grepl( "TGCG", substr(minus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE ){
                      
                      # Guide vector is created with the guide number, Protospacer sequence, Edit position, PAM sequence, Base editor
                      # Resulting vectors are combined to create 'guides'
                      
                      Guide <- c(minus_perfecteditable_NGCGdf[i,1],substr(Seq_minus_char_NGCG[minus_perfecteditable_NGCGdf[i,1]],19+j,38+j),25-(17+j),substr(minus_perfecteditable_NGCGdf[i,2], j,j+3),"NG-ABEmax/ABE8e")
                      tmp2 <- NULL
                      
                      Download <- rbind(Download,Guide[2])
                      
                      for (l in 1:(nchar(minus_perfecteditable_NGCGdf[1,2])-3))
                        
                        ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                        
                      {
                        if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                        {
                          
                          tmp <- 4+l
                          tmp2 <- rbind(tmp2,tmp)
                          #print("yeha")
                          
                        }else{
                          
                          #print("bla")
                        }
                      }
                      
                      
                      tmp2df <- data.frame(tmp2)
                      tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                      tmp3 <- data.frame(tmp3)
                      tmp4 <- data.frame(tmp3)
                      
                      t <- 0
                      #Guide[6] <- NULL
                      
                      for (b in 1:nrow(tmp3)) {
                        
                        if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                        {
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          t <- t+43
                          tmp4 <- tmp4[] + 43
                          
                        }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                        {
                          
                          Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                          tmp4 <- tmp4[] + 42
                          
                          
                        }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                        {
                          
                          
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          tmp4 <- tmp4[] + 43
                          
                        }
                        
                        
                        
                      }
                      
                      
                      guides <- rbind(guides,Guide)
                      
                      
                      
                      
                    }
                  }
                }
              
            }
            
            
          }
          
          # If check of wildtype sequence with input'-SNP' does not match, a gets assigned value of NULL
          
          else if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) != TRUE){
            
            a <- NULL
            
          }
          
          
        }
        
        # For Base editing guides on the Plus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "+" || input$SNP == "C>T" && input$`Gene Orientation2` == "-" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "A>G"){
          
          if (("G" == substr(Seq_plus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_plus_char_NGCG, start = 26, stop = 26) <- "A"
            
            # For Base editing guides on the Plus Strand:
            
            # Perfect edible window: For NG-ABEmax NGCG PAM (Editing Position: 5-7)
            # Substring for PAM is selected based on the Editing window:
            
            PerfectDownstreamPAM_NGCG <- as.character(substr(Seq_plus_char_NGCG, 40,45)) 

            ## Test if there is a suitable Downstream PAM sequence
            
            plus_perfecteditable_NGCG <- NULL
            
            for(i in 1:length(PerfectDownstreamPAM_NGCG))
            {
              if (grepl( "AGCG", PerfectDownstreamPAM_NGCG[i]) == TRUE || grepl( "GGCG", PerfectDownstreamPAM_NGCG[i]) == TRUE || grepl( "CGCG", PerfectDownstreamPAM_NGCG[i]) == TRUE || grepl( "TGCG", PerfectDownstreamPAM_NGCG[i]) == TRUE ){    
                #print(c(i,PerfectDownstreamPAM_NGCG[i]))
                tmp <- c(i,PerfectDownstreamPAM_NGCG[i])
                plus_perfecteditable_NGCG <- rbind(plus_perfecteditable_NGCG,tmp)
                
                
              }
            }
            
            if(is.null(plus_perfecteditable_NGCG) != TRUE) {
              
              
              plus_perfecteditable_NGCGdf <- data.frame(plus_perfecteditable_NGCG)
              
              plus_perfecteditable_NGCGdf$X2 <- as.character(plus_perfecteditable_NGCGdf$X2)
              
              plus_perfecteditable_NGCGdf$X1 <- as.numeric(nrow(plus_perfecteditable_NGCG))
              
              
              
              
              for (i in 1:nrow(plus_perfecteditable_NGCGdf))
                for (j in 1:(nchar(plus_perfecteditable_NGCGdf[1,2])-3))
                {
                  {
                    if (grepl( "AGCG", substr(plus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE || grepl( "GGCG", substr(plus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE || grepl( "CGCG", substr(plus_perfecteditable_NGCGdf[i,2],j,j+3)) == TRUE || grepl( "TGCG", substr(plus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE ){
                      
                      
                      
                      #print(plus_perfecteditable_NGCGdf[i,2])
                      #print(FANCA_char[perfect_editable_sense[i,1]])
                      Guide <- c(plus_perfecteditable_NGCGdf[i,1],substr(Seq_plus_char_NGCG[plus_perfecteditable_NGCGdf[i,1]],19+j,38+j),25-(17+j),substr(plus_perfecteditable_NGCGdf[i,2], j,j+3),"NG-ABEmax/ABE8e")
                      tmp2 <- NULL
                      
                      Download <- rbind(Download,Guide[2])
                      
                      for (l in 1:(nchar(plus_perfecteditable_NGCGdf[1,2])-3))
                        
                        ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                        
                      {
                        if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                        {
                          
                          tmp <- 4+l
                          tmp2 <- rbind(tmp2,tmp)
                          #print("yeha")
                          
                        }else{
                          
                          #print("bla")
                        }
                      }
                      
                      
                      tmp2df <- data.frame(tmp2)
                      tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                      tmp3 <- data.frame(tmp3)
                      tmp4 <- data.frame(tmp3)
                      
                      t <- 0
                      #Guide[6] <- NULL
                      
                      for (b in 1:nrow(tmp3)) {
                        
                        if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                        {
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          t <- t+43
                          tmp4 <- tmp4[] + 43
                          
                        }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                        {
                          
                          Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                          tmp4 <- tmp4[] + 42
                          
                          
                        }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                        {
                          
                          
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          tmp4 <- tmp4[] + 43
                          
                        }
                        
                        
                        
                      }
                      
                      
                      guides <- rbind(guides,Guide)
                      
                      
                      
                      
                    }
                  }
                }
              
            }
          }
          else{
            
            a <- NULL
            
          }
          
        } 
        
#---------------------------------------------#     
        # Guide search For Sa-ABEmax (NNGRRT PAM) depending on input'-SNP' and '-Gene Orientation':
        # For Base editing guides on the Minus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "-" || input$SNP == "C>T" && input$`Gene Orientation2` == "+" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "T>C"){
          
          # Test wildtype sequence if input'-SNP' was assigned right and if so includes patient mutation into character string:
          # If not, see code line 216
          
          if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_minus_char_NNGRRT, start = 26, stop = 26) <- "A"
            
            # patient character string is now searched for suitable PAM sequences:
            
            # Perfect edible window: For Sa-ABEmax  NNGRRT PAM (Editing Position: 3-13)
            # Substring for PAM is selected based on the Editing window:
            
            PerfectUpstreamPAM_NNGRRT <- as.character(substr(Seq_minus_char_NNGRRT, 36,49))  
            
            print(PerfectUpstreamPAM_NNGRRT)
            
            
            # 'minus_perfecteditable_NNGRRT' object is created to store found PAM candidates:
            
            minus_perfecteditable_NNGRRT <- NULL
            
            # The PerfectUpstreamPAM_NNGRRT substring is searched for all different PAM sequences
            # Matches are stored in the 'minus_perfecteditable_NNGRRT' object
            
            for(i in 1:length(PerfectUpstreamPAM_NNGRRT))
            {
              if (grepl( "GAAT", PerfectUpstreamPAM_NNGRRT[i]) == TRUE || grepl( "GAGT", PerfectUpstreamPAM_NNGRRT[i]) == TRUE || grepl( "GGGT", PerfectUpstreamPAM_NNGRRT[i]) == TRUE || grepl( "GGAT", PerfectUpstreamPAM_NNGRRT[i]) == TRUE ){    
                #print(c(i,PerfectUpstreamPAM_NNGRRT[i]))
                tmp1 <- c(i,PerfectUpstreamPAM_NNGRRT[i])
                minus_perfecteditable_NNGRRT <- rbind(minus_perfecteditable_NNGRRT,tmp1)
                
                
              }
            }
            
            # Quality check to only run code, if suitable PAM sequences are found:
            
            if(is.null(minus_perfecteditable_NNGRRT) != TRUE) {
              
              minus_perfecteditable_NNGRRTdf <- data.frame(minus_perfecteditable_NNGRRT)
              
              minus_perfecteditable_NNGRRTdf$X2 <- as.character(minus_perfecteditable_NNGRRTdf$X2)
              
              minus_perfecteditable_NNGRRTdf$X1 <- as.numeric(nrow(minus_perfecteditable_NNGRRT))
              
              
              # Guide sequences are designed based on the position of the PAM sequence in the 'minus_perfecteditable_NNGRRT' object:
              
              for (i in 1:nrow(minus_perfecteditable_NNGRRTdf))
                for (j in 1:(nchar(minus_perfecteditable_NNGRRTdf[1,2])-3))
                {
                  {
                    if (grepl( "GAAT", substr(minus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE || grepl( "GAGT", substr(minus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE || grepl( "GGGT", substr(minus_perfecteditable_NNGRRTdf[i,2],j,j+3)) == TRUE || grepl( "GGAT", substr(minus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE ){
                      
                      # Guide vector is created with the guide number, Protospacer sequence, Edit position, PAM sequence, Base editor
                      # Resulting vectors are combined to create 'guides'
                      
                      Guide <- c(minus_perfecteditable_NNGRRTdf[i,1],substr(Seq_minus_char_NGCG[minus_perfecteditable_NNGRRTdf[i,1]],13+j,32+j),31-(17+j),substr(Seq_minus_char_NGCG[minus_perfecteditable_NNGRRTdf[i,1]],33+j,38+j),"Sa-ABEmax/ABE8e")
                      tmp2 <- NULL
                      
                      Download <- rbind(Download,Guide[2])
                      
                      for (l in 1:(nchar(minus_perfecteditable_NNGRRTdf[1,2])-3))
                        
                        ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 2
                        
                      {
                        if(grepl("A", substr(Guide[2],2+l,2+l)) == TRUE)
                        {
                          
                          tmp <- 2+l
                          tmp2 <- rbind(tmp2,tmp)
                          #print("yeha")
                          
                        }else{
                          
                          #print("bla")
                        }
                      }
                      
                      
                      tmp2df <- data.frame(tmp2)
                      tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                      tmp3 <- data.frame(tmp3)
                      tmp4 <- data.frame(tmp3)
                      
                      t <- 0
                      #Guide[6] <- NULL
                      
                      for (b in 1:nrow(tmp3)) {
                        
                        if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                        {
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          t <- t+43
                          tmp4 <- tmp4[] + 43
                          
                        }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                        {
                          
                          Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                          tmp4 <- tmp4[] + 42
                          
                          
                        }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                        {
                          
                          
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          tmp4 <- tmp4[] + 43
                          
                        }
                        
                        
                        
                      }
                      
                      
                      guides <- rbind(guides,Guide)
                      
                      
                      
                    }
                  }
                }
              
            }
            
            
          }
          
          # If check of wildtype sequence with input'-SNP' does not match, a gets assigned value of NULL
          
          else if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) != TRUE){
            
            a <- NULL
            
          }
          
          
        }
        
        # For Base editing guides on the Plus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "+" || input$SNP == "C>T" && input$`Gene Orientation2` == "-" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "A>G"){
          
          if (("G" == substr(Seq_plus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_plus_char_NNGRRT, start = 26, stop = 26) <- "A"
            
            # For Base editing guides on the Plus Strand:
            
            # Perfect edible window: For Sa-ABEmax  NNGRRT PAM (Editing Position: 3-13)
            # Substring for PAM is selected based on the Editing window:
            
            PerfectDownstreamPAM_NNGRRT <- as.character(substr(Seq_plus_char_NNGRRT, 36, 49))  
            
            ## Test if there is a suitable Downstream PAM sequence
            
            plus_perfecteditable_NNGRRT <- NULL
            
            for(i in 1:length(PerfectDownstreamPAM_NNGRRT))
            {
              if (grepl( "GAAT", PerfectDownstreamPAM_NNGRRT[i]) == TRUE || grepl( "GAGT", PerfectDownstreamPAM_NNGRRT[i]) == TRUE || grepl( "GGGT", PerfectDownstreamPAM_NNGRRT[i]) == TRUE || grepl( "GGAT", PerfectDownstreamPAM_NNGRRT[i]) == TRUE ){    
                tmp <- c(i,PerfectDownstreamPAM_NNGRRT[i])
                plus_perfecteditable_NNGRRT <- rbind(plus_perfecteditable_NNGRRT,tmp)
                
                
              }
            }
            
            if(is.null(plus_perfecteditable_NNGRRT) != TRUE) {
              
              
              plus_perfecteditable_NNGRRTdf <- data.frame(plus_perfecteditable_NNGRRT)
              
              plus_perfecteditable_NNGRRTdf$X2 <- as.character(plus_perfecteditable_NNGRRTdf$X2)
              
              plus_perfecteditable_NNGRRTdf$X1 <- as.numeric(nrow(plus_perfecteditable_NNGRRT))
              
              
              
              
              for (i in 1:nrow(plus_perfecteditable_NNGRRTdf))
                for (j in 1:(nchar(plus_perfecteditable_NNGRRTdf[1,2])-3))
                {
                  {
                    if (grepl( "GAAT", substr(plus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE || grepl( "GAGT", substr(plus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE || grepl( "GGGT", substr(plus_perfecteditable_NNGRRTdf[i,2],j,j+3)) == TRUE || grepl( "GGAT", substr(plus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE ){
                      
                      Guide <- c(plus_perfecteditable_NNGRRTdf[i,1],substr(Seq_plus_char_NNGRRT[plus_perfecteditable_NNGRRTdf[i,1]],13+j,32+j),31-(17+j),substr(Seq_plus_char_NNGRRT[plus_perfecteditable_NNGRRTdf[i,1]],33+j,38+j),"Sa-ABEmax/ABE8e")
                      tmp2 <- NULL
                      
                      Download <- rbind(Download,Guide[2])
                      
                      for (l in 1:(nchar(plus_perfecteditable_NNGRRTdf[1,2])-3))
                        
                        ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 2
                        
                      {
                        if(grepl("A", substr(Guide[2],2+l,2+l)) == TRUE)
                        {
                          
                          tmp <- 2+l
                          tmp2 <- rbind(tmp2,tmp)
                          #print("yeha")
                          
                        }else{
                          
                          #print("bla")
                        }
                      }
                      
                      
                      tmp2df <- data.frame(tmp2)
                      tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                      tmp3 <- data.frame(tmp3)
                      tmp4 <- data.frame(tmp3)
                      
                      t <- 0
                      #Guide[6] <- NULL
                      
                      for (b in 1:nrow(tmp3)) {
                        
                        if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                        {
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          t <- t+43
                          tmp4 <- tmp4[] + 43
                          
                        }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                        {
                          
                          Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                          tmp4 <- tmp4[] + 42
                          
                          
                        }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                        {
                          
                          
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          tmp4 <- tmp4[] + 43
                          
                        }
                        
                        
                        
                      }
                      
                      
                      guides <- rbind(guides,Guide)
                      
                      
                      
                      
                    }
                  }
                }
              
            }
          }else{
            
            a <- NULL
            
          }
          
        }
        
#---------------------------------------------#        
        
        # Guide search For SaKKH-ABEmax (NNHRRT PAM) depending on input'-SNP' and '-Gene Orientation':
        # For Base editing guides on the Minus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "-" || input$SNP == "C>T" && input$`Gene Orientation2` == "+" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "T>C"){
          
          # Test wildtype sequence if input'-SNP' was assigned right and if so includes patient mutation into character string:
          # If not, see code line 216
          
          if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_minus_char_NNHRRT, start = 26, stop = 26) <- "A"
            
            # patient character string is now searched for suitable PAM sequences:
            
            # Perfect edible window: For SaKKH-ABEmax  NNHRRT PAM (Editing Position: 4; 6-12)
            # Substring for PAM is selected based on the Editing window:
            
            PerfectUpstreamPAM_NNHRRT <- as.character(substr(Seq_minus_char_NNHRRT, 37,48))  
            
            # 'minus_perfecteditable_NNHRRT' object is created to store found PAM candidates:
            
            minus_perfecteditable_NNHRRT <- NULL
            
            # The PerfectUpstreamPAM_NNHRRT substring is searched for all different PAM sequences
            # Matches are stored in the 'minus_perfecteditable_NNHRRT' object
            
            for(i in 1:length(PerfectUpstreamPAM_NNHRRT))
            {
              if (grepl( "GAAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "GAGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "GGGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "GGAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "TAAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "TAGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "TGGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "TGAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "AAAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "AAGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "AGGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "AGAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE ){    
                tmp1 <- c(i,PerfectUpstreamPAM_NNHRRT[i])
                minus_perfecteditable_NNHRRT <- rbind(minus_perfecteditable_NNHRRT,tmp1)
                
                
              }
            }
            
            # Quality check to only run code, if suitable PAM sequences are found:
            
            if(is.null(minus_perfecteditable_NNHRRT) != TRUE) {
              
              minus_perfecteditable_NNHRRTdf <- data.frame(minus_perfecteditable_NNHRRT)
              
              minus_perfecteditable_NNHRRTdf$X2 <- as.character(minus_perfecteditable_NNHRRTdf$X2)
              
              minus_perfecteditable_NNHRRTdf$X1 <- as.numeric(nrow(minus_perfecteditable_NNHRRT))
              
              
              # Guide sequences are designed based on the position of the PAM sequence in the 'minus_perfecteditable_NNHRRT' object:
              
              for (i in 1:nrow(minus_perfecteditable_NNHRRTdf))
                for (j in 1:(nchar(minus_perfecteditable_NNHRRTdf[1,2])-3))
                {
                  {
                    if (grepl("AGGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AAGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AGAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AAAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TGAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl("TGGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TAGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TAAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("GAAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl( "GAGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl( "GGGT", substr(minus_perfecteditable_NNHRRTdf[i,2],j,j+3)) == TRUE || grepl( "GGAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ){
                      
                      # Guide vector is created with the guide number, Protospacer sequence, Edit position, PAM sequence, Base editor
                      # Resulting vectors are combined to create 'guides'
                      
                      if(j != 2){
                        
                        Guide <- c(minus_perfecteditable_NNHRRTdf[i,1],substr(Seq_minus_char_NNHRRT[minus_perfecteditable_NNHRRTdf[i,1]],14+j,33+j),30-(17+j),substr(Seq_minus_char_NNHRRT[minus_perfecteditable_NNHRRTdf[i,1]],34+j,39+j),"SaKKH-ABEmax/ABE8e")
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(minus_perfecteditable_NNHRRTdf[1,2])-3))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 3
                          
                        {
                          if(grepl("A", substr(Guide[2],3+l,3+l)) == TRUE)
                          {
                            
                            tmp <- 3+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                      }else{
                        
                        guides <- guides
                      }
                      
                      
                    }
                    
                    
                  }
                  
                }
              
              
            }
            
          }
          # If check of wildtype sequence with input'-SNP' does not match, a gets assigned value of NULL
          
          else if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) != TRUE){
            
            a <- NULL
            
          }
        }
        
        # For Base editing guides on the Plus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "+" || input$SNP == "C>T" && input$`Gene Orientation2` == "-" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "A>G"){
          
          if (("G" == substr(Seq_plus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input" )
          {
            substr(Seq_plus_char_NNHRRT, start = 26, stop = 26) <- "A"
            
            
            
            # For Base editing guides on the Plus Strand:
            
            # Perfect edible window: For SaKKH-ABEmax  NNHRRT PAM (Editing Position: 4; 6-12)
            # Substring for PAM is selected based on the Editing window:
            
            PerfectDownstreamPAM_NNHRRT <- as.character(substr(Seq_plus_char_NNHRRT, 37, 48))  
            
            ## Test if there is a suitable Downstream PAM sequence
            
            plus_perfecteditable_NNHRRT <- NULL
            
            for(i in 1:length(PerfectDownstreamPAM_NNHRRT))
            {
              if (grepl( "GAAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "GAGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "GGGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "GGAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "TAAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "TAGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "TGGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "TGAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "AAAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "AAGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "AGGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "AGAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE ){    
                #print(c(i,PerfectDownstreamPAM_NNHRRT[i]))
                tmp <- c(i,PerfectDownstreamPAM_NNHRRT[i])
                plus_perfecteditable_NNHRRT <- rbind(plus_perfecteditable_NNHRRT,tmp)
                
                
              }
            }
            
            if(is.null(plus_perfecteditable_NNHRRT) != TRUE) {
              
              
              plus_perfecteditable_NNHRRTdf <- data.frame(plus_perfecteditable_NNHRRT)
              
              plus_perfecteditable_NNHRRTdf$X2 <- as.character(plus_perfecteditable_NNHRRTdf$X2)
              
              plus_perfecteditable_NNHRRTdf$X1 <- as.numeric(nrow(plus_perfecteditable_NNHRRT))
              
              
              
              
              for (i in 1:nrow(plus_perfecteditable_NNHRRTdf))
                for (j in 1:(nchar(plus_perfecteditable_NNHRRTdf[1,2])-3))
                {
                  {
                    if (grepl("AGGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AAGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AGAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AAAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TGAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl("TGGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TAGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TAAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("GAAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl( "GAGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl( "GGGT", substr(plus_perfecteditable_NNHRRTdf[i,2],j,j+3)) == TRUE || grepl( "GGAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ){
                      
                      if(j != 2){
                        
                        Guide <- c(plus_perfecteditable_NNHRRTdf[i,1],substr(Seq_plus_char_NNHRRT[plus_perfecteditable_NNHRRTdf[i,1]],14+j,33+j),30-(17+j),substr(Seq_plus_char_NNHRRT[plus_perfecteditable_NNHRRTdf[i,1]],34+j,39+j),"SaKKH-ABEmax/ABE8e")
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(plus_perfecteditable_NNHRRTdf[1,2])-3))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 3
                          
                        {
                          if(grepl("A", substr(Guide[2],3+l,3+l)) == TRUE)
                          {
                            
                            tmp <- 3+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                      }else{
                        
                        guides <- guides
                      }
                      
                      
                    }
                  }
                }
              
            }
            
          }else{
            
            a <- NULL
            
          }
          
        }
        
 
#---------------------------------------------#    

        # Guide search For SpG-ABEmax (NGC PAM) depending on input'-SNP' and '-Gene Orientation':
        # For Base editing guides on the Minus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "-" || input$SNP == "C>T" && input$`Gene Orientation2` == "+" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "T>C"){
          
          
          if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input" )
          {
            substr(Seq_minus_char_NGC, start = 26, stop = 26) <- "A"
            
            
            
            ## Perfect edible: SpG-ABEmax window 5-6 (PAM: NGC)
            
            PerfectUpstreamPAM_NGC <- substr(Seq_minus_char_NGC, 41,44)  
            
            
            ## Target Plus Strand C>T or Minus Strand G>A
            
            minus_perfecteditable_NGC <- NULL;
            
            for(i in 1:length(PerfectUpstreamPAM_NGC))
            {
              if (grepl( "AGC", PerfectUpstreamPAM_NGC[i]) == TRUE || grepl( "GGC", PerfectUpstreamPAM_NGC[i]) == TRUE || grepl( "CGC", PerfectUpstreamPAM_NGC[i]) == TRUE || grepl( "TGC", PerfectUpstreamPAM_NGC[i]) == TRUE ){    
                tmp <- c(i,PerfectUpstreamPAM_NGC[i])
                minus_perfecteditable_NGC <- rbind(minus_perfecteditable_NGC,tmp)
                
                
              }
            }
            
            if(is.null(minus_perfecteditable_NGC) != TRUE) {
              
              minus_perfecteditabledf_NGC <- data.frame(minus_perfecteditable_NGC)
              
              minus_perfecteditabledf_NGC$X2 <- as.character(minus_perfecteditabledf_NGC$X2)
              
              minus_perfecteditabledf_NGC$X1 <- as.numeric(nrow(minus_perfecteditable_NGC))
              
              
              for (i in 1:nrow(minus_perfecteditabledf_NGC))
                for (j in 1:(nchar(minus_perfecteditabledf_NGC[1,2])-2))
                {
                  {
                    if (grepl( "AGC", substr(minus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(minus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE || grepl( "CGC", substr(minus_perfecteditabledf_NGC[i,2],j,j+2)) == TRUE || grepl( "TGC", substr(minus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE ){
                      
                      #print(minus_perfecteditabledf_NGA[i,2])
                      #print(FANCA_char[perfect_editable_sense[i,1]])
                      Guide <- c(minus_perfecteditabledf_NGC[i,1],substr(Seq_minus_char_NGC[minus_perfecteditabledf_NGC[i,1]],20+j,39+j),26-(19+j),substr(minus_perfecteditabledf_NGC[i,2], j,j+2),"SpG-ABEmax")
                      tmp2 <- NULL
                      
                      Download <- rbind(Download,Guide[2])
                      
                      for (l in 1:(nchar(minus_perfecteditabledf_NGC[1,2])-2))
                        
                        ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                        
                      {
                        if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                        {
                          
                          tmp <- 4+l
                          tmp2 <- rbind(tmp2,tmp)
                          #print("yeha")
                          
                        }else{
                          
                          #print("bla")
                        }
                      }
                      
                      
                      tmp2df <- data.frame(tmp2)
                      tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                      tmp3 <- data.frame(tmp3)
                      tmp4 <- data.frame(tmp3)
                      
                      t <- 0
                      #Guide[6] <- NULL
                      
                      for (b in 1:nrow(tmp3)) {
                        
                        if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                        {
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          t <- t+43
                          tmp4 <- tmp4[] + 43
                          
                        }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                        {
                          
                          Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                          tmp4 <- tmp4[] + 42
                          
                          
                        }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                        {
                          
                          
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          tmp4 <- tmp4[] + 43
                          
                        }
                        
                        
                        
                      }
                      
                      
                      guides <- rbind(guides,Guide)
                      
                      
                      
                      
                    }
                  }
                }
              
            }
          }else{
            
            a <- NULL
            
          }
          
        }
        
        # For Base editing guides on the Plus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "+" || input$SNP == "C>T" && input$`Gene Orientation2` == "-" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "A>G"){
          
          
          if (("G" == substr(Seq_plus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_plus_char_NGC, start = 26, stop = 26) <- "A"
            
            
            
            ## Perfect edible: SpG-ABEmax window 5-6 (PAM: NGC)
            
            PerfectDownstreamPAM_NGC <- substr(Seq_plus_char_NGA, 41,44) 
            
            ## Target Plus Strand G>A or Minus Strand C>T (with NG-ABEmax)
            
            
            plus_perfecteditable_NGC <- NULL;
            
            for(i in 1:length(PerfectDownstreamPAM_NGC))
            {
              if (grepl( "AGC", PerfectDownstreamPAM_NGC[i]) == TRUE || grepl( "GGC", PerfectDownstreamPAM_NGC[i]) == TRUE || grepl( "CGC", PerfectDownstreamPAM_NGC[i]) == TRUE || grepl( "TGC", PerfectDownstreamPAM_NGC[i]) == TRUE ){    
                tmp <- c(i,PerfectDownstreamPAM_NGC[i])
                plus_perfecteditable_NGC <- rbind(plus_perfecteditable_NGC,tmp)
                
                
              }
            }
            
            
            if(is.null(plus_perfecteditable_NGC) != TRUE) {
              
              plus_perfecteditabledf_NGC <- data.frame(plus_perfecteditable_NGC)
              
              plus_perfecteditabledf_NGC$X2 <- as.character(plus_perfecteditabledf_NGC$X2)
              
              plus_perfecteditabledf_NGC$X1 <- as.numeric(nrow(plus_perfecteditable_NGC))
              
              
              for (i in 1:nrow(plus_perfecteditabledf_NGC))
                for (j in 1:(nchar(plus_perfecteditabledf_NGC[1,2])-2))
                {
                  {
                    if (grepl( "AGC", substr(plus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(plus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE || grepl( "CGC", substr(plus_perfecteditabledf_NGC[i,2],j,j+2)) == TRUE || grepl( "TGC", substr(plus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE ){
                      
                      #print(plus_perfecteditabledf_NGA[i,2])
                      #print(FANCA_char[perfect_editable_sense[i,1]])
                      Guide <- c(plus_perfecteditabledf_NGC[i,1],substr(Seq_plus_char_NGC[plus_perfecteditabledf_NGC[i,1]],20+j,39+j),26-(19+j),substr(plus_perfecteditabledf_NGC[i,2], j,j+2),"SpG-ABEmax")
                      tmp2 <- NULL
                      
                      Download <- rbind(Download,Guide[2])
                      
                      for (l in 1:(nchar(plus_perfecteditabledf_NGC[1,2])-2))
                        
                        ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                        
                      {
                        if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                        {
                          
                          tmp <- 4+l
                          tmp2 <- rbind(tmp2,tmp)
                          #print("yeha")
                          
                        }else{
                          
                          #print("bla")
                        }
                      }
                      
                      
                      tmp2df <- data.frame(tmp2)
                      tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                      tmp3 <- data.frame(tmp3)
                      tmp4 <- data.frame(tmp3)
                      
                      t <- 0
                      #Guide[6] <- NULL
                      
                      for (b in 1:nrow(tmp3)) {
                        
                        if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                        {
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          t <- t+43
                          tmp4 <- tmp4[] + 43
                          
                        }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                        {
                          
                          Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                          tmp4 <- tmp4[] + 42
                          
                          
                        }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                        {
                          
                          
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          tmp4 <- tmp4[] + 43
                          
                        }
                        
                        
                        
                      }
                      
                      
                      guides <- rbind(guides,Guide)
                      
                      
                      
                    }
                  }
                }
              
            }
          }else{
            
            a <- NULL
            
          }
          
          
        }
 
#---------------------------------------------#    
        
        # Guide search For SpRY-ABEmax  (NYN PAM) depending on input'-SNP' and '-Gene Orientation':
        # For Base editing guides on the Minus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "-" || input$SNP == "C>T" && input$`Gene Orientation2` == "+" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "T>C"){
          
          
          if (("G" == substr(Seq_minus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input" )
          {
            substr(Seq_minus_char_NYN, start = 26, stop = 26) <- "A"
            
            
            
            ## Perfect edible: SpRY-ABEmax window 5-7 (PAM: NYN)
            
            PerfectUpstreamPAM_NYN <- substr(Seq_minus_char_NYN, 42,44)  
            
            print(PerfectUpstreamPAM_NYN)
            
            
            ## Target Plus Strand C>T or Minus Strand G>A
            
            minus_perfecteditable_NYN <- NULL;
            
            for(i in 1:length(PerfectUpstreamPAM_NYN))
            {
              if (grepl( "C", PerfectUpstreamPAM_NYN[i]) == TRUE || grepl( "T", PerfectUpstreamPAM_NYN[i]) == TRUE || grepl( "A", PerfectUpstreamPAM_NYN[i]) == TRUE){    
                #print(c(i,PerfectUpstreamPAM_NGA[i]))
                tmp <- c(i,PerfectUpstreamPAM_NYN[i])
                minus_perfecteditable_NYN <- rbind(minus_perfecteditable_NYN,tmp)
                
                
              }
            }
            
            if(is.null(minus_perfecteditable_NYN) != TRUE) {
              
              minus_perfecteditabledf_NYN <- data.frame(minus_perfecteditable_NYN)
              
              minus_perfecteditabledf_NYN$X2 <- as.character(minus_perfecteditabledf_NYN$X2)
              
              minus_perfecteditabledf_NYN$X1 <- as.numeric(nrow(minus_perfecteditable_NYN))
              
              
              for (i in 1:nrow(minus_perfecteditabledf_NYN))
                for (j in 1:(nchar(minus_perfecteditabledf_NYN[1,2])-2))
                {
                  {
                    if (grepl( "C", substr(minus_perfecteditabledf_NYN[i,2], j,j)) == TRUE || grepl( "T", substr(minus_perfecteditabledf_NYN[i,2], j,j)) == TRUE || grepl( "A", substr(minus_perfecteditabledf_NYN[i,2], j,j)) == TRUE){
                      
                      Guide <- c(minus_perfecteditabledf_NYN[i,1],substr(Seq_minus_char_NYN[minus_perfecteditabledf_NYN[i,1]],20+j,39+j),26-(19+j),substr(Seq_minus_char_NYN[minus_perfecteditabledf_NYN[i,1]],40+j,42+j),"SpRY-ABEmax (might have lower efficiency!)")
                      
                      tmp2 <- NULL
                      
                      Download <- rbind(Download,Guide[2])
                      
                      for (l in 1:(nchar(minus_perfecteditabledf_NYN[1,2])))
                        
                        ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                        
                      {
                        if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                        {
                          
                          tmp <- 4+l
                          tmp2 <- rbind(tmp2,tmp)
                          #print("yeha")
                          
                        }else{
                          
                          #print("bla")
                        }
                      }
                      
                      
                      tmp2df <- data.frame(tmp2)
                      tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                      tmp3 <- data.frame(tmp3)
                      tmp4 <- data.frame(tmp3)
                      
                      t <- 0
                      #Guide[6] <- NULL
                      
                      for (b in 1:nrow(tmp3)) {
                        
                        if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                        {
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          t <- t+43
                          tmp4 <- tmp4[] + 43
                          
                        }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                        {
                          
                          Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                          tmp4 <- tmp4[] + 42
                          
                          
                        }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                        {
                          
                          
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          tmp4 <- tmp4[] + 43
                          
                        }
                        
                        
                        
                      }
                      
                      
                      guides <- rbind(guides,Guide)
                      
                      
                      
                    }
                  }
                }
              
            }
          }else{
            
            a <- NULL
            
          }
          
        }
        
        # For Base editing guides on the Plus Strand:
        
        if(input$SNP == "G>A" && input$`Gene Orientation2` == "+" || input$SNP == "C>T" && input$`Gene Orientation2` == "-" || input$SequenceInput2 == "Sequence input" && input$Edit3 == "A>G"){
          
          
          if (("G" == substr(Seq_plus_char, start = 26, stop = 26)) == TRUE || input$SequenceInput2 == "Sequence input")
          {
            substr(Seq_plus_char_NYN, start = 26, stop = 26) <- "A"
            
            
            
            ## Perfect edible: SpRY-ABEmax window 5-7 (PAM: NYN)
            
            PerfectDownstreamPAM_NYN <- substr(Seq_plus_char_NYN, 42,44) 
            
            print(PerfectDownstreamPAM_NYN)
            
            ## Target Plus Strand G>A or Minus Strand C>T (with NG-ABEmax)
            
            
            plus_perfecteditable_NYN <- NULL;
            
            for(i in 1:length(PerfectDownstreamPAM_NYN))
            {
              if (grepl( "C", PerfectDownstreamPAM_NGA[i]) == TRUE || grepl( "T", PerfectDownstreamPAM_NGA[i]) == TRUE|| grepl( "A", PerfectDownstreamPAM_NGA[i]) == TRUE){    
                #print(c(i,PerfectDownstreamPAM_NGA[i]))
                tmp <- c(i,PerfectDownstreamPAM_NYN[i])
                plus_perfecteditable_NYN <- rbind(plus_perfecteditable_NYN,tmp)
                
                
              }
            }
            
            
            if(is.null(plus_perfecteditable_NYN) != TRUE) {
              
              plus_perfecteditabledf_NYN <- data.frame(plus_perfecteditable_NYN)
              
              plus_perfecteditabledf_NYN$X2 <- as.character(plus_perfecteditabledf_NYN$X2)
              
              plus_perfecteditabledf_NYN$X1 <- as.numeric(nrow(plus_perfecteditable_NYN))
              
              
              for (i in 1:nrow(plus_perfecteditabledf_NYN))
                for (j in 1:(nchar(plus_perfecteditabledf_NYN[1,2])-2))
                {
                  {
                    if (grepl( "C", substr(plus_perfecteditabledf_NYN[i,2], j,j)) == TRUE || grepl( "T", substr(plus_perfecteditabledf_NYN[i,2], j,j)) == TRUE || grepl( "A", substr(plus_perfecteditabledf_NYN[i,2], j,j)) == TRUE){
                      
                      Guide <- c(plus_perfecteditabledf_NYN[i,1],substr(Seq_plus_char_NYN[plus_perfecteditabledf_NYN[i,1]],20+j,39+j),26-(19+j),substr(Seq_plus_char_NYN[plus_perfecteditabledf_NYN[i,1]],40+j,42+j),"SpRY-ABEmax (might have lower efficiency!)")
                      tmp2 <- NULL
                      
                      Download <- rbind(Download,Guide[2])
                      
                      for (l in 1:(nchar(plus_perfecteditabledf_NYN[1,2])))
                        
                        ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                        
                      {
                        if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                        {
                          
                          tmp <- 4+l
                          tmp2 <- rbind(tmp2,tmp)
                          #print("yeha")
                          
                        }else{
                          
                          #print("bla")
                        }
                      }
                      
                      
                      tmp2df <- data.frame(tmp2)
                      tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                      tmp3 <- data.frame(tmp3)
                      tmp4 <- data.frame(tmp3)
                      
                      t <- 0
                      #Guide[6] <- NULL
                      
                      for (b in 1:nrow(tmp3)) {
                        
                        if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                        {
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          t <- t+43
                          tmp4 <- tmp4[] + 43
                          
                        }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                        {
                          
                          Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                          tmp4 <- tmp4[] + 42
                          
                          
                        }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                        {
                          
                          
                          Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                          tmp4 <- tmp4[] + 43
                          
                        }
                        
                        
                        
                      }
                      
                      
                      guides <- rbind(guides,Guide)
                      
                      
                      
                      
                      
                    }
                  }
                }
              
            }
          }else{
            
            a <- NULL
            
          }
          
          
        }

#---------------------------------------------#
        ## OUTPUT GENERATION:   
        
        # Input control output: Check if object 'a' is NULL (wildtype sequence did not match input SNP) 
        #                         & and no guides have been assigned to the object 'guides':
        
        if((is.null(a) == TRUE) & (is.null(guides) == TRUE)){
          
          output$myText <- renderText({
            
            "Could not find a matching sequence, please check your Substitution!"
          })
        }
   
       
        # Reactive data.frame is created based on the 'guides' object
        # OutputDataTable is generated with the DT package; Change is marked with html tags (red, strong)
          
        if(is.null(guides) != TRUE){
          
          
          currentguidesdf <- reactive({guidesdf <- data.frame("Variant" = input$Variant2, "Protospacer" = guides[,2], "EditPos." = guides[,3], "PAM" = guides[,4], "Base Editor" = guides[,5] ) })
          
          output$mytable  <- DT::renderDataTable({ currentguidesdf() },
          options = list(dom = 't', initComplete =  JS("function(settings,json) {$('body').css({'font-family': 'Helvetica'});",
                                                       "}")),
                                                    escape = FALSE,
                                                    caption = "Output Table",
                                                    callback=JS(
                                                    'table.on("click.dt","td", function() {
                                                    var data=table.cell(this).data();
                                                    if (table.cell(this).data() > 0)
                                                    swal({title: "Edit position", text: "Edit position is defined relative to the PAM sequence. The position in the protospacer closest to the PAM has the position +20, the furthest away position +1. (Figure modified from Gaudelli et.al., 2017)",imageUrl: "Base.jpg",imageSize: "460x200"
                                                    });
                                                 })'
            
            
            
            
            
          ),
          )
          
          
          currentguidesdf2 <- reactive({guidesdf <- data.frame("Variant" = input$Variant2, "Protospacer" = Download , "EditPos." = guides[,3], "PAM" = guides[,4], "Base Editor" = guides[,5] ) })
          
          output$mytable2 <- DT::renderDataTable({ currentguidesdf2()}, caption = "Download Table",options = list(dom = 't'), escape = FALSE)
          
          
        } 
        if((is.null(a) != TRUE) & (is.null(guides) == TRUE)){
          
          output$myText<- renderText({ 
            
            "We couldnt find a Base Editing guide, 
  have you thought about Prime Editing?"
            
          })
        }
     

        # Output control ouput: Check if object 'a' is not NULL (input SNP was correctly assigned),
        #                       but no guides have been assigned to the object 'guides':
        # -> No guides could be designed for this variant, with the Base editors installed.
        
          
      }
        
#-----------------------------------------------------------------------------------------------------------------#    
      
      # Prime editing script:
      # If Prime editing is selected in the Editing strategy this script will be run: 
        
      else if (input$Editing == "Prime editing") {
      
      # Variables for later target definition and quality control are defined: 
  
        a <- 1
        b <- 1
        x <- 0
        z <- 0
        y <- 0
        u <- 0
        m <- 0
        n <- 0
        p <- 0
        v <- 0
        s <- 1
        q <- 0
        
      # Create 'guides', which the later created pegRNAs will be stored to:
        
        guides <- NULL
        
    
#-------------------------------------# 
    
      
        
      # Correcting a patient mutation using Prime Editing:
        
        # Check if there is an input in the 'Mutation' input window, but not in the 'Edit' input window:
        
        if((nchar(input$Mutation) > 2) & (nchar(input$Edit) < 2))
        {
          
          # Variables for later target definition and quality control are defined: 
          
          c <- 1
          d <- 1
          
          if((as.character(input$`Gene Orientation`) == "-") == TRUE) 
          {
            
            
            if(is.na(seqlengths(checkCompatibleSeqinfo(Target_Genome, GRanges(as.character(Chromosome),ranges = IRanges(start = as.numeric(input$Position)+as.numeric(paste(input$`Gene Orientation`,"75", sep = "")), width = 151)))[c(as.character(Chromosome))])) != TRUE)
            {
              
            # Get sequence for Plus Strand using the Bioconductor package
            # Get Sequence for Minus Strand by using reverseComplement function (bioconductor)

            Seq_plus <- getSeq(Target_Genome, GRanges(as.character(Chromosome),ranges = IRanges(start = as.numeric(input$Position)+as.numeric(paste(input$`Gene Orientation`,"75", sep = "")), width = 151)),)
            
            Seq_minus <- reverseComplement(Seq_plus)
            
            # Get reverse Sequences of both Strands using reverse function (bioconductor)
            
            Seq_plus_rev <- reverse(Seq_plus)
            Seq_minus_rev <- reverse(Seq_minus)
            
            # Convert the DnaStringSets into character strings:
            
            Seq_minus_char <- as.character(Seq_minus)
            Seq_plus_char <- as.character(Seq_plus)
            Seq_plus_rev_char <- as.character(Seq_plus_rev)
            Seq_minus_rev_char <- as.character(Seq_minus_rev)
            
            # Create DNA-character Strings for inserting the patient mutation:
            
            Seq_plus_char_n <- Seq_plus_char
            Seq_plus_rev_char_n <- Seq_plus_rev_char
            Seq_minus_char_n <- Seq_minus_char
            Seq_minus_rev_char_n <- Seq_minus_rev_char
            
            # In case of the selection of an insertion, the input sequence is converted into a DNAStringSet and 
            # then included in the DNA character string
            # 'x' and 'u' are assigned the nchar(of input sequence)
            # 'p' is assigned the value 1

            if ((substr(as.character(input$Mutation),1,3) == "ins")==TRUE){
              
              DNAchange <- DNAStringSet(substring(as.character(input$Mutation), first = 4))
              anti_DNAchange <- complement(DNAchange)
              
              DNAchange_char <- paste("[",as.character(DNAchange),"]" ,sep = "")
              
              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- paste("[",as.character(reverse(anti_DNAchange)),"]" ,sep = "")
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              
              x <- nchar(as.character(DNAchange))
              
              u <- nchar(as.character(DNAchange))
              
              x <- as.numeric(x)
              
              p <- 1
              
              Seq_plus_rev_char_n <- paste(substr(Seq_plus_rev_char, start = 1, stop = 75),as.character(anti_DNAchange), substr(Seq_plus_rev_char, start = 76, stop = 151), sep = "")
              Seq_minus_char_n <- paste(substr(Seq_minus_char, start = 1, stop = 75),as.character(DNAchange), substr(Seq_minus_char, start = 76, stop = 151), sep = "")
              
              
            }
            
            # In case of the selection of a deletion, the input sequence is converted into a DNAStringSet
            # 'z' is assigned the nchar(of input sequence)
            # 'x' is assigned the negative value of nchar(of input sequence)
            # 'm' is assigned a value

            else if ((substr(as.character(input$Mutation),1,3) == "del")==TRUE){
              
              DNAchange <- DNAStringSet(substring(as.character(input$Mutation), first = 4))
              anti_DNAchange <- complement(DNAchange)
              
              DNAchange_char <- as.character(DNAchange)
              
              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              z <-  nchar(as.character(DNAchange))
              
              x <- nchar(as.character(DNAchange))
              
              x <- as.numeric(x)* (-1)  
              
              m <- z-1
              
              # Input control check, if sequence in the 'Mutation' input is the same as wildtype sequence at the same string position:
              # If not successful, 'b' is assigned the value NULL

              if(((substr(Seq_minus, start = 76, stop = 75+z)) == (as.character(substring(as.character(input$Mutation), first = 4))))==TRUE){
                
                Seq_plus_rev_char_n <- paste(substr(Seq_plus_rev_char, start = 1, stop = 75), substr(Seq_plus_rev_char, start = 76-x, stop = 151), sep = "")
                Seq_minus_char_n <- paste(substr(Seq_minus_char, start = 1, stop = 75), substr(Seq_minus, start = 76-x, stop = 151), sep = "")
                
              }
              else{
                
                b <- NULL
              }
           
            }
            
            # In case of the selection of a substitution, the sequence before and after '>' is converted into a separated DNAStringSet
            # 'z' is assigned the nchar(of input sequence)
            # 'y' is assigned the value of the length of sequence, which is substituted in the patient
            # 'n' is assigned a value

            else if(grepl(">", as.character(input$Mutation)) == TRUE) { 
              
              y <- ((nchar(as.character(input$Mutation)))-1)/2
              n <- y-1
              
              DNAchange_mut <- DNAStringSet(substr(input$Mutation,2+y,2+(y*2)))
              anti_DNAchange_mut <- complement(DNAchange_mut)
              
              DNAchange <- DNAStringSet(substr(input$Mutation,0,y))
              anti_DNAchange <- complement(DNAchange)
              DNAchange_char <- as.character(DNAchange)
              
              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              # Input control check, if sequence in the 'Mutation' input is the same as wildtype sequence at the same string position:
              # If not successful, 'a' is assigned the value NULL
              
              if ((substr(as.character(input$Mutation),0,y) == substr(Seq_minus, start = 76, stop = 75+y)) == TRUE ){
                
                Seq_plus_rev_char_n <- Seq_plus_rev_char
                Seq_minus_char_n <- Seq_minus_char
                
                substr(Seq_plus_rev_char_n, start = 76, stop = 75+y) <- as.character(anti_DNAchange_mut)
                
                substr(Seq_minus_char_n, start = 76, stop = 75+y) <- as.character(DNAchange_mut)
                
              }
              else{
                
                a <- NULL
                
              }
            }
            
            
            # patient character string is now searched for suitable PAM sequences:

            # Substring for PAM to search for PAM in the editing range of prime editing:
            
            ## Upstream PAM for PRIME:
            # Sense Strand:
            
            primeUpstreamPAM <- substr(Seq_minus_char_n, 5,81)  
            
            # Antisense Strand:
            
            primeUpstreamPAM2 <- substr(Seq_plus_rev_char_n,71,147)

            
            ## Search for PAMs in the Sense Strand:
            
            # 'primeeditable' object is created to store found PAM candidates:

            primeeditable <- NULL;
            
            
            # The PrimeUpstreamPAM substring is searched for all different PAM sequences
            # Matches are stored in the 'primeeditable' object
            
            
            for(i in 1:length(primeUpstreamPAM))
            {
              if (grepl( "AGG", primeUpstreamPAM[i]) == TRUE || grepl( "CGG", primeUpstreamPAM[i]) == TRUE || grepl( "GGG", primeUpstreamPAM[i]) == TRUE || grepl( "TGG", primeUpstreamPAM[i]) == TRUE ){
                
                
                temp <- c(i,primeUpstreamPAM[i])
                primeeditable <- rbind(primeeditable, temp)
                
                
              }
            }
            
            
            primeeditabledf <- data.frame(primeeditable)
            
            primeeditabledf$X2 <- as.character(primeeditabledf$X2)
            
            primeeditabledf$X1 <- as.numeric(primeeditabledf$X1)
            
            
            ## Search for PAMs in the Antisense Strand:
            
            # 'primeeditable' object is created to store found PAM candidates:

            primeeditable2 <- NULL;

            # The PrimeUpstreamPAM substring is searched for all different PAM sequences
            # Matches are stored in the 'primeeditable2' object

            for(i in 1:length(primeUpstreamPAM2))
            {
              if (grepl( "GGA", primeUpstreamPAM2[i]) == TRUE || grepl( "GGC", primeUpstreamPAM2[i]) == TRUE || grepl( "GGG", primeUpstreamPAM2[i]) == TRUE || grepl( "GGT", primeUpstreamPAM2[i]) == TRUE ){
                
                
                temp <- c(i,primeUpstreamPAM2[i])
                primeeditable2 <- rbind(primeeditable2, temp)
                
                
              }
            }
            
            
            primeeditable2df <- data.frame(primeeditable2)
            
            primeeditable2df$X2 <- as.character(primeeditable2df$X2)
            
            primeeditable2df$X1 <- as.numeric(primeeditable2df$X1)
            
            
            # Guides are only designed, if deleted sequences in the patient match wildtype seq. (b <- 1) 
            # or substituted sequences in the patient match wildtype seq. (a <- 1)
            
            
            if(is.null(a) != TRUE & is.null(b) != TRUE){
              
            # Guide sequences are designed based on the position of the PAM sequence in the 'primeeditable' object:

              for (i in 1:nrow(primeeditabledf))
                for (j in 1:(nchar(primeeditabledf[1,2])-2))
                {
                  {
                    if (grepl( "AGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(primeeditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE ){
                      
                      # Check if Editing Position is included in the RT template:
                      
                      if((((76-j)+abs(x)) > as.numeric(input$RT)) != TRUE)
                      {
                        
                      # Check if Protospacer length is still 20 nt:

                      if(nchar(substr(Seq_minus_char_n[primeeditabledf[i,1]],(j-16),j+3)) == 20)
                      {
                        
                        Guide <- c("1",substr(Seq_minus_char_n[primeeditabledf[i,1]],j-16,j+3),76-j,as.character(reverse(DNAStringSet(substr(Seq_plus_rev_char[primeeditabledf[i,1]],(j+1) - as.numeric(input$PBS),((j+1) - (as.numeric(input$PBS)))+ (as.numeric(input$RT)+as.numeric(input$PBS)-1))))),substr(primeeditabledf[i,2], j,j+2), "Sense")
                        
                        
                        
                        # Modify Protospacer to enhance transcription: (Add 3'-G to Protospacer if there is no G)
                        
                        if(substr(Seq_minus_char_n[primeeditabledf[i,1]],j-16,j-16) != "G")
                        {
                          
                          Guide[2] <- paste("g",Guide[2], sep = "") 
                          
                        }
                        
                        Guide[7] <- 0
                        
                        if(substr(Guide[4],1,1) == "C")
                        {
                          
                          Guide[7] <- as.numeric(Guide[7])-28
                          
                        }
                        
                        if ((grepl("TTTTT", Guide[4])) == TRUE)
                        {                          
                          
                          Guide[7] <- as.numeric(Guide[7])-50
                          
                        }
                        
                        if (((((as.numeric(input$RT))) - ((as.numeric(Guide[3])))) > 4) != TRUE) 
                        {
                          
                        Guide[7] <- as.numeric(Guide[7])-6
                          
                        }
   

                        Guide[8] <- paste(substring(Guide[4],1,(-1)+as.numeric(input$RT)-((as.numeric(Guide[3])-1+m+n))),as.character(anti_DNAchange_char_tagged), substring(Guide[4], first = ((as.numeric(input$RT))-(as.numeric(Guide[3])-2+p))), sep = "")
                        
                        
                        Guide[7] <- as.numeric(Guide[7]) - as.numeric(Guide[3])
                        
                        guides <- rbind(guides,Guide)
                        
                      
                        }
                        }
                      
                    }
                  }
                }
              
              if(is.null(a) != TRUE & is.null(b) != TRUE){
                
                ## AntiSense strand
                
                for (i in 1:nrow(primeeditable2df))
                  for (j in 1:(nchar(primeeditable2df[1,2])-2))
                  {
                    {
                      if (grepl( "GGA", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(primeeditable2df[i,2],j,j+2)) == TRUE || grepl( "GGT", substr(primeeditable2df[i,2], j,j+2)) == TRUE ){
                        
                        if(((j+abs(x)) > as.numeric(input$RT)) != TRUE)
                        {
                          if(nchar(substr(Seq_plus_rev_char_n[primeeditabledf[i,1]],j+73,j+92)) == 20)
                          {
                          
                          Guide <- c("1",substr(Seq_plus_rev_char_n[primeeditable2df[i,1]],j+73,j+92),j-u, as.character(reverse(DNAStringSet(substr(Seq_minus_char[primeeditable2df[i,1]],((j+76-x)+(as.numeric(input$PBS)))-(as.numeric(input$RT)+as.numeric(input$PBS)),((j+76-x)+(as.numeric(input$PBS)))-1)))),as.character(reverse(DNAStringSet(substr(primeeditable2df[i,2], j,j+2)))), "Antisense")

                          
                          if(substr(Seq_plus_rev_char_n[primeeditable2df[i,1]],j+92,j+92) != "G")
                          {
                            
                            Guide[2] <- paste(Guide[2],"g", sep = "") 
                            
                          }
                          
                          Guide[2] <- as.character(reverse(DNAStringSet(Guide[2])))
                          
                          Guide[4] <- as.character(reverse(DNAStringSet(Guide[4])))
                          
                          Guide[7] <- 0
                          
                          if((substr(Guide[4],1,1) == "C") == TRUE){
                            
                            Guide[7] <- as.numeric(Guide[7])-28
                            
                          }
                          if ((grepl("TTTTT", Guide[4])) == TRUE){
                            
                            Guide[7] <- as.numeric(Guide[7])-50
                            
                          }
                          
                          if (((((as.numeric(input$RT))) - ((as.numeric(Guide[3])))) > 4) != TRUE) 
                          {
                            
                            Guide[7] <- as.numeric(Guide[7])-6
                            
                          }
                          
                          Guide[8] <- paste(substring(Guide[4],1,as.numeric(input$RT)-((as.numeric(Guide[3])+z-p))),as.character(DNAchange_char_tagged), substring(Guide[4], first = 1+((as.numeric(input$RT))-(as.numeric(Guide[3])-y-p))), sep = "")
                          
                          Guide[7] <- as.numeric(Guide[7]) - as.numeric(Guide[3])
                          
                          guides <- rbind(guides,Guide)
                          
                          }
                        }
                      }
                      
                    }
                  }
                
                
              }
              
              
              
              if (is.null(guides)!= TRUE)
              {
                
                for (m in 1:nrow(guides)) {
                  if((as.character(guides[m,6]) == "Sense") == TRUE){
                    
                    DNAchange_char_tagged <- anti_DNAchange_char_tagged
                    
                  }else{}
                  
                  
                }
                
                
                if(any(((as.numeric(guides[,3])+u-z) <= input$RT) || (as.numeric(guides[,3])+u-z >= 0)))
                {
                  guides <- subset(guides, as.numeric(guides[,3])+u-z <= input$RT)
                  
                }else{
                  
                  d <- NULL
                  c <- NULL
                  
                }
                
                if(is.null(c) || is.null(d) != TRUE)
                {
                  
                  
                  
                  currentguidesdf <- reactive({ 
                    
                    
                    guidesdf <- data.frame("Variant" = input$Variant, "Score" = as.numeric(guides[,7]), "Protospacer(Sense)" = guides[,2],"Protospacer(Antisense)" = as.character(reverseComplement(DNAStringSet(guides[,2]))) , "EditPos." = guides[,3], "Extension(coding_strand)" = guides[,8], "PAM" = guides[,5], "PAM-Strand" = guides[,6], stringsAsFactors = FALSE)
                    guidesdf[order(guidesdf$Score,decreasing = TRUE), ]
                    
                  })
                  
                }
                
                
              }
              
              
            }
            
            }else{
              
              output$myText <- renderText({
                
                "Could not find a matching sequence, please check your genomic coordinates!"
              })
              
              guides <- 1
            }
            
            
            }
            
          
          if((as.character(input$`Gene Orientation`) == "+") == TRUE)
          {
            
            if(is.na(seqlengths(checkCompatibleSeqinfo(Target_Genome, GRanges(as.character(Chromosome),ranges = IRanges(start = as.numeric(input$Position)+as.numeric(paste(input$`Gene Orientation`,"75", sep = "")), width = 151)))[c(as.character(Chromosome))])) != TRUE)
            {
            ### Get sequence for Plus Strand
            
            Seq_plus <- getSeq(Target_Genome, GRanges(as.character(Chromosome),ranges = IRanges(start = as.numeric(input$Position)-as.numeric(paste(input$`Gene Orientation`,"75", sep = "")), width = 151)),)
            
            ### Get Sequence for Minus Strand
            
            Seq_minus <- reverseComplement(Seq_plus)
            
            Seq_plus_rev <- reverse(Seq_plus)
            Seq_minus_rev <- reverse(Seq_minus)
            
            Seq_minus_char <- as.character(Seq_minus)
            Seq_plus_char <- as.character(Seq_plus)
            Seq_plus_rev_char <- as.character(Seq_plus_rev)
            Seq_minus_rev_char <- as.character(Seq_minus_rev)
            
            Seq_plus_char_n <- Seq_plus_char
            Seq_plus_rev_char_n <- Seq_plus_rev_char
            
            Seq_minus_char_n <- Seq_minus_char
            Seq_minus_rev_char_n <- Seq_minus_rev_char
            
            
            
            if ((substr(input$Mutation,1,3) == "ins")==TRUE){
              
              DNAchange <- DNAStringSet(substring(input$Mutation, first = 4))
              anti_DNAchange <- complement(DNAchange)
              
              
              DNAchange_char <- paste("[",as.character(DNAchange),"]" ,sep = "")

              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- paste("[",as.character(reverse(anti_DNAchange)),"]" ,sep = "")
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
             
              
              x <- nchar(as.character(DNAchange))
              
              u <- nchar(as.character(DNAchange))
              
              x <- as.numeric(x)
              
              Seq_plus_char_n <- paste(substr(Seq_plus_char, start = 1, stop = 75),as.character(DNAchange), substr(Seq_plus_char, start = 76, stop = 151), sep = "")
              Seq_minus_rev_char_n <- paste(substr(Seq_minus_rev_char, start = 1, stop = 75),as.character(anti_DNAchange), substr(Seq_minus_rev_char, start = 76, stop = 151), sep = "")
              
            }else if ((substr(input$Mutation,1,3) == "del")==TRUE){
              
              DNAchange <- DNAStringSet(substring(input$Mutation, first = 4))
              anti_DNAchange <- complement(DNAchange)
           
              DNAchange_char <- as.character(DNAchange)
              
              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              
              
              z <- nchar(as.character(DNAchange))
              
              
              x <- nchar(as.character(DNAchange))
              
              x <- as.numeric(x)* (-1)  
              
              
              if(((substr(Seq_plus, start = 76, stop = 75+z)) == (as.character(DNAchange)))){
                
                Seq_plus_char_n <- paste(substr(Seq_plus_char, start = 1, stop = 75), substr(Seq_plus_char, start = 76-x, stop = 151), sep = "")
                
                Seq_minus_rev_char_n <- paste(substr(Seq_minus_rev_char, start = 1, stop = 75), substr(Seq_minus_rev_char, start = 76-x, stop = 151), sep = "")
                
                
              }else{
                
                b <- NULL
                
              }
              
              
            }else if(grepl(">", input$Mutation) == TRUE) { 
              
              
              x <- 0
              y <- 0
              y <- ((nchar(input$Mutation))-1)/2
              
              
              DNAchange_mut <- DNAStringSet(substr(input$Mutation,2+y,2+(y*2)))
              anti_DNAchange_mut <- complement(DNAchange_mut)
              
              DNAchange <- DNAStringSet(substr(input$Mutation,0,y))
              anti_DNAchange <- complement(DNAchange)
              
              DNAchange_char <- as.character(DNAchange)
              
              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- as.character(anti_DNAchange)
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", reverse(anti_DNAchange_char),.noWS = "outside"),.noWS = "outside")
              
              
              if ((substr(input$Mutation,0,y) == substr(Seq_plus, start = 76, stop = 75+y))){
                
                Seq_plus_char_n <- Seq_plus_char
                Seq_minus_rev_char_n <- Seq_minus_rev_char
                
                substr(Seq_plus_char_n, start = 76, stop = 75+y) <- as.character(DNAchange_mut)
                
                substr(Seq_minus_rev_char_n, start = 76, stop = 75+y) <- as.character(anti_DNAchange_mut)
                
              }else{
                
                a <- NULL
                
              }
            }
            
            
            ## Searchig for sense PAM in Prime editing 
            
            
            
            
            ## Upstream PAM for PRIME 
            
            primeUpstreamPAM <- substr(Seq_plus_char_n, 5,81)  
            
            # smaller search area, since large deletions can otherwise not be searched for (large x):
            
            primeUpstreamPAM2 <- substr(Seq_minus_rev_char_n,71+x+y,110+x+y)
            
            
            
            
            ## prime edits possible:
            
            
            ## Sense strand
            
            primeeditable <- NULL;
            
            
            for(i in 1:length(primeUpstreamPAM))
            {
              if (grepl( "AGG", primeUpstreamPAM[i]) == TRUE || grepl( "CGG", primeUpstreamPAM[i]) == TRUE || grepl( "GGG", primeUpstreamPAM[i]) == TRUE || grepl( "TGG", primeUpstreamPAM[i]) == TRUE ){
                
                
                temp <- c(i,primeUpstreamPAM[i])
                primeeditable <- rbind(primeeditable, temp)
                
                
              }
            }
            
            
            primeeditabledf <- data.frame(primeeditable)
            
            primeeditabledf$X2 <- as.character(primeeditabledf$X2)
            
            primeeditabledf$X1 <- as.numeric(primeeditabledf$X1)
            
            
            
            
            ## AntiSense strand
            
            primeeditable2 <- NULL;
            
            
            for(i in 1:length(primeUpstreamPAM2))
            {
              if (grepl( "GGA", primeUpstreamPAM2[i]) == TRUE || grepl( "GGC", primeUpstreamPAM2[i]) == TRUE || grepl( "GGG", primeUpstreamPAM2[i]) == TRUE || grepl( "GGT", primeUpstreamPAM2[i]) == TRUE ){
                
                
                temp <- c(i,primeUpstreamPAM2[i])
                primeeditable2 <- rbind(primeeditable2, temp)
                
                
              }
            }
            
            
            primeeditable2df <- data.frame(primeeditable2)
            
            primeeditable2df$X2 <- as.character(primeeditable2df$X2)
            
            primeeditable2df$X1 <- as.numeric(primeeditable2df$X1)
            
            
            if(is.null(a) != TRUE & is.null(b) != TRUE){
              
            if(is.null(primeeditable) != TRUE){
              
              ## Sense strand
              
              for (i in 1:nrow(primeeditabledf))
                for (j in 1:(nchar(primeeditabledf[1,2])-2))
                {
                  {
                    if (grepl( "AGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(primeeditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE ){
                      
                      
                      if((nchar(substr(Seq_plus_char_n[primeeditabledf[i,1]],(j-16),j+3)) == 20))
                      {
                        
                        Guide <- c("1",substr(Seq_plus_char_n[primeeditabledf[i,1]],j-16,j+3),76-j,as.character(reverse(DNAStringSet(substr(Seq_minus_rev_char[primeeditabledf[i,1]],(j+1) - as.numeric(input$PBS),((j+1) - (as.numeric(input$PBS)))+ (as.numeric(input$RT)+as.numeric(input$PBS)-1))))), substr(primeeditabledf[i,2], j,j+2), "Sense")
                        
                        if(substr(Seq_plus_char_n[primeeditabledf[i,1]],j-16,j-16) != "G")
                        {
                          
                          Guide[2] <- paste("g",Guide[2], sep = "") 
                          
                        }
                        
                        Guide[7] <- 0
                        
                        if((substr(Guide[4],1,1) == "C") == TRUE){
                          
                          Guide[7] <- as.numeric(Guide[7])-28
                          
                        }
                        if ((grepl("TTTTT", Guide[4])) == TRUE){
                          
                          Guide[7] <- as.numeric(Guide[7])-50
                          
                        }
                        
                        if (((((as.numeric(input$RT)-z+1)) - ((as.numeric(Guide[3])))) > 4) != TRUE) 
                        {
                          
                          Guide[7] <- as.numeric(Guide[7])-6 
                          
                        }
                        
                        
                        
                        Guide[7] <- as.numeric(Guide[7])-as.numeric(Guide[3])
                        
                        Guide[8] <- paste(substring(Guide[4],1,1+as.numeric(input$RT)-((as.numeric(Guide[3])+(z)+y))),as.character(anti_DNAchange_char_tagged), substring(Guide[4], first = 2+((as.numeric(input$RT))-(as.numeric(Guide[3])))), sep = "")
                        
                        guides <- rbind(guides,Guide)
                        
                      }
                      

                      
                    } 
                  }
                }
              }
            }
            
            if((is.null(a) != TRUE) & (is.null(b) != TRUE)){
              
              if(is.null(primeeditable2) != TRUE){
              
              ## AntiSense strand
              
              for (i in 1:nrow(primeeditable2df))
                for (j in 1:(nchar(primeeditable2df[1,2])-2))
                {
                  {
                    if (grepl( "GGA", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(primeeditable2df[i,2],j,j+2)) == TRUE || grepl( "GGT", substr(primeeditable2df[i,2], j,j+2)) == TRUE ){
                      
                     
                      if(nchar(substr(Seq_minus_rev_char_n[primeeditable2df[i,1]],j+73+x+y,j+92+x+y)) == 20){
                        
                        
                        Guide <- c("1",substr(Seq_minus_rev_char_n[primeeditable2df[i,1]],j+73+x+y,j+92+x+y),j+y-z+1, as.character(reverse(DNAStringSet(substr(Seq_plus_char[primeeditable2df[i,1]],((j+76)+(as.numeric(input$PBS)))-(as.numeric(input$RT)+as.numeric(input$PBS)),((j+76)+(as.numeric(input$PBS)))-1)))),as.character(reverse(DNAStringSet(substr(primeeditable2df[i,2], j,j+2)))), "Antisense")
                        
                        
                        if(substr(Seq_minus_rev_char_n[primeeditable2df[i,1]],j+73+x+y,j+73+x+y) != "G")
                        {
                          
                          Guide[2] <- paste(Guide[2],"g", sep = "") 
                          
                        }
                        
                        Guide[2] <- as.character(reverse(DNAStringSet(Guide[2])))
                        
                        Guide[4] <- as.character(reverse(DNAStringSet(Guide[4])))
                        
                        Guide[7] <- 0
                        
                        if((substr(Guide[4],1,1) == "C") == TRUE){
                          
                          Guide[7] <- as.numeric(Guide[7])-28
                          
                        }
                        
                        if ((grepl("TTTTT", Guide[4])) == TRUE){
                          
                          Guide[7] <- as.numeric(Guide[7])-50
                          
                        }
                        
                        if (((((as.numeric(input$RT)-z)) - ((as.numeric(Guide[3])))) >= 4) != TRUE) 
                        {
                          
                          Guide[7] <- as.numeric(Guide[7])-6 
                          
                        }
                        

                        Guide[8] <- paste(substring(Guide[4],1,1+as.numeric(input$RT)-((as.numeric(Guide[3])+z))),as.character(DNAchange_char_tagged), substring(Guide[4], first = 2+((as.numeric(input$RT))-(as.numeric(Guide[3])-y))), sep = "")
                        
                        
                        Guide[7] <- as.numeric(Guide[7])-as.numeric(Guide[3])
                        
                        
                        guides <- rbind(guides,Guide)
                        
                      }
                      
                      
                    }
                  }
                }
              
              
              }
            }
            
            
            if (is.null(guides)!= TRUE)
            {
              
              
              if(any(as.numeric(guides[,3])+u-z <= input$RT) == TRUE)
              {
                guides <- subset(guides, as.numeric(guides[,3])+u-z <= input$RT)
                
              }else{
                
                d <- NULL
                c <- NULL
                
              }
              
              if(is.null(c) | is.null(d) != TRUE)
              {
                
                
                
                currentguidesdf <- reactive({ 
                  
                 
                  guidesdf <- data.frame("Variant" = input$Variant, "Score" = as.numeric(guides[,7]), "Protospacer(Sense)" = guides[,2],"Protospacer(Antisense)" = as.character(reverseComplement(DNAStringSet(guides[,2]))) , "EditPos." = guides[,3], "Extension(coding_strand)" = guides[,8], "PAM" = guides[,5], "PAM-Strand" = guides[,6], stringsAsFactors = FALSE)
                  guidesdf[order(guidesdf$Score,decreasing = TRUE), ]
                  
                })
                
              }
            }

            }else{
              
              output$myText <- renderText({
                
                "Could not find a matching sequence, please check your genomic coordinates!"
                
              })
              
              guides <- 1
            }
            
          }
         
          
          ### Make Oligo table, for ready to clone oligos 
          
          
          currentguidesdf2 <- reactive({ 
            
            
            guidesdf <- data.frame("Variant" = input$Variant, "Score" = as.numeric(guides[,7]), "Protospacer(Sense)" = paste("cacc", guides[,2],"gtttt", sep = ""),"Protospacer(Antisense)" = paste("ctctaaaac",reverseComplement(DNAStringSet(guides[,2])),sep = "") , "TargetPos." = guides[,3], "Extension(Sense)" = paste("gtgc",guides[,4],sep = ""), "Extension(Antisense)" = paste("aaaa",reverseComplement(DNAStringSet(guides[,4])), sep = ""), "PAM" = guides[,5], "PAM-Strand" = guides[,6], stringsAsFactors = FALSE)
            guidesdf[order(guidesdf$Score,decreasing = TRUE), ]
            
            
          }) 
          
          
          
          
          if(is.null(a)== TRUE){
            
            output$myText <- renderText({
              
              "Could not find a matching sequence, please check your Substitution!"
            })
          }
          else if(is.null(b)== TRUE){
            
            output$myText <- renderText({
              
              "Could not find a matching sequence, please check your Deletion!"
            })
            
          }
          else if(is.null(c) & is.null(d) == TRUE){
            
            output$myText <- renderText({
              
              "Edit is to far away, try to increase the RT length"
            })
          }
          else if(is.null(guides)== TRUE){
            
            output$myText <- renderText({
              
              "We couldn't find a Prime editing Guide!"
            })
            
          }else{
            
            output$mytable  <- DT::renderDataTable({ currentguidesdf()},
                                                   escape = FALSE,
                                                   caption = "Output Table" ,
                                                   options = list(dom = 't', initComplete =  JS("function(settings, json) {",
                                                                                                "$('body').css({'font-family': 'Helvetica'});",
                                                                                                "}")),
                                                   selection = list(mode = 'single', target = 'column'),
                                                   callback=JS(
                                                     'table.on("click.dt","td", function() {
                                                    var data=table.cell(this).data();
                                                    if (table.cell(this).data() < 0)
                                                    swal({title: "pegRNA-Score", text: "The higher the better! (calculated based on recommendation from the Liu Lab)",imageUrl: "Score.jpg",imageSize: "460x200"
                                                    });   
                                                    if (table.cell(this).data() > 0)
                                                    swal({title: "Edit position", text:  "Edit position is defined relative to the PAM sequence (Anzalone et.al. 2019)",imageUrl: "Target.jpg",imageSize: "460x200"
                                                    });                                                    
                                                    })')
                                                   
                                                   
            ) 
            
         
            output$mytable2 <- DT::renderDataTable({ currentguidesdf2()}, caption = "Ordering Table",options = list(dom = 't'))
          }
          
        }
       
#-------------------------------------# 
        
      # Installing an edit using Prime Editing:
        
        # Check if there is an input in the 'Edit' input window, but not in the 'Mutation' input window:
 
      else if ((nchar(input$Edit) > 2) & (nchar(input$Mutation) < 2) || input$SequenceInput == "Sequence input")
        {
      
        
        if((as.character(input$`Gene Orientation`) == "-") == TRUE)
        {
          c <- 1
          d <- 1
          
          if(is.na(seqlengths(checkCompatibleSeqinfo(Target_Genome, GRanges(as.character(Chromosome),ranges = IRanges(start = as.numeric(input$Position)+as.numeric(paste(input$`Gene Orientation`,"75", sep = "")), width = 151)))[c(as.character(Chromosome))])) != TRUE)
          {
          
          ### Get sequence for Plus Strand
          if ((substr(as.character(input$Edit),1,3) == "del")==TRUE){
            
            
          m <- nchar(substring(as.character(input$Edit), first = 4))-1
            
            
          Seq_plus <- getSeq(Target_Genome, GRanges(as.character(Chromosome),ranges = IRanges(start = (as.numeric(input$Position)+m)+as.numeric(paste(input$`Gene Orientation`,"75", sep = "")), width = 151)),)
          
          }else{
            
          Seq_plus <- getSeq(Target_Genome, GRanges(as.character(Chromosome),ranges = IRanges(start = as.numeric(input$Position)+as.numeric(paste(input$`Gene Orientation`,"75", sep = "")), width = 151)),)
            
            
          }
          ### Get Sequence for Minus Strand
          
          Seq_minus <- reverseComplement(Seq_plus)
          
          Seq_plus_rev <- reverse(Seq_plus)
          Seq_minus_rev <- reverse(Seq_minus)
          
          Seq_minus_char <- as.character(Seq_minus)
          Seq_plus_char <- as.character(Seq_plus)
          Seq_plus_rev_char <- as.character(Seq_plus_rev)
          Seq_minus_rev_char <- as.character(Seq_minus_rev)
          
          
          if ((substr(as.character(input$Edit),1,3) == "ins")==TRUE){
            
            DNAchange <- DNAStringSet(substring(as.character(input$Edit), first = 4))
            anti_DNAchange <- complement(DNAchange)
            
            DNAchange_char <- as.character(reverse(DNAchange))
            
            DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            anti_DNAchange_char <- as.character(anti_DNAchange)
            
            anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            p <- 1
            
            x <- nchar(as.character(DNAchange))
            u <- nchar(as.character(DNAchange))
            
            x <- as.numeric(x)
            
            Seq_plus_rev_char_n <- paste(substr(Seq_plus_rev_char, start = 1, stop = 75),as.character(anti_DNAchange), substr(Seq_plus_rev_char, start = 76, stop = 151), sep = "")
            Seq_minus_char_n <- paste(substr(Seq_minus_char, start = 1, stop = 75),as.character(DNAchange), substr(Seq_minus_char, start = 76, stop = 151), sep = "")
            
            
          }
          else if ((substr(as.character(input$Edit),1,3) == "del")==TRUE){
            
            DNAchange <- DNAStringSet(substring(as.character(input$Edit), first = 4))
            anti_DNAchange <- complement(DNAchange)
            
            DNAchange_char <- paste("[",as.character(DNAchange),"]" ,sep = "")
            
            DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            anti_DNAchange_char <- paste("[",as.character(reverse(anti_DNAchange)),"]" ,sep = "")
            
            anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            z <-  nchar(as.character(DNAchange))
            
            x <- nchar(as.character(DNAchange))
            
            x <- as.numeric(x)* (-1)  
            
            q <- 1
            
            if(((substr(Seq_minus_char, start = 76, stop = 75+z)) == (as.character(substring(as.character(input$Edit), first = 4))))==TRUE){
              
              Seq_plus_rev_char_n <- paste(substr(Seq_plus_rev_char, start = 1, stop = 75), substr(Seq_plus_rev_char, start = 76+z, stop = 151), sep = "")
              Seq_minus_char_n <- paste(substr(Seq_minus_char, start = 1, stop = 75), substr(Seq_minus, start = 76+z, stop = 151), sep = "")
              
            }else{
       
              b <- NULL
            }
            
            
          }
          else if(grepl(">", as.character(input$Edit)) == TRUE) { 
            
            
            y <- 0
            y <- ((nchar(as.character(input$Edit)))-1)/2
            
            
            DNAchange <- DNAStringSet(substr(as.character(input$Edit),2+y,2+(y*2)))
            anti_DNAchange <- complement(DNAchange)
            
            DNAchange_char <- as.character(DNAchange)
            
            DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            anti_DNAchange_char <- as.character(anti_DNAchange)
            
            anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            
            if ((substr(as.character(input$Edit),0,y) == substr(Seq_minus, start = 76, stop = 75+y)) == TRUE ){
              
              Seq_plus_rev_char_n <- Seq_plus_rev_char
              Seq_minus_char_n <- Seq_minus_char
              
              substr(Seq_plus_rev_char_n, start = 76, stop = 75+y) <- as.character(anti_DNAchange)
              
              substr(Seq_minus_char_n, start = 76, stop = 75+y) <- as.character(DNAchange)
              
            }else{
              
              a <- NULL
              
            }
          }
          
          
          ## Searchig for sense PAM in Prime editing 
          
          
          
          
          ## Upstream PAM for PRIME 
          
          primeUpstreamPAM <- substr(Seq_minus_char, 5,81)  
          
          primeUpstreamPAM2 <- substr(Seq_plus_rev_char,71,147)
          
          
          ## prime edits possible:
          
          
          ## Sense strand
          
          primeeditable <- NULL;
          
          
          for(i in 1:length(primeUpstreamPAM))
          {
            if (grepl( "AGG", primeUpstreamPAM[i]) == TRUE || grepl( "CGG", primeUpstreamPAM[i]) == TRUE || grepl( "GGG", primeUpstreamPAM[i]) == TRUE || grepl( "TGG", primeUpstreamPAM[i]) == TRUE ){
              
              
              temp <- c(i,primeUpstreamPAM[i])
              primeeditable <- rbind(primeeditable, temp)
              
              
            }
          }
          
          
          primeeditabledf <- data.frame(primeeditable)
          
          primeeditabledf$X2 <- as.character(primeeditabledf$X2)
          
          primeeditabledf$X1 <- as.numeric(primeeditabledf$X1)
          
          
          
          
          ## AntiSense strand
          
          primeeditable2 <- NULL;
          
          
          for(i in 1:length(primeUpstreamPAM2))
          {
            if (grepl( "GGA", primeUpstreamPAM2[i]) == TRUE || grepl( "GGC", primeUpstreamPAM2[i]) == TRUE || grepl( "GGG", primeUpstreamPAM2[i]) == TRUE || grepl( "GGT", primeUpstreamPAM2[i]) == TRUE ){
              
              
              temp <- c(i,primeUpstreamPAM2[i])
              primeeditable2 <- rbind(primeeditable2, temp)
              
              
            }
          }
          
          
          primeeditable2df <- data.frame(primeeditable2)
          
          primeeditable2df$X2 <- as.character(primeeditable2df$X2)
          
          primeeditable2df$X1 <- as.numeric(primeeditable2df$X1)
          
          
          
          
          if(is.null(a) != TRUE && is.null(b) != TRUE){
            
            if(is.null(primeeditable) != TRUE){
            
            ## Sense strand
            
            
            for (i in 1:nrow(primeeditabledf))
              for (j in 1:(nchar(primeeditabledf[1,2])-2))
              {
                {
                  if (grepl( "AGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(primeeditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE ){
                    

                    
                    if(nchar(substr(Seq_minus_char[primeeditabledf[i,1]],(j-16),j+3)) == 20)
                    {
                      
                      Guide <- c("1",substr(Seq_minus_char[primeeditabledf[i,1]],j-16,j+3),76-j,as.character(reverse(DNAStringSet(substr(Seq_plus_rev_char_n[primeeditabledf[i,1]],(j+1) - as.numeric(input$PBS),((j+1) - (as.numeric(input$PBS)))+ (as.numeric(input$RT)+as.numeric(input$PBS)-1))))),substr(primeeditabledf[i,2], j,j+2), "Sense")
                 

                      if(substr(Seq_minus_char[primeeditabledf[i,1]],j-16,j-16) != "G")
                      {
                        
                        Guide[2] <- paste("g",Guide[2], sep = "") 
                        
                      }
                      
                      Guide[7] <- 0
                      
                      if((substr(Guide[4],1,1) == "C") == TRUE){
                        
                        Guide[7] <- as.numeric(Guide[7])-28
                        
                      }
                      
                      if ((grepl("TTTTT", Guide[4])) == TRUE){
                        
                        Guide[7] <- as.numeric(Guide[7])-50
                        
                      }
                      
                      if (((as.numeric(input$RT) - (as.numeric(Guide[3])+y-x)) >= 4) != TRUE) 
                      {
                        
                        Guide[7] <- as.numeric(Guide[7])-6
                        
                      }
                      
                      
                      Guide[7] <- as.numeric(Guide[7]) - as.numeric(Guide[3])
                      
                      ## Add Guide[8] with the html tagged edit in red and bold on the coding strand
                      
                      Guide[8] <- paste(reverseComplement(DNAStringSet(substring(Guide[4], first = 1+z+q+y+p+((as.numeric(input$RT)-(as.numeric(Guide[3]))))))),as.character(DNAchange_char_tagged),reverseComplement(DNAStringSet(substring(Guide[4],1,as.numeric(input$RT)-((as.numeric(Guide[3])-z+u-q-p))))), sep = "")

                      guides <- rbind(guides,Guide)
                      
                        
                    }
                    
                  }
                }
              }
            }
          }
          
          if(is.null(a) != TRUE && is.null(b) != TRUE){
            
            if(is.null(primeeditable2) != TRUE){
            
            ## AntiSense strand
            
            for (i in 1:nrow(primeeditable2df))
              for (j in 1:(nchar(primeeditable2df[1,2])-2))
              {
                {
                  if (grepl( "GGA", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(primeeditable2df[i,2],j,j+2)) == TRUE || grepl( "GGT", substr(primeeditable2df[i,2], j,j+2)) == TRUE ){
            
                    if(nchar(substr(Seq_plus_rev_char[primeeditabledf[i,1]],(j+73),j+92)) == 20)
                    {
                    
                    Guide <- c("1",substr(Seq_plus_rev_char[primeeditable2df[i,1]],j+73,j+92),j, as.character(reverse(DNAStringSet(substr(Seq_minus_char_n[primeeditable2df[i,1]],((j+76)+(as.numeric(input$PBS)))-(as.numeric(input$RT)+as.numeric(input$PBS)),((j+76)+(as.numeric(input$PBS)))-1)))),as.character(reverse(DNAStringSet(substr(primeeditable2df[i,2], j,j+2)))), "Antisense")
                    
                    #Guide <- c(row.names(input),substr(Seq_minus_rev_char[primeeditable2df[i,1]],j+33,j+52),j, as.character(reverse(DNAStringSet(substr(Seq_plus_char_n[primeeditable2df[i,1]],(j+37)-as.numeric(input$PBS),((j+37)-(as.numeric(input$PBS)))+(as.numeric(input$RT)+as.numeric(input$PBS)-1))))), "NGG", "Antisense")
                    
                  
                    
                    if(substr(Seq_plus_rev_char[primeeditable2df[i,1]],j+92,j+92) != "G")
                    {
                      
                      Guide[2] <- paste(Guide[2],"g", sep = "") 
                      
                    }
                    
                    Guide[2] <- as.character(reverse(DNAStringSet(Guide[2])))
                    
                    Guide[4] <- as.character(reverse(DNAStringSet(Guide[4])))
                    
                    Guide[7] <- 0
                    
                    if((substr(Guide[4],1,1) == "C") == TRUE){
                      
                      Guide[7] <- as.numeric(Guide[7])-28
                      
                    }
                    
                    if ((grepl("TTTTT", Guide[4])) == TRUE){
                      
                      Guide[7] <- as.numeric(Guide[7])-50
                      
                    }
                    
                    if (((as.numeric(input$RT) - (as.numeric(Guide[3])+y-x)) >= 4) != TRUE) 
                    {
                      
                      Guide[7] <- as.numeric(Guide[7])-6
                      
                    }
                    
                    Guide[7] <- as.numeric(Guide[7]) - as.numeric(Guide[3])
                    
                    ## Add Guide[8] with the html tagged edit in red and bold on the coding strand
                    
                    Guide[8] <- paste(substring(Guide[4],1,as.numeric(input$RT)-((as.numeric(Guide[3])+z))),as.character(DNAchange_char_tagged), substring(Guide[4], first = 1-z+y+u+((as.numeric(input$RT)-(as.numeric(Guide[3]))))), sep = "")
                    
                    guides <- rbind(guides,Guide)
                    
 
                  }
                  }
                }
                  
                }
              }
            
            
          }
          
            
            
          if (is.null(guides)!= TRUE)
          {
            
            
            if(any(((as.numeric(guides[,3])+u+z) <= input$RT) || (as.numeric(guides[,3])+u-z >= 0)))
            {
              guides <- subset(guides, as.numeric(guides[,3])+u+z <= input$RT)
              guides <- subset(guides,(as.numeric(guides[,3])+u-z) >= 0)
              
            }else{
              
              d <- NULL
              c <- NULL
              
            }
          
            
            if(is.null(c) || is.null(d) != TRUE)
            {
          
          currentguidesdf <- reactive({ 
            
            guidesdf <- data.frame("Variant" = input$Variant,"Score" = as.numeric(guides[,7]), "Protospacer(Sense)" = guides[,2],"Protospacer(Antisense)" = as.character(reverseComplement(DNAStringSet(guides[,2]))) , "EditPos." = guides[,3], "Extension(coding_strand)" = guides[,8], "PAM" = guides[,5], "PAM-Strand" = guides[,6],  stringsAsFactors = FALSE)
            guidesdf[order(guidesdf$Score,decreasing = TRUE), ]
            
          }) 
            }
          
      
          }
          }else{
            
            output$myText <- renderText({
              
              "Could not find a matching sequence, please check your genomic coordinates!"
            })
            
            guides <- 1
          }
          
          }
        
     
        
        if((as.character(input$`Gene Orientation`) == "+") == TRUE || input$SequenceInput == "Sequence input")
        {
          c <- 1
          d <- 1
          
          
          ### Get sequence for Plus Strand
          
          if (input$SequenceInput == "Sequence input")
          {
            
            
            RIGHT = function(x,n){
              substring(x,nchar(x)-n+1)
            }

            if ((substr(input$Edit2,1,3) == "ins")==TRUE){
              
              Seq_plus <- DNAStringSet(paste(RIGHT(input$UpstreamSequence,75),substring(input$DownstreamSequence,0,75), sep = ""))

              ### Get Sequence for Minus Strand
              
              Seq_minus <- reverseComplement(Seq_plus)
              
              Seq_plus_rev <- reverse(Seq_plus)
              Seq_minus_rev <- reverse(Seq_minus)
              
              Seq_minus_char <- as.character(Seq_minus)
              Seq_plus_char <- as.character(Seq_plus)
              Seq_plus_rev_char <- as.character(Seq_plus_rev)
              Seq_minus_rev_char <- as.character(Seq_minus_rev)
              
              
              DNAchange <- DNAStringSet(substring(input$Edit2, first = 4))
              anti_DNAchange <- complement(DNAchange)
              
              DNAchange_char <- as.character(DNAchange)
              
              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              p <- 1
              
              x <- nchar(as.character(DNAchange))
              u <- nchar(as.character(DNAchange))
              
              
              Seq_plus_char_n <- paste(substr(Seq_plus_char, start = 1, stop = 75),as.character(DNAchange), substr(Seq_plus_char, start = 76, stop = 151), sep = "")
              Seq_minus_rev_char_n <- paste(substr(Seq_minus_rev_char, start = 1, stop = 75),as.character(anti_DNAchange), substr(Seq_minus_rev_char, start = 76, stop = 151), sep = "")
              
            }
            else if ((substr(input$Edit2,1,3) == "del")==TRUE){
  
              
              
              z <- nchar(substring(input$Edit2, first = 4))
              
              Seq_plus <- DNAStringSet(paste(RIGHT(input$UpstreamSequence,75), substring(input$Edit2, first = 4), substring(input$DownstreamSequence,0,76-z), sep = ""))
              
              ### Get Sequence for Minus Strand
              
              Seq_minus <- reverseComplement(Seq_plus)
              
              Seq_plus_rev <- reverse(Seq_plus)
              Seq_minus_rev <- reverse(Seq_minus)
              
              Seq_minus_char <- as.character(Seq_minus)
              Seq_plus_char <- as.character(Seq_plus)
              Seq_plus_rev_char <- as.character(Seq_plus_rev)
              Seq_minus_rev_char <- as.character(Seq_minus_rev)
              
              
              DNAchange <- DNAStringSet(substring(input$Edit2, first = 4))
              anti_DNAchange <- complement(DNAchange)
              
              
              DNAchange_char <- paste("[",as.character(DNAchange),"]" ,sep = "")
              
              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- paste("[",as.character(reverse(anti_DNAchange)),"]" ,sep = "")
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              
              x <- nchar(as.character(DNAchange))
              
              x <- as.numeric(x)* (-1)  
              
              m <- z-1
              
              
              if(((substr(Seq_plus, start = 76, stop = 75-x)) == (as.character(DNAchange)))==TRUE){
                
                Seq_plus_char_n <- paste(substr(Seq_plus_char, start = 1, stop = 75), substr(Seq_plus_char, start = 76-x, stop = 151), sep = "")
                
                Seq_minus_rev_char_n <- paste(substr(Seq_minus_rev_char, start = 1, stop = 75), substr(Seq_minus_rev_char, start = 76-x, stop = 151), sep = "")
                
                
              }else{
                
                b <- NULL
                
              }
              
              
            }
            else if(grepl(">", input$Edit2) == TRUE) { 
              
              y <- ((nchar(as.character(input$Edit2)))-1)/2
              
              Seq_plus <- DNAStringSet(paste(RIGHT(input$UpstreamSequence,75), substr(input$Edit2,0,y), substring(input$DownstreamSequence,0,75), sep = ""))
              
              ### Get Sequence for Minus Strand
              
              Seq_minus <- reverseComplement(Seq_plus)
              
              Seq_plus_rev <- reverse(Seq_plus)
              Seq_minus_rev <- reverse(Seq_minus)
              
              Seq_minus_char <- as.character(Seq_minus)
              Seq_plus_char <- as.character(Seq_plus)
              Seq_plus_rev_char <- as.character(Seq_plus_rev)
              Seq_minus_rev_char <- as.character(Seq_minus_rev)
              
              v <- 1
              
              
              DNAchange <- DNAStringSet(substr(input$Edit2,2+y,2+(y*2)))
              anti_DNAchange <- complement(DNAchange)
              
              DNAchange_char <- as.character(DNAchange)
              
              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              if ((substr(input$Edit2,0,y) == substr(Seq_plus, start = 76, stop = 75+y)) == TRUE ){
                
                Seq_plus_char_n <- Seq_plus_char
                Seq_minus_rev_char_n <- Seq_minus_rev_char
                
                substr(Seq_plus_char_n, start = 76, stop = 75+y) <- as.character(DNAchange)
                
                substr(Seq_minus_rev_char_n, start = 76, stop = 75+y) <- as.character(anti_DNAchange)
                
              }else{
                
                a <- NULL
                
              }
            }
            
          }
          else if (input$SequenceInput == "Genomic coordinates")
          {
            
            
          if(is.na(seqlengths(checkCompatibleSeqinfo(Target_Genome, GRanges(as.character(Chromosome),ranges = IRanges(start = as.numeric(input$Position)+as.numeric(paste(input$`Gene Orientation`,"75", sep = "")), width = 151)))[c(as.character(Chromosome))])) != TRUE)
          {
              
              
          Seq_plus <- getSeq(Target_Genome, GRanges(as.character(Chromosome),ranges = IRanges(start = as.numeric(input$Position)-as.numeric(paste(input$`Gene Orientation`,"75", sep = "")), width = 151)),)
          
          ### Get Sequence for Minus Strand
          
          Seq_minus <- reverseComplement(Seq_plus)
          
          Seq_plus_rev <- reverse(Seq_plus)
          Seq_minus_rev <- reverse(Seq_minus)
          
          Seq_minus_char <- as.character(Seq_minus)
          Seq_plus_char <- as.character(Seq_plus)
          Seq_plus_rev_char <- as.character(Seq_plus_rev)
          Seq_minus_rev_char <- as.character(Seq_minus_rev)
          
          if ((substr(input$Edit,1,3) == "ins")==TRUE){
            
            DNAchange <- DNAStringSet(substring(input$Edit, first = 4))
            anti_DNAchange <- complement(DNAchange)
            
            DNAchange_char <- as.character(DNAchange)
            
            DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
            
            anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            p <- 1
            
            x <- nchar(as.character(DNAchange))
            u <- nchar(as.character(DNAchange))
            
            x <- as.numeric(x)
            
            Seq_plus_char_n <- paste(substr(Seq_plus_char, start = 1, stop = 75),as.character(DNAchange), substr(Seq_plus_char, start = 76, stop = 151), sep = "")
            Seq_minus_rev_char_n <- paste(substr(Seq_minus_rev_char, start = 1, stop = 75),as.character(anti_DNAchange), substr(Seq_minus_rev_char, start = 76, stop = 151), sep = "")
            
          }
          else if ((substr(input$Edit,1,3) == "del")==TRUE){
            
            DNAchange <- DNAStringSet(substring(input$Edit, first = 4))
            
            anti_DNAchange <- complement(DNAchange)
            
            
            DNAchange_char <- paste("[",as.character(DNAchange),"]" ,sep = "")
            
            DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            anti_DNAchange_char <- paste("[",as.character(reverse(anti_DNAchange)),"]" ,sep = "")
            
            anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            z <- nchar(as.character(DNAchange))
            
            x <- nchar(as.character(DNAchange))
            
            x <- as.numeric(x)* (-1)  
            
            m <- z-1
            
            
            if(((substr(Seq_plus, start = 76, stop = 75-x)) == (as.character(DNAchange)))==TRUE){
              
              Seq_plus_char_n <- paste(substr(Seq_plus_char, start = 1, stop = 75), substr(Seq_plus_char, start = 76-x, stop = 151), sep = "")
              
              Seq_minus_rev_char_n <- paste(substr(Seq_minus_rev_char, start = 1, stop = 75), substr(Seq_minus_rev_char, start = 76-x, stop = 151), sep = "")
              
              
            }else{
              
              b <- NULL
            
            }
            
            
          }
          else if(grepl(">", input$Edit) == TRUE) { 
            
            
       
            y <- ((nchar(input$Edit))-1)/2
            
            v <- 1
            
            
            DNAchange <- DNAStringSet(substr(input$Edit,2+y,2+(y*2)))
            anti_DNAchange <- complement(DNAchange)
            
            DNAchange_char <- as.character(DNAchange)
            
            DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
            
            anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
            
            if ((substr(input$Edit,0,y) == substr(Seq_plus, start = 76, stop = 75+y)) == TRUE ){
              
              Seq_plus_char_n <- Seq_plus_char
              Seq_minus_rev_char_n <- Seq_minus_rev_char
              
              substr(Seq_plus_char_n, start = 76, stop = 75+y) <- as.character(DNAchange)
              
              substr(Seq_minus_rev_char_n, start = 76, stop = 75+y) <- as.character(anti_DNAchange)
              
            }else{
              
              a <- NULL

            }
          }
          
          }else{
  
            s <- NULL
            
          }
            
          }
          
          
          ## Searchig for sense PAM in Prime editing 
          
          
          if(is.null(s) != TRUE)
          {
          
          ## Upstream PAM for PRIME 
          
          primeUpstreamPAM <- substr(Seq_plus_char, 5,81)  
          
          primeUpstreamPAM2 <- substr(Seq_minus_rev_char,71,147)
          
          
          
          
          ## prime edits possible:
          
          
          ## Sense strand
          
          primeeditable <- NULL;
          
          
          for(i in 1:length(primeUpstreamPAM))
          {
            if (grepl( "AGG", primeUpstreamPAM[i]) == TRUE || grepl( "CGG", primeUpstreamPAM[i]) == TRUE || grepl( "GGG", primeUpstreamPAM[i]) == TRUE || grepl( "TGG", primeUpstreamPAM[i]) == TRUE ){
              
              
              temp <- c(i,primeUpstreamPAM[i])
              primeeditable <- rbind(primeeditable, temp)
              
              
            }
          }
          
          
          primeeditabledf <- data.frame(primeeditable)
          
          primeeditabledf$X2 <- as.character(primeeditabledf$X2)
          
          primeeditabledf$X1 <- as.numeric(primeeditabledf$X1)
          
          
          
          
          ## AntiSense strand
          
          primeeditable2 <- NULL;
          
          
          for(i in 1:length(primeUpstreamPAM2))
          {
            if (grepl( "GGA", primeUpstreamPAM2[i]) == TRUE || grepl( "GGC", primeUpstreamPAM2[i]) == TRUE || grepl( "GGG", primeUpstreamPAM2[i]) == TRUE || grepl( "GGT", primeUpstreamPAM2[i]) == TRUE ){
              
              
              temp <- c(i,primeUpstreamPAM2[i])
              primeeditable2 <- rbind(primeeditable2, temp)
              
              
            }
          }
          
          
          primeeditable2df <- data.frame(primeeditable2)
          
          primeeditable2df$X2 <- as.character(primeeditable2df$X2)
          
          primeeditable2df$X1 <- as.numeric(primeeditable2df$X1)
          
          
          if(is.null(a) != TRUE & is.null(b) != TRUE){
            
            if(is.null(primeeditable) != TRUE){
            
            
            ## Sense strand
            
            for (i in 1:nrow(primeeditabledf))
              for (j in 1:(nchar(primeeditabledf[1,2])-2))
              {
                {
                  if (grepl( "AGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(primeeditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE ){

                      
                    if((nchar(substr(Seq_plus_char[primeeditabledf[i,1]],(j-16),j+3)) == 20))
                    {
                     
                      Guide <- c("1",substr(Seq_plus_char[primeeditabledf[i,1]],j-16,j+3),76-j,as.character(reverse(DNAStringSet(substr(Seq_minus_rev_char_n[primeeditabledf[i,1]],(j+1) - as.numeric(input$PBS),((j+1) - (as.numeric(input$PBS)))+ (as.numeric(input$RT)+as.numeric(input$PBS)-1))))),substr(primeeditabledf[i,2], j,j+2), "Sense")
                      
                      if(substr(Seq_plus_char[primeeditabledf[i,1]],j-16,j-16) != "G")
                      {
                        
                        Guide[2] <- paste("g",Guide[2], sep = "") 
                        
                      }
                      
                      Guide[7] <- 0
                      
                      if((substr(Guide[4],1,1) == "C") == TRUE){
                        
                        Guide[7] <- as.numeric(Guide[7])-28
                        
                      }
                      
                      if ((grepl("TTTTT", Guide[4])) == TRUE){
                        
                        Guide[7] <- as.numeric(Guide[7])-50
                        
                      }
                      
                      if (((as.numeric(input$RT) - (as.numeric(Guide[3])+y+x)) >= 4) != TRUE)
                      {
                        
                        Guide[7] <- as.numeric(Guide[7])-6
                        
                      }
                      
                      
                      Guide[7] <- as.numeric(Guide[7])-as.numeric(Guide[3])
                      Guide[8] <- paste(reverseComplement(DNAStringSet(substring(Guide[4], first = 2+((as.numeric(input$RT)-(as.numeric(Guide[3]))))))),as.character(DNAchange_char_tagged),reverseComplement(DNAStringSet(substring(Guide[4],1,1+as.numeric(input$RT)-((as.numeric(Guide[3])+y+x+z))))), sep = "")

                      guides <- rbind(guides,Guide)
                      
                    }
                  
                    
                  } 
                }
              }
            }
          }
          
          if(is.null(a) != TRUE & is.null(b) != TRUE){
            
            if(is.null(primeeditable2) != TRUE){
            
              ## AntiSense strand
            
            for (i in 1:nrow(primeeditable2df))
              for (j in 1:(nchar(primeeditable2df[1,2])-2))
              {
                {
                  if (grepl( "GGA", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(primeeditable2df[i,2],j,j+2)) == TRUE || grepl( "GGT", substr(primeeditable2df[i,2], j,j+2)) == TRUE ){
                    
                    if((j+u-z) >= 0)
                       {
                      
                    if(nchar(substr(Seq_minus_rev_char[primeeditable2df[i,1]],j+73,j+92)) == 20){
                      
                      Guide <- c("1",substr(Seq_minus_rev_char[primeeditable2df[i,1]],j+73,j+92),j, as.character(reverse(DNAStringSet(substr(Seq_plus_char_n[primeeditable2df[i,1]],((j+76+x)+(as.numeric(input$PBS)))-(as.numeric(input$RT)+as.numeric(input$PBS)),((j+76+x)+(as.numeric(input$PBS)))-1)))),as.character(reverse(DNAStringSet(substr(primeeditable2df[i,2], j,j+2)))), "Antisense")
                      
                      Guide[2] <- as.character(reverse(DNAStringSet(Guide[2])))
                      
                      Guide[4] <- as.character(reverse(DNAStringSet(Guide[4])))
                      
                      if((substr(Seq_minus_rev_char[primeeditable2df[i,1]],j+92,j+92) == "G") == FALSE)
                      {
                        
                        Guide[2] <- paste("g",Guide[2], sep = "") 
                        
                      }
                      
                     
                      
                      Guide[7] <- 0
                      
                      if((substr(Guide[4],1,1) == "C") == TRUE){
                        
                        Guide[7] <- as.numeric(Guide[7])-28
                        
                      }
                      
                      if ((grepl("TTTTT", Guide[4])) == TRUE){
                        
                        Guide[7] <- as.numeric(Guide[7])-50
                        
                      }
                      
                      if (((as.numeric(input$RT) - (as.numeric(Guide[3])+y-x)) >= 4) != TRUE)   
                      {
                        
                        Guide[7] <- as.numeric(Guide[7])-6
                        
                      }
                      
                      
                      
                      Guide[7] <- as.numeric(Guide[7])-as.numeric(Guide[3])
                      
                      Guide[8] <-paste(substring(Guide[4],1,1+as.numeric(input$RT)-((as.numeric(Guide[3])+v+x+z+p-m))),as.character(DNAchange_char_tagged), substring(Guide[4], first = 2+((as.numeric(input$RT)-(as.numeric(Guide[3])+p-m-y+v)))), sep = "")
                      
                      
                      
                      
                      guides <- rbind(guides,Guide)
                      
                    }
                    }
                  }

                }
              }

            }
          }
          
          if (is.null(guides)!= TRUE)
          {
            if(any((as.numeric(guides[,3])+u+z) <= input$RT) == TRUE)
            {
              guides <- subset(guides, as.numeric(guides[,3])+u+z <= input$RT)

            }else{
              
              d <- NULL
              c <- NULL
              
            }

            if(is.null(c) || is.null(d) != TRUE)
            {
          
          
          currentguidesdf <- reactive({ 
            
            guidesdf <- data.frame("Variant" = input$Variant, "Score" = as.numeric(guides[,7]),"Protospacer(Sense)" = guides[,2],"Protospacer(Antisense)" = as.character(reverseComplement(DNAStringSet(guides[,2]))) , "EditPos." = guides[,3], "Extension(coding.strand)" = guides[,8], "PAM" = guides[,5], "PAM-Strand" = guides[,6], stringsAsFactors = FALSE)
            guidesdf[order(guidesdf$Score,decreasing = TRUE), ]
            
          })
          
            }
            
          
          }
          
          
          
          
          }
          else if(is.null(s) == TRUE)
          {
            print("wrong coordinates")
          }
        }
          


        
        if(is.null(a) == TRUE){
          
          output$myText <- renderText({
            
            "Could not find a matching sequence, please check your Substitution!"
          })
        }
        else if(is.null(b) == TRUE){
          
          output$myText <- renderText({
            
            "Could not find a matching sequence, please check your Deletion!"
          })
          
        }
        else if(is.null(c) == TRUE && is.null(d) == TRUE ){
          
          output$myText <- renderText({
            
            "Edit is too far away, you might increase the RT length"
          })
        }
        else if(is.null(guides)== TRUE & is.null(s) != TRUE){
          
          output$myText <- renderText({
            
            "We couldn't find a Prime editing Guide!"
          })
          
          
        }
        else if(is.null(s) == TRUE){

          output$myText <- renderText({
            
            "Could not find a matching sequence, please check your genomic coordinates!"
          })

        
        }else{

          ### Make Oligo table, for ready to clone oligos 
          
          
          currentguidesdf2 <- reactive({ 
            
            
            guidesdf <- data.frame("Variant" = input$Variant, "Protospacer(Sense)" = paste("cacc", guides[,2],"gtttt", sep = ""),"Protospacer(Antisense)" = paste("ctctaaaac",reverseComplement(DNAStringSet(guides[,2])),sep = "") , "EditPos." = guides[,3], "Extension(Sense)" = paste("gtgc",guides[,4],sep = ""), "Extension(Antisense)" = paste("aaaa",reverseComplement(DNAStringSet(guides[,4])), sep = ""), "PAM" = guides[,5], "PAM-Strand" = guides[,6], "Score" = as.numeric(guides[,7]), stringsAsFactors = FALSE)
            guidesdf[order(guidesdf$Score,decreasing = TRUE), ]
            
          }) 
          
        
          output$mytable  <- DT::renderDataTable({ currentguidesdf()},
                                                 escape = FALSE,
                                                 caption = "Output Table" ,
                                                 options = list(dom = 't', initComplete =  JS("function(settings, json) {",
                                                 "$('body').css({'font-family': 'Helvetica'});",
                                                 "}")),
                                                  selection = list(mode = 'single', target = 'column'),
                                                  callback=JS(
                                                    'table.on("click.dt","td", function() {
                                                    var data=table.cell(this).data();
                                                    if (table.cell(this).data() < 0)
                                                    swal({title: "pegRNA-Score", text: "The higher the better! (calculated based on recommendation from the Liu Lab)",imageUrl: "Score.jpg",imageSize: "460x200"
                                                    });   
                                                    if (table.cell(this).data() > 0)
                                                    swal({title: "Edit position", text:  "Edit position is defined relative to the PAM sequence (Anzalone et.al. 2019)",imageUrl: "Target.jpg",imageSize: "460x200"
                                                    });                                                    
 
                                                    })')
                                               
                                                  
                                                 ) 
          
          
      
          
          
          
        
          output$mytable2 <- DT::renderDataTable({ currentguidesdf2()}, caption = "Ordering Table",options = list(dom = 't'))
          
          
        }
        
        
      }
        
      
        
        output$downloadData <- downloadHandler(
          filename = function() {
            paste('pegRNA Oligos for cloning-', Sys.Date(), '.csv', sep='')
          },
          content = function(con) {
            write.csv(data.frame("Variant" = input$Variant, "Protospacer(Sense)" = paste("cacc", guides[,2],"gtttt", sep = ""),"Protospacer(Antisense)" = paste("ctctaaaac",reverseComplement(DNAStringSet(guides[,2])),sep = "") , "TargetPos." = guides[,3], "Extension(Sense)" = paste("gtgc",guides[,4],sep = ""), "Extension(Antisense)" = paste("aaaa",reverseComplement(DNAStringSet(guides[,4])), sep = ""), "PAM" = guides[,5], "PAM-Strand" = guides[,6], "Score" = as.numeric(guides[,7]), "PBS" = input$PBS, "RTT" = input$RT, stringsAsFactors = FALSE), con, sep = ",", row.names=FALSE)
          }
        )
        

      } 
      }
      
#-----------------------------------------------------------------------------------------------------------------#    
#-----------------------------------------------------------------------------------------------------------------#  
    
    ## Start of a Multi Sample Run, if this option is selected in the User Interface:
      
      if(input$Mode == 'Multi Sample Run') {
        

        withProgress(message = 'Generating Oligos', value = 0, {
        
        inFile <- input$file1
        
        inFile <- read.csv(inFile$datapath, header = TRUE, sep = ";", stringsAsFactors = FALSE)
        
        if(input$Editing == "Base editing")
        {
          
          guides <- NULL
          
          Download <- NULL
          
          output$downloadData <- downloadHandler(
            filename = function() {
              paste('Base editing guides-', Sys.Date(), '.csv', sep='')
            },
            content = function(con) {
              write.csv(data.frame("Variant" = guides[,6],"Protospacer" = Download, "EditPos." = guides[,3], "PAM" = guides[,4], "Base Editor" = guides[,5]), con, sep = ",", row.names=FALSE)
            }
          )

          
          
          for (k in 1:nrow(inFile))
          {
            if(input$Genomes != "Rice (MSU7)")
            {
              
              Chromosome <- as.character(paste("chr",inFile[k,]$Chromosome, sep = ""))
                                            
              
            }
            else if(input$Genomes == "Rice (MSU7)")
            {
              
              Chromosome <- as.character(paste("Chr",inFile[k,]$Chromosome, sep = ""))
              
            }
            
          
          ### Get sequence for Plus Strand
          
          Seq_plus <- getSeq(Target_Genome, GRanges(Chromosome,ranges = IRanges(start = c(as.numeric(inFile[k,]$GenomicLocation)-25), width = 51)),)
          
          ### Get Sequence for Minus Strand
          
          Seq_minus <- reverseComplement(Seq_plus)
          
          Guide <- NULL
          
          
          ## Insert the Patient mutation in the Plus Strand at string position 26
          ## Code for G>A on Plus Strand and C>T on Minus Strand
          
          Seq_plus_char <- as.character(Seq_plus)
          Seq_plus_char_cntrl <- as.character(Seq_plus)
          
          
          Seq_minus_char <- as.character(Seq_minus)
          Seq_minus_char_cntrl <- as.character(Seq_minus)
          
          
          # Create additional character strings, in which the mutation found in the patient will be included:
          # '-NGG' is for ABEmax with a NGG PAM; '-NGA' for NG-ABEmax/VRQR-ABEmax with a NGA PAM; '-NGCG' for NG-ABEmax;
          # '-NNGRRT' is for Sa-ABEmax; '-NNHRRT' is for SaKKH-ABEmax
          
          Seq_minus_char_NGG <- Seq_minus_char
          
          Seq_minus_char_NGA <- Seq_minus_char
          
          Seq_minus_char_NGCG <- Seq_minus_char
          
          Seq_minus_char_NNGRRT <- Seq_minus_char
          
          Seq_minus_char_NNHRRT <- Seq_minus_char
          
          Seq_minus_char_NGC <- Seq_minus_char
          
          Seq_minus_char_NYN <- Seq_minus_char
          
          
          Seq_plus_char_NGG <- Seq_plus_char
          
          Seq_plus_char_NGA <- Seq_plus_char
          
          Seq_plus_char_NGCG <- Seq_plus_char
          
          Seq_plus_char_NNGRRT <- Seq_plus_char
          
          Seq_plus_char_NNHRRT <- Seq_plus_char
          
          Seq_plus_char_NGC <- Seq_plus_char
          
          Seq_plus_char_NYN <- Seq_plus_char
          
          
          
          
          incProgress(1/nrow(inFile), detail = paste("Variant", k))
          
          
          # Guide search For ABEmax (NGG PAM) depending on input'-SNP' and '-GeneOrientation':
          # For Base editing guides on the Minus Strand:
          
          if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "-")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "+"))){
            
            # Test wildtype sequence if inFile[k,]'-SNP' was assigned right and if so includes patient mutation into character string:
            # If not, see code line 216
            
            if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) == TRUE)
            {
              substr(Seq_minus_char_NGG, start = 26, stop = 26) <- "A"
              
              # patient character string is now searched for suitable PAM sequences:
              
              # Perfect edible window: For ABEmax (Editing Position: 6-8)
              # Substring for PAM is selected based on the Editing window:
              
              PerfectUpstreamPAM <- as.character(substr(Seq_minus_char_NGG, 39,43))  
              
              # 'minus_perfecteditable' object is created to store found PAM candidates:
              
              minus_perfecteditable <- NULL
              
              # The PerfectUpstreamPAM substring is searched for all different PAM sequences
              # Matches are stored in the 'minus_perfecteditable' object
              
              for(i in 1:length(PerfectUpstreamPAM))
              {
                if (grepl( "AGG", PerfectUpstreamPAM[i]) == TRUE || grepl( "GGG", PerfectUpstreamPAM[i]) == TRUE || grepl( "CGG", PerfectUpstreamPAM[i]) == TRUE || grepl( "TGG", PerfectUpstreamPAM[i]) == TRUE ){    
                  #print(c(i,PerfectUpstreamPAM[i]))
                  tmp1 <- c(i,PerfectUpstreamPAM[i])
                  minus_perfecteditable <- rbind(minus_perfecteditable,tmp1)
                  
                  
                }
              }
              
              # Quality check to only run code, if suitable PAM sequences are found:
              
              if(is.null(minus_perfecteditable) != TRUE) {
                
                minus_perfecteditabledf <- data.frame(minus_perfecteditable)
                
                minus_perfecteditabledf$X2 <- as.character(minus_perfecteditabledf$X2)
                
                minus_perfecteditabledf$X1 <- as.numeric(nrow(minus_perfecteditable))
                
                
                # Guide sequences are designed based on the position of the PAM sequence in the 'minus_perfecteditable' object:
                
                for (i in 1:nrow(minus_perfecteditabledf))
                  for (j in 1:(nchar(minus_perfecteditabledf[1,2])-2))
                  {
                    {
                      if (grepl( "AGG", substr(minus_perfecteditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(minus_perfecteditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(minus_perfecteditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(minus_perfecteditabledf[i,2], j,j+2)) == TRUE ){
                        
                        # Guide vector is created with the guide number, Protospacer sequence, Edit Position, PAM sequence, Base editor
                        # Resulting vectors are combined to create 'guides'
                        
                        Guide <- c(row.names(inFile[k,]),substr(Seq_minus_char_NGG[minus_perfecteditabledf[i,1]],18+j,37+j),26-(17+j), substr(minus_perfecteditabledf[i,2], j,j+2),"ABEmax/ABE8e",inFile[k,]$Variant)
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(minus_perfecteditabledf[1,2])-2))
                          
                        {
                          if(grepl("A", substr(Guide[2],5+l,5+l)) == TRUE)
                          {
                            
                            tmp <- 5+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        
                        nrow(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                        
                      }
                    }
                  }
                
              }
              
              
            }
            
            # If check of wildtype sequence with inFile[k,]'-SNP' does not match, a gets assigned value of NULL
            
            else if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) != TRUE){
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
            }
            
            
          }
          
          # For Base editing guides on the Plus Strand:
          
          else if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "+")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "-"))){
            
            if (("G" == substr(Seq_plus_char_cntrl, start = 26, stop = 26)) == TRUE )
            {
              substr(Seq_plus_char_NGG, start = 26, stop = 26) <- "A"
              
              
              
              # For Base editing guides on the Plus Strand:
              
              # Perfect edible window: For ABEmax (Editing Position: 6-8)
              # Substring for PAM is selected based on the Editing window:
              
              PerfectDownstreamPAM <- as.character(substr(Seq_plus_char_NGG, 39,43))  
              
              
              ## Test if there is a suitable Downstream PAM sequence
              
              plus_perfecteditable <- NULL
              
              for(i in 1:length(PerfectDownstreamPAM))
              {
                if (grepl( "AGG", PerfectDownstreamPAM[i]) == TRUE || grepl( "GGG", PerfectDownstreamPAM[i]) == TRUE || grepl( "CGG", PerfectDownstreamPAM[i]) == TRUE || grepl( "TGG", PerfectDownstreamPAM[i]) == TRUE ){    
                  #print(c(i,PerfectDownstreamPAM[i]))
                  tmp <- c(i,PerfectDownstreamPAM[i])
                  plus_perfecteditable <- rbind(plus_perfecteditable,tmp)
                  
                  
                }
              }
              
              if(is.null(plus_perfecteditable) != TRUE) {
                
                
                plus_perfecteditabledf <- data.frame(plus_perfecteditable)
                
                plus_perfecteditabledf$X2 <- as.character(plus_perfecteditabledf$X2)
                
                plus_perfecteditabledf$X1 <- as.numeric(nrow(plus_perfecteditable))
                
                
                
                
                for (i in 1:nrow(plus_perfecteditabledf))
                  for (j in 1:(nchar(plus_perfecteditabledf[1,2])-2))
                  {
                    {
                      if (grepl( "AGG", substr(plus_perfecteditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(plus_perfecteditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(plus_perfecteditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(plus_perfecteditabledf[i,2], j,j+2)) == TRUE ){
                        
                        
                        
                        #print(plus_perfecteditabledf[i,2])
                        #print(FANCA_char[perfect_editable_sense[i,1]])
                        Guide <- c(row.names(inFile[k,]),substr(Seq_plus_char_NGG[plus_perfecteditabledf[i,1]],18+j,37+j),26-(17+j),substr(plus_perfecteditabledf[i,2], j,j+2),"ABEmax/ABE8e",inFile[k,]$Variant)
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(plus_perfecteditabledf[1,2])-2))
                          
                        {
                          if(grepl("A", substr(Guide[2],5+l,5+l)) == TRUE)
                          {
                            
                            tmp <- 5+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        
                        nrow(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                      }
                    }
                  }
                
              }
            }else{

              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
              
            }
            
          }
          
          #---------------------------------------------#        
          
          # Guide search For NG-ABEmax, VRQR-ABEmax (NGA PAM) depending on inFile[k,]'-SNP' and '-GeneOrientation':
          # For Base editing guides on the Minus Strand:
          
          if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "-")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "+"))){
            
            
            if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) == TRUE )
            {
              substr(Seq_minus_char_NGA, start = 26, stop = 26) <- "A"
              
              
              
              ## Perfect edible: NG-ABEmax window 4-6 (PAM: NGA)
              
              PerfectUpstreamPAM_NGA <- substr(Seq_minus_char_NGA, 41,45)  
              
              
              ## Target Plus Strand C>T or Minus Strand G>A
              
              minus_perfecteditable_NGA <- NULL;
              
              for(i in 1:length(PerfectUpstreamPAM_NGA))
              {
                if (grepl( "AGA", PerfectUpstreamPAM_NGA[i]) == TRUE || grepl( "GGA", PerfectUpstreamPAM_NGA[i]) == TRUE || grepl( "CGA", PerfectUpstreamPAM_NGA[i]) == TRUE || grepl( "TGA", PerfectUpstreamPAM_NGA[i]) == TRUE ){    
                  #print(c(i,PerfectUpstreamPAM_NGA[i]))
                  tmp <- c(i,PerfectUpstreamPAM_NGA[i])
                  minus_perfecteditable_NGA <- rbind(minus_perfecteditable_NGA,tmp)
                  
                  
                }
              }
              
              if(is.null(minus_perfecteditable_NGA) != TRUE) {
                
                minus_perfecteditabledf_NGA <- data.frame(minus_perfecteditable_NGA)
                
                minus_perfecteditabledf_NGA$X2 <- as.character(minus_perfecteditabledf_NGA$X2)
                
                minus_perfecteditabledf_NGA$X1 <- as.numeric(nrow(minus_perfecteditable_NGA))
                
                
                for (i in 1:nrow(minus_perfecteditabledf_NGA))
                  for (j in 1:(nchar(minus_perfecteditabledf_NGA[1,2])-2))
                  {
                    {
                      if (grepl( "AGA", substr(minus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE || grepl( "GGA", substr(minus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE || grepl( "CGA", substr(minus_perfecteditabledf_NGA[i,2],j,j+2)) == TRUE || grepl( "TGA", substr(minus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE ){
                        
                        #print(minus_perfecteditabledf_NGA[i,2])
                        #print(FANCA_char[perfect_editable_sense[i,1]])
                        Guide <- c(row.names(inFile[k,]),substr(Seq_minus_char_NGA[minus_perfecteditabledf_NGA[i,1]],20+j,39+j),26-(19+j),substr(minus_perfecteditabledf_NGA[i,2], j,j+2),"NG-ABEmax/ABE8e, VRQR-ABEmax",inFile[k,]$Variant)
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(minus_perfecteditabledf_NGA[1,2])-2))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 3
                          
                        {
                          if(grepl("A", substr(Guide[2],3+l,3+l)) == TRUE)
                          {
                            
                            tmp <- 3+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                      }
                    }
                  }
                
              }
            }else{
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
            }
            
          }
          
          # For Base editing guides on the Plus Strand:
          
          else if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "+")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "-"))){
            
            
            if (("G" == substr(Seq_plus_char_cntrl, start = 26, stop = 26)) == TRUE )
            {
              substr(Seq_plus_char_NGA, start = 26, stop = 26) <- "A"
              
              
              
              ## Perfect edible: NG-ABEmax window 4-6 (PAM: NGA)
              
              PerfectDownstreamPAM_NGA <- substr(Seq_plus_char_NGA, 41,45) 
              
              ## Target Plus Strand G>A or Minus Strand C>T (with NG-ABEmax)
              
              
              plus_perfecteditable_NGA <- NULL;
              
              for(i in 1:length(PerfectDownstreamPAM_NGA))
              {
                if (grepl( "AGA", PerfectDownstreamPAM_NGA[i]) == TRUE || grepl( "GGA", PerfectDownstreamPAM_NGA[i]) == TRUE || grepl( "CGA", PerfectDownstreamPAM_NGA[i]) == TRUE || grepl( "TGA", PerfectDownstreamPAM_NGA[i]) == TRUE ){    
                  #print(c(i,PerfectDownstreamPAM_NGA[i]))
                  tmp <- c(i,PerfectDownstreamPAM_NGA[i])
                  plus_perfecteditable_NGA <- rbind(plus_perfecteditable_NGA,tmp)
                  
                  
                }
              }
              
              
              if(is.null(plus_perfecteditable_NGA) != TRUE) {
                
                plus_perfecteditabledf_NGA <- data.frame(plus_perfecteditable_NGA)
                
                plus_perfecteditabledf_NGA$X2 <- as.character(plus_perfecteditabledf_NGA$X2)
                
                plus_perfecteditabledf_NGA$X1 <- as.numeric(nrow(plus_perfecteditable_NGA))
                
                
                for (i in 1:nrow(plus_perfecteditabledf_NGA))
                  for (j in 1:(nchar(plus_perfecteditabledf_NGA[1,2])-2))
                  {
                    {
                      if (grepl( "AGA", substr(plus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE || grepl( "GGA", substr(plus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE || grepl( "CGA", substr(plus_perfecteditabledf_NGA[i,2],j,j+2)) == TRUE || grepl( "TGA", substr(plus_perfecteditabledf_NGA[i,2], j,j+2)) == TRUE ){
                        
                        #print(plus_perfecteditabledf_NGA[i,2])
                        #print(FANCA_char[perfect_editable_sense[i,1]])
                        Guide <- c(row.names(inFile[k,]),substr(Seq_plus_char_NGA[plus_perfecteditabledf_NGA[i,1]],20+j,39+j),26-(19+j),substr(plus_perfecteditabledf_NGA[i,2], j,j+2),"NG-ABEmax/ABE8e, VRQR-ABEmax",inFile[k,]$Variant)
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(plus_perfecteditabledf_NGA[1,2])-2))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 3
                          
                        {
                          if(grepl("A", substr(Guide[2],3+l,3+l)) == TRUE)
                          {
                            
                            tmp <- 3+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        
                        nrow(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                        
                      }
                    }
                  }
                
              }
            }else{
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
            }
            
            
          }

          #---------------------------------------------#
          
          # Guide search For NG-ABEmax (NGCG PAM) depending on inFile[k,]'-SNP' and '-GeneOrientation':
          # For Base editing guides on the Minus Strand:
          
          if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "-")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "+"))){
            
            # Test wildtype sequence if inFile[k,]'-SNP' was assigned right and if so includes patient mutation into character string:
            # If not, see code line 216
            
            if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) == TRUE)
            {
              substr(Seq_minus_char_NGCG, start = 26, stop = 26) <- "A"
              
              # patient character string is now searched for suitable PAM sequences:
              
              # Perfect edible window: For NG-ABEmax NGCG PAM (Editing Position: 5-7)
              # Substring for PAM is selected based on the Editing window:
              
              PerfectUpstreamPAM_NGCG <- as.character(substr(Seq_minus_char_NGCG, 40,45))  
              
              # 'minus_perfecteditable' object is created to store found PAM candidates:
              
              minus_perfecteditable_NGCG <- NULL
              
              # The PerfectUpstreamPAM_NGCG substring is searched for all different PAM sequences
              # Matches are stored in the 'minus_perfecteditable_NGCG' object
              
              for(i in 1:length(PerfectUpstreamPAM_NGCG))
              {
                if (grepl( "AGCG", PerfectUpstreamPAM_NGCG[i]) == TRUE || grepl( "GGCG", PerfectUpstreamPAM_NGCG[i]) == TRUE || grepl( "CGCG", PerfectUpstreamPAM_NGCG[i]) == TRUE || grepl( "TGCG", PerfectUpstreamPAM_NGCG[i]) == TRUE ){    
                  tmp1 <- c(i,PerfectUpstreamPAM_NGCG[i])
                  minus_perfecteditable_NGCG <- rbind(minus_perfecteditable_NGCG,tmp1)
                  
                  
                }
              }
              
              # Quality check to only run code, if suitable PAM sequences are found:
              
              if(is.null(minus_perfecteditable_NGCG) != TRUE) {
                
                minus_perfecteditable_NGCGdf <- data.frame(minus_perfecteditable_NGCG)
                
                minus_perfecteditable_NGCGdf$X2 <- as.character(minus_perfecteditable_NGCGdf$X2)
                
                minus_perfecteditable_NGCGdf$X1 <- as.numeric(nrow(minus_perfecteditable_NGCG))
                
                
                # Guide sequences are designed based on the position of the PAM sequence in the 'minus_perfecteditable_NGCG' object:
                
                for (i in 1:nrow(minus_perfecteditable_NGCGdf))
                  for (j in 1:(nchar(minus_perfecteditable_NGCGdf[1,2])-3))
                  {
                    {
                      if (grepl( "AGCG", substr(minus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE || grepl( "GGCG", substr(minus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE || grepl( "CGCG", substr(minus_perfecteditable_NGCGdf[i,2],j,j+3)) == TRUE || grepl( "TGCG", substr(minus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE ){
                        
                        # Guide vector is created with the guide number, Protospacer sequence, Edit Position, PAM sequence, Base editor
                        # Resulting vectors are combined to create 'guides'
                        
                        Guide <- c(row.names(inFile[k,]),substr(Seq_minus_char_NGCG[minus_perfecteditable_NGCGdf[i,1]],19+j,38+j),25-(17+j),substr(minus_perfecteditable_NGCGdf[i,2], j,j+3),"NG-ABEmax/ABE8e",inFile[k,]$Variant)
                        
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(minus_perfecteditable_NGCGdf[1,2])-3))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                          
                        {
                          if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                          {
                            
                            tmp <- 4+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                        
                      }
                    }
                  }
                
              }
              
              
            }
            
            # If check of wildtype sequence with inFile[k,]'-SNP' does not match, a gets assigned value of NULL
            
            else if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) != TRUE){
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
              
            }
            
            
          }
          
          # For Base editing guides on the Plus Strand:
          
          else if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "+")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "-"))){
            
            if (("G" == substr(Seq_plus_char_cntrl, start = 26, stop = 26)) == TRUE)
            {
              substr(Seq_plus_char_NGCG, start = 26, stop = 26) <- "A"
              
              # For Base editing guides on the Plus Strand:
              
              # Perfect edible window: For NG-ABEmax NGCG PAM (Editing Position: 5-7)
              # Substring for PAM is selected based on the Editing window:
              
              PerfectDownstreamPAM_NGCG <- as.character(substr(Seq_plus_char_NGCG, 40,45)) 
              
              ## Test if there is a suitable Downstream PAM sequence
              
              plus_perfecteditable_NGCG <- NULL
              
              for(i in 1:length(PerfectDownstreamPAM_NGCG))
              {
                if (grepl( "AGCG", PerfectDownstreamPAM_NGCG[i]) == TRUE || grepl( "GGCG", PerfectDownstreamPAM_NGCG[i]) == TRUE || grepl( "CGCG", PerfectDownstreamPAM_NGCG[i]) == TRUE || grepl( "TGCG", PerfectDownstreamPAM_NGCG[i]) == TRUE ){    
                  #print(c(i,PerfectDownstreamPAM_NGCG[i]))
                  tmp <- c(i,PerfectDownstreamPAM_NGCG[i])
                  plus_perfecteditable_NGCG <- rbind(plus_perfecteditable_NGCG,tmp)
                  
                  
                }
              }
              
              if(is.null(plus_perfecteditable_NGCG) != TRUE) {
                
                
                plus_perfecteditable_NGCGdf <- data.frame(plus_perfecteditable_NGCG)
                
                plus_perfecteditable_NGCGdf$X2 <- as.character(plus_perfecteditable_NGCGdf$X2)
                
                plus_perfecteditable_NGCGdf$X1 <- as.numeric(nrow(plus_perfecteditable_NGCG))
                
                
                
                
                for (i in 1:nrow(plus_perfecteditable_NGCGdf))
                  for (j in 1:(nchar(plus_perfecteditable_NGCGdf[1,2])-3))
                  {
                    {
                      if (grepl( "AGCG", substr(plus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE || grepl( "GGCG", substr(plus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE || grepl( "CGCG", substr(plus_perfecteditable_NGCGdf[i,2],j,j+3)) == TRUE || grepl( "TGCG", substr(plus_perfecteditable_NGCGdf[i,2], j,j+3)) == TRUE ){
                        
                        
                        
                        #print(plus_perfecteditable_NGCGdf[i,2])
                        #print(FANCA_char[perfect_editable_sense[i,1]])
                        Guide <- c(row.names(inFile[k,]),substr(Seq_plus_char_NGCG[plus_perfecteditable_NGCGdf[i,1]],19+j,38+j),25-(17+j),substr(plus_perfecteditable_NGCGdf[i,2], j,j+3),"NG-ABEmax/ABE8e",inFile[k,]$Variant)
                        
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(plus_perfecteditable_NGCGdf[1,2])-3))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                          
                        {
                          if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                          {
                            
                            tmp <- 4+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                        
                        
                      }
                    }
                  }
                
              }
            }
            else{
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
              
            }
            
          } 
          
          #---------------------------------------------#     
          
          # Guide search For Sa-ABEmax (NNGRRT PAM) depending on inFile[k,]'-SNP' and '-GeneOrientation':
          # For Base editing guides on the Minus Strand:
          
          if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "-")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "+"))){
            
            # Test wildtype sequence if inFile[k,]'-SNP' was assigned right and if so includes patient mutation into character string:
            # If not, see code line 216
            
            if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) == TRUE)
            {
              substr(Seq_minus_char_NNGRRT, start = 26, stop = 26) <- "A"
              
              # patient character string is now searched for suitable PAM sequences:
              
              # Perfect edible window: For Sa-ABEmax  NNGRRT PAM (Editing Position: 3-13)
              # Substring for PAM is selected based on the Editing window:
              
              PerfectUpstreamPAM_NNGRRT <- as.character(substr(Seq_minus_char_NNGRRT, 36,49))  

              # 'minus_perfecteditable_NNGRRT' object is created to store found PAM candidates:
              
              minus_perfecteditable_NNGRRT <- NULL
              
              # The PerfectUpstreamPAM_NNGRRT substring is searched for all different PAM sequences
              # Matches are stored in the 'minus_perfecteditable_NNGRRT' object
              
              for(i in 1:length(PerfectUpstreamPAM_NNGRRT))
              {
                if (grepl( "GAAT", PerfectUpstreamPAM_NNGRRT[i]) == TRUE || grepl( "GAGT", PerfectUpstreamPAM_NNGRRT[i]) == TRUE || grepl( "GGGT", PerfectUpstreamPAM_NNGRRT[i]) == TRUE || grepl( "GGAT", PerfectUpstreamPAM_NNGRRT[i]) == TRUE ){    
                  #print(c(i,PerfectUpstreamPAM_NNGRRT[i]))
                  tmp1 <- c(i,PerfectUpstreamPAM_NNGRRT[i])
                  minus_perfecteditable_NNGRRT <- rbind(minus_perfecteditable_NNGRRT,tmp1)
                  
                  
                }
              }
              
              # Quality check to only run code, if suitable PAM sequences are found:
              
              if(is.null(minus_perfecteditable_NNGRRT) != TRUE) {
                
                minus_perfecteditable_NNGRRTdf <- data.frame(minus_perfecteditable_NNGRRT)
                
                minus_perfecteditable_NNGRRTdf$X2 <- as.character(minus_perfecteditable_NNGRRTdf$X2)
                
                minus_perfecteditable_NNGRRTdf$X1 <- as.numeric(nrow(minus_perfecteditable_NNGRRT))
                
                
                # Guide sequences are designed based on the position of the PAM sequence in the 'minus_perfecteditable_NNGRRT' object:
                
                for (i in 1:nrow(minus_perfecteditable_NNGRRTdf))
                  for (j in 1:(nchar(minus_perfecteditable_NNGRRTdf[1,2])-3))
                  {
                    {
                      if (grepl( "GAAT", substr(minus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE || grepl( "GAGT", substr(minus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE || grepl( "GGGT", substr(minus_perfecteditable_NNGRRTdf[i,2],j,j+3)) == TRUE || grepl( "GGAT", substr(minus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE ){
                        
                        # Guide vector is created with the guide number, Protospacer sequence, Edit Position, PAM sequence, Base editor
                        # Resulting vectors are combined to create 'guides'
                        
                        Guide <- c(row.names(inFile[k,]),substr(Seq_minus_char_NGCG[minus_perfecteditable_NNGRRTdf[i,1]],13+j,32+j),31-(17+j),substr(Seq_minus_char_NGCG[minus_perfecteditable_NNGRRTdf[i,1]],33+j,38+j),"Sa-ABEmax/ABE8e",inFile[k,]$Variant)
                        
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(minus_perfecteditable_NNGRRTdf[1,2])-3))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 2
                          
                        {
                          if(grepl("A", substr(Guide[2],2+l,2+l)) == TRUE)
                          {
                            
                            tmp <- 2+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                      }
                    }
                  }
                
              }
              
              
            }
            
            # If check of wildtype sequence with inFile[k,]'-SNP' does not match, a gets assigned value of NULL
            
            else if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) != TRUE){
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
            }
            
            
          }
          
          # For Base editing guides on the Plus Strand:
          
          else if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "+")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "-"))){
            
            if (("G" == substr(Seq_plus_char_cntrl, start = 26, stop = 26)) == TRUE )
            {
              substr(Seq_plus_char_NNGRRT, start = 26, stop = 26) <- "A"
              
              # For Base editing guides on the Plus Strand:
              
              # Perfect edible window: For Sa-ABEmax  NNGRRT PAM (Editing Position: 3-13)
              # Substring for PAM is selected based on the Editing window:
              
              PerfectDownstreamPAM_NNGRRT <- as.character(substr(Seq_plus_char_NNGRRT, 36, 49))  
              
              ## Test if there is a suitable Downstream PAM sequence
              
              plus_perfecteditable_NNGRRT <- NULL
              
              for(i in 1:length(PerfectDownstreamPAM_NNGRRT))
              {
                if (grepl( "GAAT", PerfectDownstreamPAM_NNGRRT[i]) == TRUE || grepl( "GAGT", PerfectDownstreamPAM_NNGRRT[i]) == TRUE || grepl( "GGGT", PerfectDownstreamPAM_NNGRRT[i]) == TRUE || grepl( "GGAT", PerfectDownstreamPAM_NNGRRT[i]) == TRUE ){    
                  tmp <- c(i,PerfectDownstreamPAM_NNGRRT[i])
                  plus_perfecteditable_NNGRRT <- rbind(plus_perfecteditable_NNGRRT,tmp)
                  
                  
                }
              }
              
              if(is.null(plus_perfecteditable_NNGRRT) != TRUE) {
                
                
                plus_perfecteditable_NNGRRTdf <- data.frame(plus_perfecteditable_NNGRRT)
                
                plus_perfecteditable_NNGRRTdf$X2 <- as.character(plus_perfecteditable_NNGRRTdf$X2)
                
                plus_perfecteditable_NNGRRTdf$X1 <- as.numeric(nrow(plus_perfecteditable_NNGRRT))
                
                
                
                
                for (i in 1:nrow(plus_perfecteditable_NNGRRTdf))
                  for (j in 1:(nchar(plus_perfecteditable_NNGRRTdf[1,2])-3))
                  {
                    {
                      if (grepl( "GAAT", substr(plus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE || grepl( "GAGT", substr(plus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE || grepl( "GGGT", substr(plus_perfecteditable_NNGRRTdf[i,2],j,j+3)) == TRUE || grepl( "GGAT", substr(plus_perfecteditable_NNGRRTdf[i,2], j,j+3)) == TRUE ){
                        
                        Guide <- c(row.names(inFile[k,]),substr(Seq_plus_char_NNGRRT[plus_perfecteditable_NNGRRTdf[i,1]],13+j,32+j),31-(17+j),substr(Seq_plus_char_NNGRRT[plus_perfecteditable_NNGRRTdf[i,1]],33+j,38+j),"Sa-ABEmax/ABE8e",inFile[k,]$Variant)
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(plus_perfecteditable_NNGRRTdf[1,2])-3))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 2
                          
                        {
                          if(grepl("A", substr(Guide[2],2+l,2+l)) == TRUE)
                          {
                            
                            tmp <- 2+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                        
                      }
                    }
                  }
                
              }
            }else{

              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
              
            }
            
          }
          
          #---------------------------------------------#        
          
          # Guide search For SaKKH-ABEmax (NNHRRT PAM) depending on inFile[k,]'-SNP' and '-GeneOrientation':
          # For Base editing guides on the Minus Strand:
          
          if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "-")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "+"))){
            
            # Test wildtype sequence if inFile[k,]'-SNP' was assigned right and if so includes patient mutation into character string:
            # If not, see code line 216
            
            if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) == TRUE)
            {
              substr(Seq_minus_char_NNHRRT, start = 26, stop = 26) <- "A"
              
              # patient character string is now searched for suitable PAM sequences:
              
              # Perfect edible window: For SaKKH-ABEmax  NNHRRT PAM (Editing Position: 4; 6-12)
              # Substring for PAM is selected based on the Editing window:
              
              PerfectUpstreamPAM_NNHRRT <- as.character(substr(Seq_minus_char_NNHRRT, 37,48))  
              
              # 'minus_perfecteditable_NNHRRT' object is created to store found PAM candidates:
              
              minus_perfecteditable_NNHRRT <- NULL
              
              # The PerfectUpstreamPAM_NNHRRT substring is searched for all different PAM sequences
              # Matches are stored in the 'minus_perfecteditable_NNHRRT' object
              
              for(i in 1:length(PerfectUpstreamPAM_NNHRRT))
              {
                if (grepl( "GAAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "GAGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "GGGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "GGAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "TAAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "TAGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "TGGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "TGAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "AAAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "AAGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "AGGT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE || grepl( "AGAT", PerfectUpstreamPAM_NNHRRT[i]) == TRUE ){    
                  tmp1 <- c(i,PerfectUpstreamPAM_NNHRRT[i])
                  minus_perfecteditable_NNHRRT <- rbind(minus_perfecteditable_NNHRRT,tmp1)
                  
                  
                }
              }
              
              # Quality check to only run code, if suitable PAM sequences are found:
              
              if(is.null(minus_perfecteditable_NNHRRT) != TRUE) {
                
                minus_perfecteditable_NNHRRTdf <- data.frame(minus_perfecteditable_NNHRRT)
                
                minus_perfecteditable_NNHRRTdf$X2 <- as.character(minus_perfecteditable_NNHRRTdf$X2)
                
                minus_perfecteditable_NNHRRTdf$X1 <- as.numeric(nrow(minus_perfecteditable_NNHRRT))
                
                
                # Guide sequences are designed based on the position of the PAM sequence in the 'minus_perfecteditable_NNHRRT' object:
                
                for (i in 1:nrow(minus_perfecteditable_NNHRRTdf))
                  for (j in 1:(nchar(minus_perfecteditable_NNHRRTdf[1,2])-3))
                  {
                    {
                      if (grepl("AGGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AAGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AGAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AAAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TGAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl("TGGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TAGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TAAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("GAAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl( "GAGT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl( "GGGT", substr(minus_perfecteditable_NNHRRTdf[i,2],j,j+3)) == TRUE || grepl( "GGAT", substr(minus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ){
                        
                        # Guide vector is created with the guide number, Protospacer sequence, Edit Position, PAM sequence, Base editor
                        # Resulting vectors are combined to create 'guides'
                        
                        if(j != 8){
                          Guide <- c(row.names(inFile[k,]),substr(Seq_minus_char_NNHRRT[minus_perfecteditable_NNHRRTdf[i,1]],14+j,33+j),30-(17+j),substr(Seq_minus_char_NNHRRT[minus_perfecteditable_NNHRRTdf[i,1]],34+j,39+j),"SaKKH-ABEmax/ABE8e",inFile[k,]$Variant)
                          tmp2 <- NULL
                          
                          Download <- rbind(Download,Guide[2])
                          
                          for (l in 1:(nchar(minus_perfecteditable_NNHRRTdf[1,2])-3))
                            
                            ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 3
                            
                          {
                            if(grepl("A", substr(Guide[2],3+l,3+l)) == TRUE)
                            {
                              
                              tmp <- 3+l
                              tmp2 <- rbind(tmp2,tmp)
                              #print("yeha")
                              
                            }else{
                              
                              #print("bla")
                            }
                          }
                          
                          
                          tmp2df <- data.frame(tmp2)
                          tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                          tmp3 <- data.frame(tmp3)
                          tmp4 <- data.frame(tmp3)
                          
                          t <- 0
                          #Guide[6] <- NULL
                          
                          for (b in 1:nrow(tmp3)) {
                            
                            if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                            {
                              Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                              t <- t+43
                              tmp4 <- tmp4[] + 43
                              
                            }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                            {
                              
                              Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                              tmp4 <- tmp4[] + 42
                              
                              
                            }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                            {
                              
                              
                              Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                              tmp4 <- tmp4[] + 43
                              
                            }
                            
                            
                            
                          }
                          
                          
                          guides <- rbind(guides,Guide)
                          
                        }else{
                          
                          guides <- guides
                        }
                        
                        
                      }
                      
                      
                    }
                    
                  }
                
                
              }
              
            }
            # If check of wildtype sequence with inFile[k,]'-SNP' does not match, a gets assigned value of NULL
            
            else if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) != TRUE){
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
              
            }
          }

          # For Base editing guides on the Plus Strand:
          
          else if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "+")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "-"))){
            
            if (("G" == substr(Seq_plus_char_cntrl, start = 26, stop = 26)) == TRUE )
            {
              substr(Seq_plus_char_NNHRRT, start = 26, stop = 26) <- "A"
              
              
              
              # For Base editing guides on the Plus Strand:
              
              # Perfect edible window: For SaKKH-ABEmax  NNHRRT PAM (Editing Position: 4; 6-12)
              # Substring for PAM is selected based on the Editing window:
              
              PerfectDownstreamPAM_NNHRRT <- as.character(substr(Seq_plus_char_NNHRRT, 37, 48))  
              
              ## Test if there is a suitable Downstream PAM sequence
              
              plus_perfecteditable_NNHRRT <- NULL
              
              for(i in 1:length(PerfectDownstreamPAM_NNHRRT))
              {
                if (grepl( "GAAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "GAGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "GGGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "GGAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "TAAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "TAGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "TGGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "TGAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "AAAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "AAGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "AGGT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE || grepl( "AGAT", PerfectDownstreamPAM_NNHRRT[i]) == TRUE ){    
                  #print(c(i,PerfectDownstreamPAM_NNHRRT[i]))
                  tmp <- c(i,PerfectDownstreamPAM_NNHRRT[i])
                  plus_perfecteditable_NNHRRT <- rbind(plus_perfecteditable_NNHRRT,tmp)
                  
                  
                }
              }
              
              if(is.null(plus_perfecteditable_NNHRRT) != TRUE) {
                
                
                plus_perfecteditable_NNHRRTdf <- data.frame(plus_perfecteditable_NNHRRT)
                
                plus_perfecteditable_NNHRRTdf$X2 <- as.character(plus_perfecteditable_NNHRRTdf$X2)
                
                plus_perfecteditable_NNHRRTdf$X1 <- as.numeric(nrow(plus_perfecteditable_NNHRRT))
                
                
                
                
                for (i in 1:nrow(plus_perfecteditable_NNHRRTdf))
                  for (j in 1:(nchar(plus_perfecteditable_NNHRRTdf[1,2])-3))
                  {
                    {
                      if (grepl("AGGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AAGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AGAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("AAAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TGAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl("TGGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TAGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("TAAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ||grepl("GAAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl( "GAGT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE || grepl( "GGGT", substr(plus_perfecteditable_NNHRRTdf[i,2],j,j+3)) == TRUE || grepl( "GGAT", substr(plus_perfecteditable_NNHRRTdf[i,2], j,j+3)) == TRUE ){
                        
                         if(j != 8){ 
                        
                          Guide <- c(row.names(inFile[k,]),substr(Seq_plus_char_NNHRRT[plus_perfecteditable_NNHRRTdf[i,1]],14+j,33+j),30-(17+j),substr(Seq_plus_char_NNHRRT[plus_perfecteditable_NNHRRTdf[i,1]],34+j,39+j),"SaKKH-ABEmax/ABE8e",inFile[k,]$Variant)
                          tmp2 <- NULL
                         
                          
                          Download <- rbind(Download,Guide[2])
                          
                          for (l in 1:(nchar(plus_perfecteditable_NNHRRTdf[1,2])-3))
                            
                            ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 3
                            
                          {
                            if(grepl("A", substr(Guide[2],3+l,3+l)) == TRUE)
                            {
                              
                              tmp <- 3+l
                              tmp2 <- rbind(tmp2,tmp)
                              #print("yeha")
                              
                            }else{
                              
                              #print("bla")
                            }
                          }
                          
                          
                          tmp2df <- data.frame(tmp2)
                          tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                          tmp3 <- data.frame(tmp3)
                          tmp4 <- data.frame(tmp3)
                          
                          t <- 0
                          
                          #Guide[6] <- NULL
                          
                          for (b in 1:nrow(tmp3)) {
                            
                            if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                            {
                              Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                              t <- t+43
                              tmp4 <- tmp4[] + 43
                              
                            }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                            {
                              
                              Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                              tmp4 <- tmp4[] + 42
                              
                              
                            }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                            {
                              
                              
                              Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                              tmp4 <- tmp4[] + 43
                              
                            }
                            
                            
                            
                          }
                          
                          
                          guides <- rbind(guides,Guide)
                          
                         }else{
                           
                           guides <- guides
                         }
                        
                        
                        
                      }
                    }
                  }
                
              }
              
            }else{
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
            }
            
          }
          
          #---------------------------------------------#    
          
          # Guide search For SpG-ABEmax (NGC PAM) depending on inFile[k,]'-SNP' and '-GeneOrientation':
          # For Base editing guides on the Minus Strand:
          
          if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "-")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "+"))){
            
            
            if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) == TRUE )
            {
              substr(Seq_minus_char_NGC, start = 26, stop = 26) <- "A"
              
              
              
              ## Perfect edible: SpG-ABEmax window 5-6 (PAM: NGC)
              
              PerfectUpstreamPAM_NGC <- substr(Seq_minus_char_NGC, 41,44)  
              
              
              ## Target Plus Strand C>T or Minus Strand G>A
              
              minus_perfecteditable_NGC <- NULL;
              
              for(i in 1:length(PerfectUpstreamPAM_NGC))
              {
                if (grepl( "AGC", PerfectUpstreamPAM_NGC[i]) == TRUE || grepl( "GGC", PerfectUpstreamPAM_NGC[i]) == TRUE || grepl( "CGC", PerfectUpstreamPAM_NGC[i]) == TRUE || grepl( "TGC", PerfectUpstreamPAM_NGC[i]) == TRUE ){    
                  tmp <- c(i,PerfectUpstreamPAM_NGC[i])
                  minus_perfecteditable_NGC <- rbind(minus_perfecteditable_NGC,tmp)
                  
                  
                }
              }
              
              if(is.null(minus_perfecteditable_NGC) != TRUE) {
                
                minus_perfecteditabledf_NGC <- data.frame(minus_perfecteditable_NGC)
                
                minus_perfecteditabledf_NGC$X2 <- as.character(minus_perfecteditabledf_NGC$X2)
                
                minus_perfecteditabledf_NGC$X1 <- as.numeric(nrow(minus_perfecteditable_NGC))
                
                
                for (i in 1:nrow(minus_perfecteditabledf_NGC))
                  for (j in 1:(nchar(minus_perfecteditabledf_NGC[1,2])-2))
                  {
                    {
                      if (grepl( "AGC", substr(minus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(minus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE || grepl( "CGC", substr(minus_perfecteditabledf_NGC[i,2],j,j+2)) == TRUE || grepl( "TGC", substr(minus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE ){
                        
                        #print(minus_perfecteditabledf_NGA[i,2])
                        #print(FANCA_char[perfect_editable_sense[i,1]])
                        Guide <- c(row.names(inFile[k,]),substr(Seq_minus_char_NGC[minus_perfecteditabledf_NGC[i,1]],20+j,39+j),26-(19+j),substr(minus_perfecteditabledf_NGC[i,2], j,j+2),"SpG-ABEmax",inFile[k,]$Variant)
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(minus_perfecteditabledf_NGC[1,2])-2))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                          
                        {
                          if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                          {
                            
                            tmp <- 4+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                      }
                    }
                  }
                
              }
            }else{

              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
            }
            
          }
          
          # For Base editing guides on the Plus Strand:
          
          else if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "+")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "-"))){
            
            
            if (("G" == substr(Seq_plus_char_cntrl, start = 26, stop = 26)) == TRUE )
            {
              substr(Seq_plus_char_NGC, start = 26, stop = 26) <- "A"
              
              
              
              ## Perfect edible: SpG-ABEmax window 5-6 (PAM: NGC)
              
              PerfectDownstreamPAM_NGC <- substr(Seq_plus_char_NGA, 41,44) 
              
              ## Target Plus Strand G>A or Minus Strand C>T (with NG-ABEmax)
              
              
              plus_perfecteditable_NGC <- NULL;
              
              for(i in 1:length(PerfectDownstreamPAM_NGC))
              {
                if (grepl( "AGC", PerfectDownstreamPAM_NGC[i]) == TRUE || grepl( "GGC", PerfectDownstreamPAM_NGC[i]) == TRUE || grepl( "CGC", PerfectDownstreamPAM_NGC[i]) == TRUE || grepl( "TGC", PerfectDownstreamPAM_NGC[i]) == TRUE ){    
                  tmp <- c(i,PerfectDownstreamPAM_NGC[i])
                  plus_perfecteditable_NGC <- rbind(plus_perfecteditable_NGC,tmp)
                  
                  
                }
              }
              
              
              if(is.null(plus_perfecteditable_NGC) != TRUE) {
                
                plus_perfecteditabledf_NGC <- data.frame(plus_perfecteditable_NGC)
                
                plus_perfecteditabledf_NGC$X2 <- as.character(plus_perfecteditabledf_NGC$X2)
                
                plus_perfecteditabledf_NGC$X1 <- as.numeric(nrow(plus_perfecteditable_NGC))
                
                
                for (i in 1:nrow(plus_perfecteditabledf_NGC))
                  for (j in 1:(nchar(plus_perfecteditabledf_NGC[1,2])-2))
                  {
                    {
                      if (grepl( "AGC", substr(plus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(plus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE || grepl( "CGC", substr(plus_perfecteditabledf_NGC[i,2],j,j+2)) == TRUE || grepl( "TGC", substr(plus_perfecteditabledf_NGC[i,2], j,j+2)) == TRUE ){
                        
                        #print(plus_perfecteditabledf_NGA[i,2])
                        #print(FANCA_char[perfect_editable_sense[i,1]])
                        Guide <- c(row.names(inFile[k,]),substr(Seq_plus_char_NGC[plus_perfecteditabledf_NGC[i,1]],20+j,39+j),26-(19+j), substr(plus_perfecteditabledf_NGC[i,2], j,j+2),"SpG-ABEmax",inFile[k,]$Variant)
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(plus_perfecteditabledf_NGC[1,2])-2))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                          
                        {
                          if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                          {
                            
                            tmp <- 4+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                      }
                    }
                  }
                
              }
            }else{
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
            }
            
            
          }
          
          #---------------------------------------------#    
          
          # Guide search For SpRY-ABEmax  (NYN PAM) depending on inFile[k,]'-SNP' and '-GeneOrientation':
          # For Base editing guides on the Minus Strand:
          
          if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "-")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "+"))){
            
            
            if (("G" == substr(Seq_minus_char_cntrl, start = 26, stop = 26)) == TRUE )
            {
              substr(Seq_minus_char_NYN, start = 26, stop = 26) <- "A"
              
              
              
              ## Perfect edible: SpRY-ABEmax window 5-6 (PAM: NYN)
              
              PerfectUpstreamPAM_NYN <- substr(Seq_minus_char_NYN, 41,43)  
              
              ## Target Plus Strand C>T or Minus Strand G>A
              
              minus_perfecteditable_NYN <- NULL;
              
              for(i in 1:length(PerfectUpstreamPAM_NYN))
              {
                if (grepl( "C", PerfectUpstreamPAM_NYN[i]) == TRUE || grepl( "T", PerfectUpstreamPAM_NYN[i]) == TRUE || grepl( "A", PerfectUpstreamPAM_NYN[i]) == TRUE){    
                  #print(c(i,PerfectUpstreamPAM_NGA[i]))
                  tmp <- c(i,PerfectUpstreamPAM_NYN[i])
                  minus_perfecteditable_NYN <- rbind(minus_perfecteditable_NYN,tmp)
                  
                  
                }
              }
              
              if(is.null(minus_perfecteditable_NYN) != TRUE) {
                
                minus_perfecteditabledf_NYN <- data.frame(minus_perfecteditable_NYN)
                
                minus_perfecteditabledf_NYN$X2 <- as.character(minus_perfecteditabledf_NYN$X2)
                
                minus_perfecteditabledf_NYN$X1 <- as.numeric(nrow(minus_perfecteditable_NYN))
                
                
                for (i in 1:nrow(minus_perfecteditabledf_NYN))
                  for (j in 1:(nchar(minus_perfecteditabledf_NYN[1,2])-2))
                  {
                    {
                      if (grepl( "C", substr(minus_perfecteditabledf_NYN[i,2], j,j)) == TRUE || grepl( "T", substr(minus_perfecteditabledf_NYN[i,2], j,j)) == TRUE || grepl( "A", substr(minus_perfecteditabledf_NYN[i,2], j,j)) == TRUE){
                        
                        Guide <- c(row.names(inFile[k,]),substr(Seq_minus_char_NYN[minus_perfecteditabledf_NYN[i,1]],20+j,39+j),26-(19+j),substr(Seq_minus_char_NYN[minus_perfecteditabledf_NYN[i,1]],40+j,42+j),"SpRY-ABEmax (might have lower efficiency!)",inFile[k,]$Variant)
                        
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(minus_perfecteditabledf_NYN[1,2])))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                          
                        {
                          if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                          {
                            
                            tmp <- 4+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                      }
                    }
                  }
                
              }
              
              if(is.null(Guide) == TRUE){
                
                guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA","have you thought about using Prime Editing?","There is no Guide for this Variant!","",inFile[k,]$Variant))
                Download <- rbind(Download,"NANANA")
                next
                
                
              }

            }else{
              
              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
            }

          }

          # For Base editing guides on the Plus Strand:
          
          else if(((inFile[k,]$SNP == "G>A") & (inFile[k,]$GeneOrientation == "+")) | ((inFile[k,]$SNP == "C>T") & (inFile[k,]$GeneOrientation == "-"))){
            
            
            if (("G" == substr(Seq_plus_char_cntrl, start = 26, stop = 26)) == TRUE )
            {
              substr(Seq_plus_char_NYN, start = 26, stop = 26) <- "A"
              
              
              
              ## Perfect edible: SpRY-ABEmax window 5-7 (PAM: NYN)
              
              PerfectDownstreamPAM_NYN <- substr(Seq_plus_char_NYN, 41,43) 
              

              ## Target Plus Strand G>A or Minus Strand C>T (with NG-ABEmax)
              
              
              plus_perfecteditable_NYN <- NULL;
              
              for(i in 1:length(PerfectDownstreamPAM_NYN))
              {
                if (grepl( "C", PerfectDownstreamPAM_NGA[i]) == TRUE || grepl( "T", PerfectDownstreamPAM_NGA[i]) == TRUE|| grepl( "A", PerfectDownstreamPAM_NGA[i]) == TRUE){    
                  #print(c(i,PerfectDownstreamPAM_NGA[i]))
                  tmp <- c(i,PerfectDownstreamPAM_NYN[i])
                  plus_perfecteditable_NYN <- rbind(plus_perfecteditable_NYN,tmp)
                  
                  
                }
              }
              
              
              if(is.null(plus_perfecteditable_NYN) != TRUE) {
                
                plus_perfecteditabledf_NYN <- data.frame(plus_perfecteditable_NYN)
                
                plus_perfecteditabledf_NYN$X2 <- as.character(plus_perfecteditabledf_NYN$X2)
                
                plus_perfecteditabledf_NYN$X1 <- as.numeric(nrow(plus_perfecteditable_NYN))
                
                
                for (i in 1:nrow(plus_perfecteditabledf_NYN))
                  for (j in 1:(nchar(plus_perfecteditabledf_NYN[1,2])-2))
                  {
                    {
                      if (grepl( "C", substr(plus_perfecteditabledf_NYN[i,2], j,j)) == TRUE || grepl( "T", substr(plus_perfecteditabledf_NYN[i,2], j,j)) == TRUE || grepl( "A", substr(plus_perfecteditabledf_NYN[i,2], j,j)) == TRUE){
                        
                        Guide <- c(row.names(inFile[k,]),substr(Seq_plus_char_NYN[plus_perfecteditabledf_NYN[i,1]],20+j,39+j),26-(19+j),substr(Seq_plus_char_NYN[plus_perfecteditabledf_NYN[i,1]],40+j,42+j),"SpRY-ABEmax (might have lower efficiency!)",inFile[k,]$Variant)
                        tmp2 <- NULL
                        
                        Download <- rbind(Download,Guide[2])
                        
                        for (l in 1:(nchar(plus_perfecteditabledf_NYN[1,2])))
                          
                          ## Make substring of the specific nucleotide (lowest Edit position -1) <- here 4
                          
                        {
                          if(grepl("A", substr(Guide[2],4+l,4+l)) == TRUE)
                          {
                            
                            tmp <- 4+l
                            tmp2 <- rbind(tmp2,tmp)
                            #print("yeha")
                            
                          }else{
                            
                            #print("bla")
                          }
                        }
                        
                        
                        tmp2df <- data.frame(tmp2)
                        tmp3 <- tmp2df[order(tmp2df$tmp2,decreasing =  FALSE),]
                        tmp3 <- data.frame(tmp3)
                        tmp4 <- data.frame(tmp3)
                        
                        t <- 0
                        #Guide[6] <- NULL
                        
                        for (b in 1:nrow(tmp3)) {
                          
                          if(as.numeric(tmp3[b,]) < as.numeric(Guide[3]))
                          {
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            t <- t+43
                            tmp4 <- tmp4[] + 43
                            
                          }else if(as.numeric(tmp3[b,]) == as.numeric(Guide[3]))
                          {
                            
                            Guide[2] <- paste(substr(Guide[2],1,as.numeric(as.numeric(Guide[3])+t-1)),"<strong><font color='red'>A</font></strong>",substring(Guide[2], first = as.numeric(as.numeric(Guide[3])+t+1)), sep = "")
                            tmp4 <- tmp4[] + 42
                            
                            
                          }else if (as.numeric(tmp3[b,]) > as.numeric(Guide[3]))
                          {
                            
                            
                            Guide[2] <- paste(substr(Guide[2],1,tmp4[b,]-1),"<strong><font color='blue'>A</font></strong>",substring(Guide[2], first = tmp4[b,]+1),sep = "")
                            tmp4 <- tmp4[] + 43
                            
                          }
                          
                          
                          
                        }
                        
                        
                        guides <- rbind(guides,Guide)
                        
                        
                        
                        
                      }
                    }
                  }
                
              }
              
              if(is.null(Guide) == TRUE){
                
                guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA","have you thought about using Prime Editing?","There is no Guide for this Variant!","",inFile[k,]$Variant))
                Download <- rbind(Download,"NANANA")
                next
                
                
              }
       
              
            }else{

              guides <- rbind(guides,c(row.names(inFile[k,]),"NANANA", "Could not find target sequence, check substitution!","","",inFile[k,]$Variant))
              Download <- rbind(Download,"NANANA")
              next
            }
            
            
          }

          #---------------------------------------------#   
 
            
          
          
          }
    
          
          
          currentguidesdf <- reactive({guidesdf <- data.frame("Variant" = guides[,6],"Protospacer" = guides[,2], "EditPos." = guides[,3], "PAM" = guides[,4], "Base Editor" = guides[,5] )})


          output$mytable  <- DT::renderDataTable({ currentguidesdf()
                                                  },
                                                  escape = FALSE,
                                                  callback=JS(
                                                  'table.on("click.dt","td", function() {
                                                    var data=table.cell(this).data();
                                                    if (table.cell(this).data() > 0)
                                                    swal({title: "Edit position", text:"Edit position is defined relative to the PAM sequence. The protospacer position closest to the PAM is position +20, the furthest away +1. (Figure modified from Gaudelli et.al., 2017)",imageUrl: "Base.jpg",imageSize: "460x200"
                                                    });
                                                 })'
            
            
            
            
            
          ),
          )
     
          
          
         
          
          output$myText <- renderText({

            
            total <- data.frame(inFile$Variant)
            targetable <- data.frame("Baseeditor" = guides[,5], "Variant" = guides[,6], stringsAsFactors = FALSE)
            targetable <- subset(targetable, nchar(targetable$Baseeditor) > 3)
            targetable <- data.frame(unique(targetable$Variant))
            
            
            x <- nrow(targetable)
            y <- nrow(total)
            
            
            paste(x, "out of", y, "Variants could be targeted")
            
            
          })
          

        }
        
        else if (input$Editing == "Prime editing") 
          
        {
          
          output$downloadData <- downloadHandler(
            filename = function() {
              paste('pegRNA oligos for cloning-', Sys.Date(), '.csv', sep='')
            },
            content = function(con) {
              write.csv(data.frame("Variant" = inFile[as.numeric(guides[,1]),]$Variant, "Protospacer(Sense)" = paste("cacc", guides[,2],"gtttt", sep = ""),"Protospacer(Antisense)" = paste("ctctaaaac",reverseComplement(DNAStringSet(guides[,2])),sep = "") , "EditPos." = guides[,3], "Extension(Sense)" = paste("gtgc",guides[,4],sep = ""), "Extension(Antisense)" = paste("aaaa",reverseComplement(DNAStringSet(guides[,4])), sep = ""), "PAM" = guides[,5], "PAM-Strand" = guides[,6], "Score" = as.numeric(guides[,7]), "PBS" =  inFile[as.numeric(guides[,1]),]$PBS, "RTT" = inFile[as.numeric(guides[,1]),]$RTT, stringsAsFactors = FALSE), con, sep = ",", row.names=FALSE)
            }
          )
          
      
          guides <- NULL
          
          for (k in 1:nrow(inFile))
          {
            a <- 1
            b <- 1
            x <- 0
            z <- 0
            y <- 0
            u <- 0
            m <- 0
            n <- 0
            p <- 0
            v <- 0
            w <- 0
            q <- 0
            
            NoGuide <- NULL
            
            if(input$Genomes != "Rice (MSU7)")
            {
              
              Chromosome <- as.character(paste("chr",inFile[k,]$Chromosome, sep = ""))
              
              
            }
            else if(input$Genomes == "Rice (MSU7)")
            {
              
              Chromosome <- as.character(paste("Chr",inFile[k,]$Chromosome, sep = ""))
              
            }
            
            inFile$RT <- inFile$RTT
            
            incProgress(1/nrow(inFile), detail = paste("Variant", k))  
            
            
            if((inFile[k,]$GeneOrientation == "-")== TRUE)
            {
              
              if(is.na(seqlengths(checkCompatibleSeqinfo(Target_Genome, GRanges(as.character(paste("chr",inFile[k,]$Chromosome, sep = "")),ranges = IRanges(start = as.numeric(inFile[k,]$GenomicLocation)+as.numeric(paste(inFile[k,]$GeneOrientation,"75", sep = "")), width = 151)))[c(as.character(Chromosome))])) != TRUE)
              {
                
              ### Get sequence for Plus Strand
                
              
              Seq_plus <- getSeq(Target_Genome, GRanges(as.character(paste("chr",inFile[k,]$Chromosome, sep = "")),ranges = IRanges(start = as.numeric(inFile[k,]$GenomicLocation)+as.numeric(paste(inFile[k,]$GeneOrientation,"75", sep = "")), width = 151)),)
              
              ### Get Sequence for Minus Strand
              
              Seq_minus <- reverseComplement(Seq_plus)
              
              Seq_plus_rev <- reverse(Seq_plus)
              Seq_minus_rev <- reverse(Seq_minus)
              
              Seq_minus_char <- as.character(Seq_minus)
              Seq_plus_char <- as.character(Seq_plus)
              Seq_plus_rev_char <- as.character(Seq_plus_rev)
              Seq_minus_rev_char <- as.character(Seq_minus_rev)
              
              
              if ((substr(inFile[k,]$Edit,1,3) == "ins")==TRUE){
              
                DNAchange <- DNAStringSet(substring(inFile[k,]$Edit, first = 4))
                anti_DNAchange <- complement(DNAchange)
                
                DNAchange_char <- as.character(DNAchange)
                
                DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
                
                anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
                
                anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
                
                
                x <- nchar(as.character(DNAchange))
                
                u <- nchar(as.character(DNAchange))
              
                #x <- as.numeric(x)
                
                p <- 1
              
                Seq_plus_rev_char_n <- paste(substr(Seq_plus_rev_char, start = 1, stop = 76),as.character(anti_DNAchange), substr(Seq_plus_rev_char, start = 77, stop = 151), sep = "")
                Seq_minus_char_n <- paste(substr(Seq_minus_char, start = 1, stop = 75),as.character(DNAchange), substr(Seq_minus_char, start = 76, stop = 151), sep = "")
             
                
              }
              else if ((substr(inFile[k,]$Edit,1,3) == "del")==TRUE){
              
               DNAchange <- DNAStringSet(substring(inFile[k,]$Edit, first = 4))
               anti_DNAchange <- complement(DNAchange)
               
               DNAchange_char <- paste("[",as.character(DNAchange),"]" ,sep = "")
               
               DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
               
               anti_DNAchange_char <- paste("[",as.character(reverse(anti_DNAchange)),"]" ,sep = "")
               
               anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
               
              
               x <- nchar(as.character(DNAchange))
               
               z <- nchar(as.character(DNAchange))
              
               x <- as.numeric(x)* (-1)  
               
               m <- z-1
               
               q <- 1
              
              
                if(((substr(Seq_minus, start = 76, stop = 75-x)) == (as.character(substring(inFile[k,]$Edit, first = 4))))==TRUE){
                
                  Seq_plus_rev_char_n <- paste(substr(Seq_plus_rev_char, start = 1, stop = 75), substr(Seq_plus_rev_char, start = 76-x, stop = 151), sep = "")
                  Seq_minus_char_n <- paste(substr(Seq_minus_char, start = 1, stop = 75), substr(Seq_minus, start = 76-x, stop = 151), sep = "")
                
                }else{
                
                
                  NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","CHECK DEL!",-999,"")
                  guides <- rbind(guides, NoGuide)
                  next
                
                }
              
              
            }
              else if(grepl(">", inFile[k,]$Edit) == TRUE) { 
              
              
              y <- as.numeric(((nchar(inFile[k,]$Edit))-1)/2)
              n <- y-1
              
              
              DNAchange <- DNAStringSet(substr(inFile[k,]$Edit,2+y,2+(y*2)))
              anti_DNAchange <- complement(DNAchange)
              
              DNAchange_char <- as.character(DNAchange)
              
              DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              anti_DNAchange_char <- as.character(anti_DNAchange)
              
              anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
              
              
                if ((substr(inFile[k,]$Edit,0,y) == substr(Seq_minus, start = 76, stop = 75+y)) == TRUE ){
                
                  Seq_plus_rev_char_n <- Seq_plus_rev_char
                  Seq_minus_char_n <- Seq_minus_char
                
                  substr(Seq_plus_rev_char_n, start = 76, stop = 75+y) <- as.character(anti_DNAchange)
                
                  substr(Seq_minus_char_n, start = 76, stop = 75+y) <- as.character(DNAchange)
                
                }else{
                
                  NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAA","no PAM","CHECK SUB!",-999,"")
                  guides <- rbind(guides, NoGuide)
                  next
                }
            }

        
            ## Searchig for sense PAM in Prime editing 
            
            
            
            
            ## Upstream PAM for PRIME 
            
            primeUpstreamPAM <- substr(Seq_minus_char, 5,81)  
            
            primeUpstreamPAM2 <- substr(Seq_plus_rev_char,71,147)
            
            
            ## prime edits possible:
            
            
            ## Sense strand
            
            primeeditable <- NULL;
            
            
            for(i in 1:length(primeUpstreamPAM))
            {
              if (grepl( "AGG", primeUpstreamPAM[i]) == TRUE || grepl( "CGG", primeUpstreamPAM[i]) == TRUE || grepl( "GGG", primeUpstreamPAM[i]) == TRUE || grepl( "TGG", primeUpstreamPAM[i]) == TRUE ){
                
                
                temp <- c(i,primeUpstreamPAM[i])
                primeeditable <- rbind(primeeditable, temp)
                
                
              }
            }
            
            
            primeeditabledf <- data.frame(primeeditable)
            
            primeeditabledf$X2 <- as.character(primeeditabledf$X2)
            
            primeeditabledf$X1 <- as.numeric(primeeditabledf$X1)
            
            
            
            
            ## AntiSense strand
            
            primeeditable2 <- NULL;
            
            
            for(i in 1:length(primeUpstreamPAM2))
            {
              if (grepl( "GGA", primeUpstreamPAM2[i]) == TRUE || grepl( "GGC", primeUpstreamPAM2[i]) == TRUE || grepl( "GGG", primeUpstreamPAM2[i]) == TRUE || grepl( "GGT", primeUpstreamPAM2[i]) == TRUE ){
                
                
                temp <- c(i,primeUpstreamPAM2[i])
                primeeditable2 <- rbind(primeeditable2, temp)
                
                
              }
            }
            
            
            primeeditable2df <- data.frame(primeeditable2)
            
            primeeditable2df$X2 <- as.character(primeeditable2df$X2)
            
            primeeditable2df$X1 <- as.numeric(primeeditable2df$X1)
            
            
            
            
            if(is.null(primeeditable) != TRUE){
              
              
              ## Sense strand
              
              
              
              for (i in 1:nrow(primeeditabledf))
                for (j in 1:(nchar(primeeditabledf[1,2])-2))
                {
                  {
                    if (grepl( "AGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(primeeditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE ){
                      
                      
                      if(nchar(substr(Seq_minus_char_n[primeeditabledf[i,1]],(j-16),j+3)) == 20)
                      {
                        if((76-j+u+y) <= as.numeric(inFile[k,]$RT))
                        {
 
                        Guide <- c(row.names(inFile[k,]),substr(Seq_minus_char[primeeditabledf[i,1]],j-16,j+3),76-j,as.character(reverse(DNAStringSet(substr(Seq_plus_rev_char_n[primeeditabledf[i,1]],(j+1) - as.numeric(inFile[k,]$PBS),((j+1) - (as.numeric(inFile[k,]$PBS)))+ (as.numeric(inFile[k,]$RT)+as.numeric(inFile[k,]$PBS)-1))))),substr(primeeditabledf[i,2], j,j+2), "Sense")
                        
                        if(substr(Seq_minus_char_n[primeeditabledf[i,1]],j-16,j-16) != "G")
                        {
                          
                          Guide[2] <- paste("g",Guide[2], sep = "") 
                          
                        }
                        
                        Guide[7] <- 0
                        
                        
                        
                        if((substr(Guide[4],1,1) == "C") == TRUE){
                          
                          Guide[7] <- as.numeric(Guide[7])-28
                          
                        }
                        
                        if ((grepl("TTTTT", Guide[4])) == TRUE){
                          
                          Guide[7] <- as.numeric(Guide[7])-50
                          
                        }
                        
                        if (((as.numeric(inFile[k,]$RT) - (as.numeric(Guide[3])+y+u)) >= 4) != TRUE) 
                        {
                          
                          Guide[7] <- as.numeric(Guide[7])-6
                          
                        }
                        
                        
                        
                        Guide[7] <- as.numeric(Guide[7]) - as.numeric(Guide[3])
                        
                        Guide[8] <- paste(reverseComplement(DNAStringSet(substring(Guide[4], first = 1+y+q+((as.numeric(inFile[k,]$RT)-(as.numeric(Guide[3]))))))),as.character(DNAchange_char_tagged),reverseComplement(DNAStringSet(substring(Guide[4],1,as.numeric(inFile[k,]$RT)-((as.numeric(Guide[3])+u-q))))), sep = "")

                        guides <- rbind(guides,Guide)
                        
                        
                        }else if ((76-j+u+y) > as.numeric(inFile[k,]$RT))
                        {
                          
                          NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","Edit too far away, try to increase the RT length!",-999,"")
                          
                          next
                          
                        }
                      }
                      
                    }
                  }
                }
            }
            
            if(is.null(primeeditable2) != TRUE){
              
              ## AntiSense strand
              
              for (i in 1:nrow(primeeditable2df))
                for (j in 1:(nchar(primeeditable2df[1,2])-2))
                {
                  {
                    if (grepl( "GGA", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(primeeditable2df[i,2],j,j+2)) == TRUE || grepl( "GGT", substr(primeeditable2df[i,2], j,j+2)) == TRUE ){
                      
                      if((j+u+y) <= as.numeric(inFile[k,]$RT))
                      {
                        if((j-z) >= 0)
                        {
                      
                         Guide <- c(row.names(inFile[k,]),substr(Seq_plus_rev_char[primeeditable2df[i,1]],j+73,j+92),j, as.character(reverse(DNAStringSet(substr(Seq_minus_char_n[primeeditable2df[i,1]],((j+76+x)+(as.numeric(inFile[k,]$PBS)))-(as.numeric(inFile[k,]$RT)+as.numeric(inFile[k,]$PBS)),((j+76+x)+(as.numeric(inFile[k,]$PBS)))-1)))),as.character(reverse(DNAStringSet(substr(primeeditable2df[i,2], j,j+2)))), "Antisense")

                         Guide[2] <- as.character(reverse(DNAStringSet(Guide[2])))
                         
                         Guide[4] <- as.character(reverse(DNAStringSet(Guide[4])))
                         
                         
                         if((substr(Seq_plus_rev_char[primeeditable2df[i,1]],j+73,j+73) != "G") == TRUE)
                         {
                           
                           Guide[2] <- paste("g",Guide[2], sep = "") 
                           
                         }
                  
                        Guide[7] <- 0
                        
                        
                        if((substr(Guide[4],1,1) == "C") == TRUE){
                          
                          Guide[7] <- as.numeric(Guide[7])-28
                          
                        }
                        
                        if ((grepl("TTTTT", Guide[4])) == TRUE){
                          
                          Guide[7] <- as.numeric(Guide[7])-50
                          
                        }
                        
                        if (((as.numeric(inFile[k,]$RT) - (as.numeric(Guide[3])+y+u)) >= 4) != TRUE) 
                        {
                          
                          Guide[7] <- as.numeric(Guide[7])-6
                          
                        }
                        
                        
                        Guide[7] <- as.numeric(Guide[7]) - as.numeric(Guide[3])
                        Guide[8] <- paste(substring(Guide[4],1,as.numeric(inFile[k,]$RT)-((as.numeric(Guide[3])+x))),as.character(DNAchange_char_tagged), substring(Guide[4], first = 1+z+y+((as.numeric(inFile[k,]$RT)-(as.numeric(Guide[3]))))), sep = "")
                        guides <- rbind(guides,Guide)
                        
                        
                        }
                      
                      else if ((j-z) < 0)
                      {
                        NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","No Prime Editing guide found!",-999,"")
                        
                        next
                      }

                        
                    }else if ((j+u+y) > as.numeric(inFile[k,]$RT))
                      {
                      
                      NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","Edit too far away, try to increase the RT length!",-999,"")
                      
                      
                      next
                      
                    }
                      
                    }
                }
              
              
            }


            }
            }
              else{
              
              NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","CHECK GENOME COORDINATES!",-999,"")
              
              next
              
            }
            }
            if((inFile[k,]$GeneOrientation == "+") == TRUE)
            {
              
              if(is.na(seqlengths(checkCompatibleSeqinfo(Target_Genome, GRanges(as.character(paste("chr",inFile[k,]$Chromosome, sep = "")),ranges = IRanges(start = as.numeric(inFile[k,]$GenomicLocation)+as.numeric(paste(inFile[k,]$GeneOrientation,"75", sep = "")), width = 151)))[c(as.character(Chromosome))])) != TRUE)
              {
              ### Get sequence for Plus Strand
              
              Seq_plus <- getSeq(Target_Genome, GRanges(as.character(paste("chr",inFile[k,]$Chromosome, sep = "")),ranges = IRanges(start = as.numeric(inFile[k,]$GenomicLocation)-as.numeric(paste(inFile[k,]$GeneOrientation,"75", sep = "")), width = 151)),)
              
              ### Get Sequence for Minus Strand
              
              Seq_minus <- reverseComplement(Seq_plus)
              
              Seq_plus_rev <- reverse(Seq_plus)
              Seq_minus_rev <- reverse(Seq_minus)
              
              Seq_minus_char <- as.character(Seq_minus)
              Seq_plus_char <- as.character(Seq_plus)
              Seq_plus_rev_char <- as.character(Seq_plus_rev)
              Seq_minus_rev_char <- as.character(Seq_minus_rev)

              if ((substr(inFile[k,]$Edit,1,3) == "ins")==TRUE){
                
                DNAchange <- DNAStringSet(substring(inFile[k,]$Edit, first = 4))
                anti_DNAchange <- complement(DNAchange)
                
                
                DNAchange_char <- as.character(DNAchange)
                
                DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
                
                anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
                
                anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
                
                p <- 1
                
                x <- nchar(as.character(DNAchange))
                u <- nchar(as.character(DNAchange))
                
                x <- as.numeric(x)
                
                Seq_plus_char_n <- paste(substr(Seq_plus_char, start = 1, stop = 76),as.character(DNAchange), substr(Seq_plus_char, start = 77, stop = 151), sep = "")
                Seq_minus_rev_char_n <- paste(substr(Seq_minus_rev_char, start = 1, stop = 75),as.character(anti_DNAchange), substr(Seq_minus_rev_char, start = 76, stop = 151), sep = "")
                
              }
              else if ((substr(inFile[k,]$Edit,1,3) == "del")==TRUE){
                
                DNAchange <- DNAStringSet(substring(inFile[k,]$Edit, first = 4))
                anti_DNAchange <- complement(DNAchange)
        
                DNAchange_char <- paste("[",as.character(DNAchange),"]" ,sep = "")
                
                DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
                
                anti_DNAchange_char <- paste("[",as.character(reverse(anti_DNAchange)),"]" ,sep = "")
                
                anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
                
                z <- nchar(as.character(DNAchange))
                
                x <- nchar(as.character(DNAchange))
                
                x <- as.numeric(x)* (-1)  
                
                m <- z-1
                
                w <- 1
                
                
                if(((substr(Seq_plus, start = 76, stop = 75+z)) == (as.character(DNAchange)))==TRUE){
                  
                  Seq_plus_char_n <- paste(substr(Seq_plus_char, start = 1, stop = 75), substr(Seq_plus_char, start = 76+z, stop = 151), sep = "")
                  
                  Seq_minus_rev_char_n <- paste(substr(Seq_minus_rev_char, start = 1, stop = 75), substr(Seq_minus_rev_char, start = 76+z, stop = 151), sep = "")
                  
                  
                }else{
                  
                  
                  NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","CHECK DEL!",-999,"")
                  guides <- rbind(guides, NoGuide)
                  next
                  
                }
                
                
              }
              else if(grepl(">", inFile[k,]$Edit) == TRUE) { 
                
                
                y <- 0
                y <- ((nchar(inFile[k,]$Edit))-1)/2
                v <- 1
                
                
                DNAchange <- DNAStringSet(substr(inFile[k,]$Edit,2+y,2+(y*2)))
                anti_DNAchange <- complement(DNAchange)
                
                DNAchange_char <- as.character(DNAchange)
                
                DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", DNAchange_char,.noWS = "outside"),.noWS = "outside")
                
                anti_DNAchange_char <- as.character(reverse(anti_DNAchange))
                
                anti_DNAchange_char_tagged <- tags$strong(tags$span(style="color:red", anti_DNAchange_char,.noWS = "outside"),.noWS = "outside")
 
                if ((substr(inFile[k,]$Edit,0,y) == substr(Seq_plus, start = 76, stop = 75+y)) == TRUE ){
                  
                  Seq_plus_char_n <- Seq_plus_char
                  Seq_minus_rev_char_n <- Seq_minus_rev_char
                  
                  substr(Seq_plus_char_n, start = 76, stop = 75+y) <- as.character(DNAchange)
                  
                  substr(Seq_minus_rev_char_n, start = 76, stop = 75+y) <- as.character(anti_DNAchange)
                  
                }else{
                  
                  NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAA","no PAM","CHECK SUB!",-999,"")
                  guides <- rbind(guides, NoGuide)
                  next
                }
              }
              
              
              ## Searchig for sense PAM in Prime editing 
              
              
              
              
              ## Upstream PAM for PRIME 
              
              primeUpstreamPAM <- substr(Seq_plus_char, 5,81)  
              
              primeUpstreamPAM2 <- substr(Seq_minus_rev_char,71,147)

              ## prime edits possible:
              
              
              ## Sense strand
              
              primeeditable <- NULL;
              
              
              for(i in 1:length(primeUpstreamPAM))
              {
                if (grepl( "AGG", primeUpstreamPAM[i]) == TRUE || grepl( "CGG", primeUpstreamPAM[i]) == TRUE || grepl( "GGG", primeUpstreamPAM[i]) == TRUE || grepl( "TGG", primeUpstreamPAM[i]) == TRUE ){
                  
                  
                  temp <- c(i,primeUpstreamPAM[i])
                  primeeditable <- rbind(primeeditable, temp)
                  
                  
                }
              }
              
              
              primeeditabledf <- data.frame(primeeditable)
              
              primeeditabledf$X2 <- as.character(primeeditabledf$X2)
              
              primeeditabledf$X1 <- as.numeric(primeeditabledf$X1)
              
              
              
              
              ## AntiSense strand
              
              primeeditable2 <- NULL;
              
              
              for(i in 1:length(primeUpstreamPAM2))
              {
                if (grepl( "GGA", primeUpstreamPAM2[i]) == TRUE || grepl( "GGC", primeUpstreamPAM2[i]) == TRUE || grepl( "GGG", primeUpstreamPAM2[i]) == TRUE || grepl( "GGT", primeUpstreamPAM2[i]) == TRUE ){
                  
                  
                  temp <- c(i,primeUpstreamPAM2[i])
                  primeeditable2 <- rbind(primeeditable2, temp)
                  
                  
                }
              }
              
              
              primeeditable2df <- data.frame(primeeditable2)
              
              primeeditable2df$X2 <- as.character(primeeditable2df$X2)
              
              primeeditable2df$X1 <- as.numeric(primeeditable2df$X1)
          
              
              if(is.null(primeeditable) != TRUE){
                
                
                ## Sense strand
             
                for (i in 1:nrow(primeeditabledf))
                  for (j in 1:(nchar(primeeditabledf[1,2])-2))
                  {
                    {
                      if (grepl( "AGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE || grepl( "CGG", substr(primeeditabledf[i,2],j,j+2)) == TRUE || grepl( "TGG", substr(primeeditabledf[i,2], j,j+2)) == TRUE ){
                        
                        if((nchar(substr(Seq_plus_char[primeeditabledf[i,1]],(j-16),j+3)) == 20))
                        {
                          if(as.numeric(76-j)+u-p+y <= as.numeric(inFile[k,]$RT))
                          {

                          Guide <- c(row.names(inFile[k,]),substr(Seq_plus_char[primeeditabledf[i,1]],j-16,j+3),76-j,as.character(reverse(DNAStringSet(substr(Seq_minus_rev_char_n[primeeditabledf[i,1]],(j+1) - as.numeric(inFile[k,]$PBS),((j+1) - (as.numeric(inFile[k,]$PBS)))+ (as.numeric(inFile[k,]$RT)+as.numeric(inFile[k,]$PBS)-1))))),substr(primeeditabledf[i,2], j,j+2), "Sense")
                          
                          if(substr(Seq_plus_char[primeeditabledf[i,1]],j-16,j-16) != "G")
                          {
                            
                            Guide[2] <- paste("g",Guide[2], sep = "") 
                            
                          }
                          
                          Guide[7] <- 0
                          
                          
                          if((substr(Guide[4],1,1) == "C") == TRUE){
                            
                            Guide[7] <- as.numeric(Guide[7])-28
                            
                          }
                          
                          if ((grepl("TTTTT", Guide[4])) == TRUE){
                            
                            Guide[7] <- as.numeric(Guide[7])-50
                            
                          }
                          
                          if (((as.numeric(inFile[k,]$RT) - (as.numeric(Guide[3])+y+u)) >= 4) != TRUE) 
                          {
                            
                            Guide[7] <- as.numeric(Guide[7])-6
                            
                          }
                          
                          
                          
                          
                          Guide[7] <- as.numeric(Guide[7])-as.numeric(Guide[3])
                          
                          
                          
                          Guide[8] <- paste(reverseComplement(DNAStringSet(substring(Guide[4], first = 2+((as.numeric(inFile[k,]$RT)-(as.numeric(Guide[3]))))))),as.character(DNAchange_char_tagged), reverseComplement(DNAStringSet(substring(Guide[4],1,1+as.numeric(inFile[k,]$RT)-((as.numeric(Guide[3])+y+x+z))))), sep = "")

                          
                          guides <- rbind(guides,Guide)
                          
                          }else if(as.numeric(76-j)+u-p+y > as.numeric(inFile[k,]$RT))
                          {
                            NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","Edit too far away, try to increase the RT length!",-999,"")
                            
                            next
                            
                          }
                        }
                        
                      } 
                    }
                  }
              }
              
              if(is.null(primeeditable2) != TRUE){
                
                ## AntiSense strand
                
                for (i in 1:nrow(primeeditable2df))
                  for (j in 1:(nchar(primeeditable2df[1,2])-2))
                  {
                    {
                      if (grepl( "GGA", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGG", substr(primeeditable2df[i,2], j,j+2)) == TRUE || grepl( "GGC", substr(primeeditable2df[i,2],j,j+2)) == TRUE || grepl( "GGT", substr(primeeditable2df[i,2], j,j+2)) == TRUE ){
                        
                        
                        if(nchar(substr(Seq_minus_rev_char[primeeditable2df[i,1]],j+73,j+92)) == 20){
                          
                          if((as.numeric(j)+u-p+y <= as.numeric(inFile[k,]$RT)))
                          {
                            
                          if((as.numeric(j)-z) >= 0)
                          {
             
                          Guide <- c(row.names(inFile[k,]),substr(Seq_minus_rev_char[primeeditable2df[i,1]],j+73,j+92),j,as.character(reverse(DNAStringSet(substr(Seq_plus_char_n[primeeditable2df[i,1]],((j+76+x)+(as.numeric(inFile[k,]$PBS)))-(as.numeric(inFile[k,]$RT)+as.numeric(inFile[k,]$PBS)),((j+76+x)+(as.numeric(inFile[k,]$PBS)))-1)))),as.character(reverse(DNAStringSet(substr(primeeditable2df[i,2], j,j+2)))), "Antisense")

                          
                          
                          Guide[2] <- as.character(reverse(DNAStringSet(Guide[2])))
                          
                          Guide[4] <- as.character(reverse(DNAStringSet(Guide[4])))
                          
                          
                          if((substr(Seq_minus_rev_char[primeeditable2df[i,1]],j+92,j+92) != "G") == TRUE)
                          {
                            
                            Guide[2] <- paste("g",Guide[2], sep = "") 
                            
                          }

                          Guide[7] <- 0
                          
                          if((substr(Guide[4],1,1) == "C") == TRUE){
                            
                            Guide[7] <- as.numeric(Guide[7])-28
                            
                          }
                          
                          if ((grepl("TTTTT", Guide[4])) == TRUE){
                            
                            Guide[7] <- as.numeric(Guide[7])-50
                            
                          }
                          
                          if (((as.numeric(inFile[k,]$RT) - (as.numeric(Guide[3])+y+u)) >= 4) != TRUE) 
                          {
                            
                            Guide[7] <- as.numeric(Guide[7])-6
                            
                          }
                          
                          
                          
                          Guide[7] <- as.numeric(Guide[7])-as.numeric(Guide[3])
                          
                          Guide[8] <-paste(substring(Guide[4],1,1+as.numeric(inFile[k,]$RT)-((as.numeric(Guide[3])+v+x+z+p-m))),as.character(DNAchange_char_tagged), substring(Guide[4], first = 2+((as.numeric(inFile[k,]$RT)-(as.numeric(Guide[3])+p-m+v-y)))), sep = "")

                          
                          guides <- rbind(guides,Guide)
                          
                            }else if((as.numeric(j)-z) < 0)
                            {
                              NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","Edit too far away, try to increase the RT length!",-999,"")
                              
                              
                              next
                              
                              
                            }
                          
                          }else if(as.numeric(j)+u-p+y > inFile[k,]$RT)
                          {
                            NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","Edit too far away, try to increase the RT length!",-999,"")
                            
                            next
                            
                            
                          }
                        }
                        
                      }
                    }
                  }
                
                
              }
 
              }
              else{
                
                NoGuide <- c(row.names(inFile[k,]),as.character("AAAAAAAAAAAAAAAAA"),0,"AAAAA","no PAM","CHECK GENOME COORDINATES!",-999,"")
                
                next
                
              }
              
              
            }
            
            
            guides <- rbind(guides,NoGuide)
            
          }

      
          ## Change Orientation of Protospacer and 3' extension

          
          currentguidesdf <- reactive({ 
            
            if (input$expert != 1) {
              
              guidesdf <- data.frame("Variant" = inFile[as.numeric(guides[,1]),]$Variant, "Protospacer(Sense)" = guides[,2],"Protospacer(Antisense)" = as.character(reverseComplement(DNAStringSet(guides[,2]))), "EditPos." = guides[,3], "Extension(coding strand)" = guides[,8], "PAM" = guides[,5], "PAM-Strand" = guides[,6], "Score" = as.numeric(guides[,7]), stringsAsFactors = FALSE)
              
              df.agg <- aggregate(Score ~ Variant, guidesdf, max)
              guidesdf$id  <- 1:nrow(guidesdf)
              df.max <- join(df.agg, guidesdf)
              df.max <- df.max[!duplicated(df.max[,"Variant"]),]
              df.max[order(df.max$id), ]
              
            }
            else{
              
              guidesdf <- data.frame("Variant" = inFile[as.numeric(guides[,1]),]$Variant, "Protospacer(Sense)" = guides[,2],"Protospacer(Antisense)" = as.character(reverseComplement(DNAStringSet(guides[,2]))), "EditPos." = guides[,3], "Extension(Sense)" = guides[,8], "PAM" = guides[,5], "PAM-Strand" = guides[,6], "Score" = as.numeric(guides[,7]), stringsAsFactors = FALSE)
              
              
            }
            
          }) 
          
          
          #currentguidesdf2 <- reactive({
          #  
          #  guidesdf2 <- data.frame("Variant" = inFile[as.numeric(guides[,1]),]$Variant, "Protospacer(Sense)" = paste("cacc", guides[,2],"gtttt", sep = ""),"Protospacer(Antisense)" = paste("ctctaaaac",reverseComplement(DNAStringSet(guides[,2])),sep = "") , "TargetPos." = guides[,3], "Extension(Sense)" = paste("gtgc",guides[,4],sep = ""), "Extension(Antisense)" = paste("aaaa",reverseComplement(DNAStringSet(guides[,4])), sep = ""), "PAM" = guides[,5], "PAM-Strand" = guides[,6], "Score" = as.numeric(guides[,7]), stringsAsFactors = FALSE)
          #  
          #  
          #})
          #
          #output$mytable2 <- DT::renderDataTable ({ currentguidesddf2()})
          
          output$mytable  <- DT::renderDataTable({ currentguidesdf()}, 
                                                 escape = FALSE,
                                                 option = list(initComplete =  JS("function(settings, json) {",
                                                                                  "$('body').css({'font-family': 'Helvetica'});",
                                                                                  "}")),
                                                 selection = list(mode = 'single', target = 'column'),
                                                 callback=JS(
                                                   'table.on("click.dt","td", function() {
                                                    var data=table.cell(this).data();
                                                    if (table.cell(this).data() < 0)
                                                    swal({title: "pegRNA-Score", text: "The higher the better! (calculated based on recommendation from the Liu Lab)",imageUrl: "Score.jpg",imageSize: "460x200"
                                                    });   
                                                    if (table.cell(this).data() > 0)
                                                    swal({title: "Edit position", text:  "Edit position is defined relative to the PAM sequence (Anzalone et.al. 2019)",imageUrl: "Target.jpg",imageSize: "460x200"
                                                    });                                                    
 
                                                    })'
                                                   
                                                 )
          
                                                 )
          
          output$myText <- renderText({
            
            guidesdf <- data.frame("Variant" = inFile[as.numeric(guides[,1]),]$Variant, "Protospacer(Sense)" = paste("cacc", guides[,2],"gtttt", sep = ""),"Protospacer(Antisense)" = paste("ctctaaaac",reverseComplement(DNAStringSet(guides[,2])),sep = "") , "EditPos." = guides[,3], "Extension(Sense)" = paste("gtgc",guides[,4],sep = ""), "Extension(Antisense)" = paste("aaaa",reverseComplement(DNAStringSet(guides[,4])), sep = ""), "PAM" = guides[,5], "PAM-Strand" = guides[,6], "Score" = as.numeric(guides[,7]), stringsAsFactors = FALSE)
            df.agg <- aggregate(Score ~ Variant, guidesdf, max)
            guidesdf$id  <- 1:nrow(guidesdf)
            df.max <- join(df.agg, guidesdf)
            df.max <- df.max[!duplicated(df.max[,"Variant"]),]
            df.max[order(df.max$id), ]
            
            x <- nrow(subset(df.max, as.numeric(df.max$Score) > -999))
            y <- nrow(inFile)
            
            
            paste(x, "out of", y, "Variants could be targeted")
            
            
          })
          
          
        }
        
        
            
      }) 
        
      
        
      }

         updateButton(session, "search",disabled = TRUE)

         
         shinyjs::runjs("window.scrollTo({
                        top: 0,
                        right: 200,
                        behavior: 'smooth'
                        })")
         
         shinyjs::showElement(id = "download")
         output$download <- renderUI({
           tagList(
           downloadButton('downloadData', "Download Results"),
           
           )
           })
       
      
    
    })
  

 
})
  