# Load packages ----
conf=read.table("EigenGWAS.conf", as.is = T)

library(shiny)
library(bsplus)

# Source helpers ----
source('helper.R')

unzip = "unzip -o -q plink.zip"
plink2 = "./plink_mac --allow-extra-chr"
if(length(grep("linux",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("linux")
  zip::unzip("plink.zip","plink_linux",overwrite = T)
  plink2 = "./plink_linux --allow-extra-chr"
} else if(length(grep("apple",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("apple")
  zip::unzip("plink.zip","plink_mac",overwrite = T)
  plink2 = "./plink_mac --allow-extra-chr"
  #  system("git rev-list head --max-count 1 > gitTag.txt")
} else {
  print("windows")
  zip::unzip("plink.zip","plink_win.exe",overwrite = T)
  plink2 = "plink_win.exe --allow-extra-chr"
}

options(shiny.maxRequestSize=conf[1,2]*1024^2, shiny.launch.browser=T)
gTag=read.table("gitTag.txt")

# Define UI for EigenGWAS Application
ui <- fluidPage(
  theme = "style.css",
  div(style = "padding: 1px 0px; width: '100%'",
      titlePanel(
        title = "",
        windowTitle = "EigenGWAS"
      )
  ),
  navbarPage(
    title = div(
      span(
        #img(src = "logo.png"),
        #strong("EigenGWAS"),
        HTML("<input type=button style='font-size:30px;border:0;height:35px' value=EigenGWAS onclick=\"window.history.go(-1)\">"),
        style = "position: relative; top: 30%; transform: translateY(-50%);"
      )
    ),
    id = "inNavbar",
    tabPanel(
      title = "Data Input",
      value = "datainput",
      fluidRow(
        column(
          6,
          fileInput(
            "file_input",
            paste0('Source files (.bim, .bed, .fam) [< ', conf[1,2],' MB]'),
            multiple = TRUE,
            accept = c("bed", "fam", "bim")
          ) %>%
            shinyInput_label_embed(
              icon("question-circle") %>%
                bs_embed_popover(
                  title = "The chromosome index (the first column of the .bim file) must be numeric", content = "", placement = "right"
                )
            )
        ),
        column(
          6,
          numericInput('autosome', "Autosome number", value=22)
        )
      ),
      hr(),
      fluidRow(
        column(
          2,
          radioButtons(
            'bred',
            'Population type',
            choices = list('Outbred' = 'outbred', 'Inbred' = 'inbred'),
            selected = 'outbred',
            inline = T
          ) %>%
            shinyInput_label_embed(
              icon("question-circle") %>%
                bs_embed_popover(
                  title = "INBRED is chosen if your sample has homogenous genome, otherwise choose OUTBRED", content = "", placement = "right"
                )
            )
        ),
        column(
          3,
          sliderInput(
            'maf_cut',
            'MAF threshold',
            value = 0.01, min = 0.01, max = 0.1, step = 0.01
          ) %>%
            shinyInput_label_embed(
              icon("question-circle") %>%
                bs_embed_popover(
                  title = "Marker with allele frequency lower than the given MAF threshold will be filtered out", content = "", placement = "right"
                )
            )
        ),
        column(
          3,
          sliderInput(
            'espace',
            'Eigen space',
            value = 2, min = 1, max = 5, step = 1
          ) %>%
            shinyInput_label_embed(
              icon("question-circle") %>%
                bs_embed_popover(
                  title = "This specifies the first few eigenvectors that will be scanned to obtain the loci under selection", content = "", placement = "right"
                )
            )
        ),
        column(
          3,
          sliderInput(
            'proportion',
            'Marker sampling proportion',
            value = 1, min = 0.2, max = 1,step=0.2
          ) %>%
            shinyInput_label_embed(
              icon("question-circle") %>%
                bs_embed_popover(
                  title = "A proprotion of the whole genome markers are sampled for quick GRM", content = "", placement = "right"
                )
            )
        )
      ),
      fluidRow(
        column(
          4,
          actionButton(
            'run',
            'EigenGWAS, Go!'
          )
        )
      )
    ),
    tabPanel(
      title = "Visualization",
      value = "visualization",
      fluidRow(
        column(12, 
               mainPanel(
                 tabsetPanel(
                   id = 'EgFunctions',
                   type = 'tabs',
                   tabPanel(
                     'MAF',
                     plotOutput('freq'),
                     htmlOutput('maf_note')
                   ),
                   tabPanel(
                     'GRM',
                     plotOutput('grm'),
                     htmlOutput('grs_note')
                   ),
                   tabPanel(
                     "Eigenanalysis",
                     mainPanel(
                       plotOutput("PCA")
                     ),
                     sidebarPanel(
                       sliderInput("x", "xLab", 1, 8, 1, step = 1),
                       sliderInput("y", "yLab", 1, 8, 2, step = 1)
                     )
                   ),          					
                   tabPanel(
                     'Eigenvalue',
                     plotOutput('eigenvalue'),
                     htmlOutput('eigenvalue_note')
                   ),
                   tabPanel(
                     'EigenGWAS',
                     tabPanel(
                       'EigenGWAS visualization',
                       mainPanel(plotOutput('eigengwas'),width = 6),
                       mainPanel(plotOutput('qqplot'),width = 6),
                       mainPanel(textOutput('manhattan_note'),width = 6), # Do not use htmlOutput(), it would make renderImage() run twice for no reason.
                       mainPanel(textOutput('qq_note'),width = 6),
                       sidebarPanel(
                         sliderInput(
                           'EigenGWASPlot_espace',
                           "Eigen Space", 
                           min = 1, max = 5, value = 1, dragRange=TRUE, step = 1),width = 4
                       )
                     )
                   ),
                   tabPanel(
                     "Top Hits",
                     sidebarPanel(
                       sliderInput(
                         'EigenGWASDT_espace',
                         "Eigen Space", 
                         min = 1, max = 5, value = 1, dragRange=TRUE, step = 1),width = 4
                     ),
                     mainPanel(dataTableOutput('topsnp'),width = 8)
                   )
                 )
               )
        )
      ),
      hr(),
      fluidRow(
        column(
          12,
          mainPanel(
            column(3,
                   downloadButton(
                     'eReport', 
                     'Summary Report'
                   )
            ),
            column(3,
                   downloadButton(
                     'fReport', 
                     'Full Report'
                   )
            )
          )
        )
      )
    ),
    # tabPanel(
    #   title = "Power",
    #   value = "power",
    #   fluidPage(
    #     # Application title
    #     titlePanel("EigenGWAS Power Calculator"),
    #     hr(),
    #     # Sidebar with a slider input for number of bins
    #     sidebarLayout(
    #       sidebarPanel(
    #         #        numericInput("n", "Sample size",
    #         #                     value = 500, min=100),
    #         numericInput("m", "Marker number",
    #                      value = 10000, min=1),
    #         sliderInput("w1",
    #                     "Subpop 1 proportion",
    #                     min = 0.05, max = 0.95, value = 0.5),
    #         sliderInput("p1", "Frequency at population 1",
    #                     min = 0.01, max=0.99, value=0.3),
    #         sliderInput("p2", "Frequency at population 2",
    #                     min = 0.01, max=0.99, value=0.5),
    #         selectInput("alpha", "Alpha",
    #                     choices = c("0.001", "0.005", "0.01", "0.05", "0.1"), 
    #                     selected="0.05"),
    #         selectInput("beta", "Beta",
    #                     choices = c("0.5", "0.6", "0.7", "0.8", "0.9"),
    #                     selected="0.8"),
    #         actionButton(
    #           'powerupdate',
    #           'Update power',
    #           icon = icon("refresh")
    #         )
    #         #submitButton(text="Update power", icon("refresh"))
    #       ),
    #       
    #       # Show a plot of the generated distribution
    #       mainPanel(
    #         plotOutput("distPlot")
    #       )
    #     )
    #   )
    # ),
    tabPanel(
      title = "About",
      value = "about",
      #tags$br(),
      tags$h3("Source Code"),
      tags$p(HTML("For EigenGWAS core algorithm implementation and R Shiny code in this web tool please refer to")),
      tags$p(HTML("<a href=\"https://github.com/GuoanQi1996/EigenGWASRshiny\" target=\"_blank\">GitHub repository: EigenGWASRShiny.</a>")),
      tags$p(HTML("For EigenGWAS implementation in JAVA, GEAR (GEnetic Analysis Repository), please refer to")),
      tags$p(HTML("<a href=\"https://github.com/gc5k/GEAR\" target=\"_blank\">GitHub repository: GEAR.</a>")),
      tags$br(),
      tags$h3("Citation"),
      tags$p(HTML("<a href=\"https://www.nature.com/articles/hdy201625\" target=\"_blank\">Chen, G.B. et al, EigenGWAS: finding loci under selection through genome-wide association studies of eigenvectors in structured populations, Heredity, 2016, 117:51-61.</a>")),
      tags$br(),
      tags$p(HTML("<a href=\"https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13370\" target=\"_blank\">Guo-An Qi et al, EigenGWAS: An online visualizing and interactive application for detecting genomic signatures of natural selection, <i>Molecular Ecology Resource</i>, 2021, <i>00</i>:1-13.</a>")),      tags$br(),
      tags$p(HTML(paste("Git version:", gTag[1,1])))
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  #Plot on the web
  currentFile <- reactive({
    withProgress(message="EigenGWAS:", value=0, {
      incProgress(1/3, detail = paste0(" check filesets ..."))
      FileLoad=0
      str=""
      if(length(which(grepl("*.bed", input$file_input$name)))  != 1) {
        str=paste(str, "No bed file found.") 
      } else {
        FileLoad=FileLoad+1
      }
      
      if(length(which(grepl("*.bim", input$file_input$name)))  != 1) {
        str=paste(str, "\nNo bim file found.")
      } else {
        FileLoad=FileLoad+1
      }
      
      if(length(which(grepl("*.fam", input$file_input$name)))  != 1) {
        str=paste(str, "\nNo fam file found.")
      } else {
        FileLoad=FileLoad+1
      }
      
      if (FileLoad < 3) {
        showNotification(str, duration = 5, type = "error")
        return()
      } else if (FileLoad > 3) {
        showNotification("More than 3 files selected", duration = 5, type="error")
      }
      
      idx=grep(".bed$", input$file_input$datapath)
      if (length(idx)==1) {
        rt=substr(input$file_input$datapath[idx], 1, nchar(input$file_input$datapath[idx])-4)
      }
      for (i in 1:3) {
        if (i != idx) {
          f1 = input$file_input$datapath[i]
          tl = substr(f1, nchar(f1)-2, nchar(f1))
          file.symlink(f1, paste0(rt, ".", tl))
        }
      }
      
      incProgress(1/3, detail = paste0(" check chromosome ..."))
      froot = substr(input$file_input$datapath[idx], 1, nchar(input$file_input$datapath[idx])-4)
      get_chr = read.table(paste0(froot,'.bim'),header=F,colClasses = c("character","NULL","NULL","NULL","NULL","NULL"))
      if (length(which(is.na(as.numeric(get_chr[,1]))))>0){
        showNotification("The chromosome index in the .bim file must be numeric!", duration = 5, type="error")
        stop("The chromosome index in the .bim file must be numeric! Refresh to continue.")
      }
      
      incProgress(1/3, detail = paste0(" data pre-filter ..."))
      if(as.numeric(input$maf_cut) == 0) {
        frootMAF=froot
      } else {
        frootMAF=paste0(froot, "_t") 
        frqMAF=paste(plink2, "--allow-no-sex --bfile", froot, "--autosome-num", input$autosome, "--maf", input$maf_cut, "--make-bed --out", frootMAF)
        system(frqMAF)
      }

      return (frootMAF)
    })
  })
  
  mark <- gsub('[-: ]','',as.character(Sys.time()))
  
  observeEvent(input$run, {
    updateNavbarPage(session, "inNavbar", selected = "visualization")
    updateTabsetPanel(session, "EgFunctions","EigenGWAS")
    
    froot <<- currentFile()
    withProgress(message="EigenGWAS:", value=0, {
      incProgress(1/1, detail = paste0(" collecting information ..."))
      PC = input$espace
      
      nn<-nrow(read.table(paste0(froot, ".fam"), as.is = T, header = F, colClasses = c("character","NULL","NULL","NULL","NULL","NULL")))
      mm<-nrow(read.table(paste0(froot,'.bim'), as.is = T, header=F, colClasses = c("character","NULL","NULL","NULL","NULL","NULL")))
    })
    
    #plink
    withProgress(message="EigenGWAS:", value=0, {
      time1 = proc.time()
      n = 4+PC
      pcRun=10
      incProgress(1/n, detail = paste0(" estimating freq ..."))
      frqCmd=paste(plink2, "--bfile", froot, "--autosome-num", input$autosome, "--freq --out", froot)
      system(frqCmd)
      
      incProgress(2/n, detail = paste0(" making grm and conduct PCA ..."))
      if (as.numeric(input$proportion)==1){
        grmCmd=paste(plink2, "--allow-no-sex --bfile", froot, "--autosome-num", input$autosome, "--make-grm-gz --pca", pcRun, "--out", froot)
      } else {
        grmCmd=paste(plink2, "--allow-no-sex --bfile", froot, "--autosome-num", input$autosome, "--make-grm-gz --pca", pcRun, "--out", froot,"--thin",as.numeric(input$proportion))
      }
      system(grmCmd)

      #EigenGWAS
      for(i in 1:PC) {
        incProgress(1/n, detail = paste0(" scanning eigen space ", i))
        outRoot=paste0(froot, ".", i)
        liCmd=paste0(plink2, " --allow-no-sex --linear --bfile ", froot, " --autosome-num ", input$autosome, " --pheno ", froot, ".eigenvec --mpheno ", i," --out ", outRoot)
        system(liCmd)
      }
      incProgress(1/n, detail = paste0(" finishing EigenGWAS."))
      time2 = proc.time()
      time = (time2-time1)[3][[1]]

      print(paste0(froot,' takes ',time,' seconds to finish the EigenGWAS analysis.'))
    })
    
    sc=ifelse(input$bred == 'inbred', 2, 1)
  
    withProgress(message="EigenGWAS complete! Visualizing:", value=0, {
      PC = input$espace
      n=4+PC
      
      incProgress(1/n, detail = paste0(" collecting MAF info  ... "))
      fq=read.table(paste0(froot, ".frq"), as.is = T, header = T)
      
      incProgress(1/n, detail = paste0(" collecting PCA info ... "))
      evalF=read.table(paste0(froot, ".eigenval"), as.is = T)
      pcF=read.table(paste0(froot, ".eigenvec"), as.is = T)
      evalF = as.numeric(evalF[,1])
      names(evalF) = c(1:pcRun)
      
      incProgress(1/n, detail = paste0(" collecting GRM info ... "))
      gz=gzfile(paste0(froot, ".grm.gz"))
      grm=read.table(gz, as.is = T)
      diagnol = grm[grm[,1]==grm[,2],4]
      off_diagnol = grm[grm[,1]!=grm[,2], 4]
      Ne=-1/mean(off_diagnol/sc)
      Me=1/var(off_diagnol/sc)
      
      incProgress(PC/n, detail = paste0(" collecting Eigenscan info ... "))
      GC = array(0, dim=PC)
      tophit = data.frame()
      for (i in 1:PC){
        eg = read.table(paste0(froot, ".", i, ".assoc.linear"), as.is = T, header = T, colClasses = c("numeric","character","numeric","NULL","NULL","NULL","NULL","numeric","numeric"))
        GC[i] = qchisq(median(eg$P, na.rm = T), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
        
        eg = eg[which(!is.na(eg$P)),]
        eg$Praw = eg$P
        gc = GC[i]
        eg$P = pchisq((eg$STAT)^2/gc, 1, lower.tail = F)
        eg = eg[which(!is.na(eg$P)),]
        assign(paste0("EgResDT",i),eg)

        eg_sig = eg[eg$P<=0.05,]
        eg_no_sig = eg[eg$P>0.05,]
        eg_com = rbind(eg_sig,eg_no_sig[sample(1:nrow(eg_no_sig),0.6*nrow(eg_no_sig)),])
        assign(paste0("EgResPlot",i),eg_com)
        
        eg$Espace = rep(i,dim(eg)[1])
        tophit = rbind(tophit, eg[order(eg$P),][c(1:10),c("Espace","CHR","SNP","BP","P","Praw")])
      }
      tophit[,c(5,6)] = format(tophit[,c(5,6)],digits = 4)
      
      egc=matrix(c(evalF[1:PC]/sc, GC), PC, 2, byrow = F)
      rownames(egc)=seq(1, PC)
      incProgress(1/n, detail = paste0(" Analysis complete. "))
      
      output$freq <- renderPlot({
        hist(fq$MAF, main="Minor allele frequency", xlab="MAF", xlim=c(0, 0.5), breaks = 50)
      })
      
      output$maf_note <- renderText({
        paste0("Minor allele frequencies for the ", mm," markers included for analysis.")
      })
      
      output$PCA <- renderPlot({
        layout(matrix(1:2, 1, 2))
        barplot(evalF/sc, border = F, main="Eigenvalue",ylim = c(0,max(evalF/sc)*1.2),xlab = 'Eigenspace')
        abline(h=1, lty=2, col="black")
        plot(main=paste0("Eigen space ", input$x, " vs ", input$y), pcF[,input$x+2], pcF[,input$y+2], xlab=paste0("Eigen space ", input$x), ylab=paste0("Eigen space ", input$y), bty='n', pch=16, cex=0.5,
             col=ifelse(pcF[,3]<0, "red", "blue"))
        abline(v=0, h=0, col="grey", lty=2)
      })
      
      output$grm <- renderPlot({
        layout(matrix(1:2, 1, 2))
        hist(off_diagnol/sc, main="Pairwise relatedness", xlab="Relatedness score", breaks = 50,freq = F)
        x=seq(range(off_diagnol/sc)[1],range(off_diagnol/sc)[2],length.out = 100)
        y=dnorm(x,-1/nn,sqrt(1/Me))
        lines(x,y,col="red",lwd=2)
        
        Ne=format(Ne, digits=3, nsmall=2)
        Me=format(Me, digits=3, nsmall=2)
        
        legen1 = bquote(italic(n)[e]==.(Ne)(.(nn)))
        legen2 = bquote(italic(m)[e]==.(Me)(.(mm)))
        legend("topright", legend = c(as.expression(legen1), as.expression(legen2)),bty ='n')
        
        hist(diagnol/sc, main="Self-relatedness", xlab="Relatedness score", breaks = 15)
      })
      
      output$grs_note <- renderText({
        "Relatedness score is defined as the pairwise relatedness for any pair of individuals as measured over genome-wide markers. It is often employed for the estimation for additive genetic variance, see VanRaden (<a href=\"https://www.sciencedirect.com/science/article/pii/S0022030208709901\" target=\"_blank\"><i>J Dairy Sci, 2008, 91:4414-4423</i></a>) for more details. <i>n</i><sub>e</sub> is the effective sample size. If the samples are related to each other much, <i>n</i><sub>e</sub> would be smaller than the real sample size (in brackets). <i>m</i><sub>e</sub> is the effective number of markers. When markers are in linkage equilibrium, <i>m</i><sub>e</sub>=<i>m</i>, the number of markers in study (in brackets). Of note, when the sample has experienced, recent, strong selection, <i>m</i><sub>e</sub> can be very small, say less than 0.01<i>m</i>; however, it can be of demographic factors possible."
      })
      
      output$eigenvalue <- renderPlot({
        barplot(t(egc), beside = T, border = F, xlab="Eigen Space", ylim=c(0,max(egc)+2))
        abline(h=1, lty=2, lwd=2)
        legend("topright", legend = c("Eigenvalue", expression(paste(lambda[gc]))), pch=15, col=c("black", "grey"), bty='n')
      })
      
      output$eigenvalue_note <- renderText({
        "Eigenvalue follows a mixture distribution <i>&eta;</i><i>F</i><sub>st,s</sub>+(1-<i>&eta;</i>)<i>F</i><sub>st,d</sub>, and <i>&lambda;</i><sub>GC</sub> is proportional to <i>F</i><sub>st,d</sub>. If Eigenvalue is far larger than <i>&lambda;</i><sub>GC</sub>, it indicates the presence of selection sweep of the sample."
      })

      output$eigengwas <- renderImage({
        width = session$clientData$output_eigengwas_width
        height = session$clientData$output_eigengwas_height
        pcIdx=input$EigenGWASPlot_espace[1]
        if (pcIdx > PC) {
          showNotification(paste0("Eigen space ", pcIdx, " hasn't been scanned yet."), type = "warning")
          return()
        }
        print(paste0("EG ePC ", pcIdx))
        
        EigenRes = get(paste0("EgResPlot",pcIdx))
        EigenRes$P[which(EigenRes$P==0)] = 1e-300
        man_cache = manhattanCache(pvalues = EigenRes$P,chr = EigenRes$CHR,pos = EigenRes$BP,ismlog10 = FALSE)
        png(filename = paste0(tempdir(),'/EgE',pcIdx,'.png'),width = width,height = height)
        manhattan(man_cache,significant = 0.05/nrow(eg),cex = 0.4,colorSet = c("grey","darkblue"),title=paste0("ePC",pcIdx))
        dev.off()
        
        outImage = paste0(tempdir(),'/EgE',pcIdx,'.png')
        list(src = outImage,
             width = width,
             height = height
        )
      }, deleteFile = FALSE)
      
      output$manhattan_note <- renderText({
        "Manhattan plot for EigenGWAS scanning on the chosen eigenvector. The p-value threshold, grey line, is set at genome-wide control for type I error rate of 0.05 (Bonferroni correction)."
      })
      
      output$qq_note <- renderText({
        "QQ plot for the test statistics for the chosen eigenvector. The grey points are for the test statistics without technical correction for genome-wide drift, and the dark ones with technical correction for genome-wide drift."
      })
      
      output$qqplot <- renderImage({
        width = session$clientData$output_qqplot_width
        height = session$clientData$output_qqplot_height
        pcIdx=input$EigenGWASPlot_espace[1]
        
        EigenRes = get(paste0("EgResDT",pcIdx))
        EigenRes$P[which(EigenRes$P==0)] = 1e-300
        EigenRes$Praw[which(EigenRes$Praw==0)] = 1e-300
        qq_cache_raw = qqPlotCache(EigenRes$Praw)
        qq_cache_gc = qqPlotCache(EigenRes$P)
        #QQplot
        png(filename = paste0(tempdir(),'/QQE',pcIdx,'.png'),width = width, height = height)
        qqPlotQ(qq_cache_raw,ci.level = NULL, cex = 0.4,title=paste0("ePC",pcIdx))
        qqPlotQ(qq_cache_gc,ci.level = NULL,newplot = F, cex = 0.4)
        dev.off()
        
        outImage = paste0(tempdir(),'/QQE',pcIdx,'.png')
        list(src = outImage,
             width = width,
             height = height
      )}, deleteFile = FALSE)
      
      output$topsnp <- renderDataTable({
        pcIdx = input$EigenGWASDT_espace[1]
        if (pcIdx > input$espace) {
          showNotification(paste0("eSpace ", pcIdx, " hasn't been scanned yet."), type = "warning")
          return()
        }
        
        outputDT = get(paste0("EgResDT",pcIdx))
        outputDT = outputDT[,c(1:4,6,5)]
        colnames(outputDT[c(5,6)]) = c("P","Pgc")
        outputDT[order(outputDT$P),][1:100,]
      },searchDelay = 1000)
      
      output$eReport <- downloadHandler(
        # For PDF output, change this to "report.pdf"
        # filename = function(){
        #   paste0('summary_report',mark,'.html')
        # },
        filename = 'summaryReport.html',
        content = function(file) {
          # summary_file = file.path(tempdir(),paste0('summary_report',mark,'.html'))
          # file.copy(summary_file,file)
          isRunable = TRUE
          if (!is.null(input$file_input$datapath)) {
            if ( !file.exists(paste0(froot, ".eigenval"))) {
              isRunable = FALSE
            } else if (!file.exists(paste0(froot, ".frq"))) {
              isRunable = FALSE
            } else if (!file.exists(paste0(froot, ".eigenvec"))) {
              isRunable = FALSE
            } else {
              for(i in 1:input$espace) {
                if (!file.exists(paste0(froot, ".", i, ".assoc.linear"))) {
                  isRunable = FALSE
                }
              }
            }
          } else {
            isRunable = FALSE
          }
          if(!isRunable) {
            showNotification("No results generated yet", duration = 5, type = "message")
            return()
          }
          
          tempReport <- file.path(tempdir(), "egReport.Rmd")
          file.copy("egReport.Rmd", tempReport, overwrite = TRUE)
          
          # Set up parameters to pass to Rmd document
          index=grep(".bed$", input$file_input$name)
          params <- list(froot = froot,
                         uploadfile = substr(input$file_input$name[index], 1, nchar(input$file_input$name[index])-4),
                         proportion = input$proportion,
                         espace = input$espace, 
                         sc = ifelse(input$bred == 'inbred', 2, 1),
                         pcut = 0.05,
                         nn = nn,
                         mm = mm,
                         ne = Ne,
                         me = Me,
                         GC = GC,
                         offDiag = off_diagnol,
                         Diag = diagnol,
                         hitDT = tophit,
                         width = session$clientData$output_eigengwas_width,
                         height = session$clientData$output_eigengwas_height)
          withProgress(message="Processing: ", value=0, {
            incProgress(1, detail = paste0("Please wait for a moment. This gonna takes some time in case of large marker number."))
            rmarkdown::render(tempReport, output_file = file,
                              params = params,
                              envir = new.env(parent = globalenv()))
          })
        }
      )
      
      output$fReport <- downloadHandler(
        filename = function(){
          paste0("FullReports.zip")
        },
        content = function(file) {
          pcIdx=input$EigenGWASPlot_espace[1]
          if (pcIdx > input$espace) {
            showNotification(paste0("eSpace ", pcIdx, " hasn't been scanned yet."), type = "warning")
            return()
          }
          
          withProgress(message="Processing: ", value=0, {
            incProgress(1, detail = paste0("Please wait for a moment. This gonna takes some time in case of large marker number."))
            froot = currentFile()
            files = NULL
            for (i in 1:input$espace){
              fname=paste0('FullReport.E',i,'.txt')
              EigenRes = get(paste0("EgResDT",i))
              write.table(EigenRes,fname,quote=F,col.names = T,row.names = F)
              if(i==1){
                zip::zip(file,fname)
              } else {
                zip::zip_append(file,fname)
              }
            }
          })
        })
    })
  })
  
 # observeEvent(input$powerupdate, {
  #  output$distPlot <- renderPlot({
  #    m=input$m
  #    alpha=as.numeric(input$alpha)
  #    pcut=alpha/m
  #    chiT=qchisq(pcut, 1, lower.tail = F)
      
  #    n=c(100,  200,  500,   1000,  1500,  2000,
  #        5000, 7500, 10000, 15000, 20000, 50000)
  #    PW=matrix(0, 2, length(n))
      
  #    w1=input$w1
  #    w2=1-w1
      
  #    p1=input$p1
  #   h1=2*p1*(1-p1)
  #   p2=input$p2
  #    h2=2*p2*(1-p2)
      
  #    p=w1*p1+w2*p2
  #    H=w1*h1+w2*h2
      
   #   for(i in 1:length(n)) {
  #      ncpA=4*n[i]*w1*w2*(p1-p2)^2/(2*p*(1-p))
  #      ncpD=n[i]*w1*w2*(h1-h2)^2/H
        
  #      PW[1,i]=pchisq(chiT, 1, ncp=ncpA, lower.tail = F)
  #      PW[2,i]=pchisq(chiT, 1, ncp=ncpD, lower.tail = F)
  #    }
  #    colnames(PW)=n
  #    par(las=2)
  #    barplot(PW, beside = T, border = F, ylab="Statistical power", xlab="Sample size")
  #    abline(h=as.numeric(input$beta), lty=2, col="grey")
  #    legend("topleft", legend=c("Add", "Dom"), pch=15, col=c("black", "grey"), bty='n')
  #  })
 # })
}

# Run the application
shinyApp(ui = ui, server = server)
