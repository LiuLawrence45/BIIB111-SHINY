rm(list=ls())
library(shiny)
library(survminer)
library(survival)
library(readxl)
library(plyr)
library(readxl)
library(drc)
library(dr4pl)
library(reshape2)
library(stringr)
library(shinythemes)
path = 'NSR-REP-02_Covance_ELISPOT_PROD_02JUL2020_v1.xlsx'

data = data.frame(read_excel(path))

# exclude cv.s.t<80%
data.0=data[data$Analyte.ID=='% Cell Viability',]
data.0=data.0[data.0$Lab.Result.Value>80,]
data.0=data.0[,c(2,13)]
colnames(data.0)=c('subject','time.point')

data$assay=substr(data$Analyte.ID,1,nchar(data$Analyte.ID)-6)
data=data[,c(2,13,17,19,23)]
colnames(data)=c('subject','time.point','analyte','value','assay')
data=data[substr(data$assay,1,3)!='CD3',]
data=data[(grep('Rep', data$analyte, fixed = T)),]
data=data[-grep('Human',data$analyte,fixed=T),]
data$value=as.numeric(data$value)
data=data[!is.na(data$value),]

data=merge(data,data.0,by=c('subject','time.point'))


####define cvs

n=length(unique(data$assay))
id=unique(data$assay)
id.color=c('Green','pink','chocolate1','blue','black')
cvs=data.frame()

for (i in 1:n){
  temp=data[data$assay==id[i],]
  temp$ID=paste0(temp$subject, temp$time.point,temp$assay)
  mean.temp = tapply(temp$value, temp$ID,mean,na.rm=T)
  sd.temp= tapply(temp$value, temp$ID,sd,na.rm=T)
  criteria.temp= sd.temp^2/mean.temp
  cv.temp=sd.temp/mean.temp*100
  data.temp=data.frame(mean=mean.temp,cv=cv.temp,criteria=criteria.temp)
  data.temp$color=id.color[i]
  cvs=rbind(data.temp,cvs)
}
cvs$id=row.names(cvs)



########UI
ui <- fluidPage(theme = shinytheme("slate"),
  
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("BIIB111 "),
      # tags$hr(),
      # # fileInput("file1", "Choose CSV File :",
      # #           accept = c(
      # #             "text/csv",
      # #             "text/comma-separated-values,text/plain",
      # #             ".csv")
      # # ),
      # #checkboxInput("header", "Header", TRUE),
      # tags$hr(),
      # 
      # 
      # submitButton("Submit", icon("refresh")),
      
      selectInput("sample", "Choose a sample:",
                  data$subject
      ),
      sliderInput(inputId = "bins",
                  label = "Number of bins:",
                  min = 1,
                  max = 50,
                  value = 30),
      
      
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  
                  tabPanel("Plot", plotOutput("Plot")
                  ),
                  tabPanel("Dataset",
                    # checkboxInput("checkbox", label = "Display Time Point", value = TRUE),
                    # checkboxInput("checkbox", label = "Display Analyte", value = TRUE),
                    # checkboxInput("checkbox", label = "Display Assay", value = TRUE),hr(),
                    # fluidRow(column(3, verbatimTextOutput("value"))),

                  tableOutput("contents"),
                          ),
                  tabPanel("Summary Data",
                           plotOutput("hist"),hr(),
                           plotOutput("cvMean"),
                           )
                  
      )
      
    )
    
  )
)




server <- function(input, output) {
  output$contents <- renderTable({
    data1 = data[data$subject==input$sample,]
    data1 <- as.data.frame(lapply(data1, as.character))
    data1
    
  
  })
  
  output$cvMean <- renderPlot({
    par(mar=c(5,5,5,5))

    k=14
    rgx=range(cvs$mean)
    rg1=rgx[1]-(rgx[2]-rgx[1])/k
    rg2=rgx[2]+(rgx[2]-rgx[1])/k
    rgy=range(cvs$cv)
    rg3=rgy[1]-(rgy[2]-rgy[1])/k
    rg4=rgy[2]+(rgy[2]-rgy[1])/k
    plot(cvs$mean, cvs$cv, xlim=c(rg1, rg2), ylim=c(rg3, rg4), xlab="Mean",
         ylab = 'CV', col=cvs$color, pch=19, cex=0.2)
    lines(smooth.spline(cvs$mean, cvs$cv,df=4), col='blue', lty=1, lwd=2)
    abline(h=50,lty=2)
    title("BIIB111")
    legend('topright', id, col=id.color, pch=19, cex=0.5)
    length(cvs[cvs$cv<50,]$cv)/length(cvs$cv)
    
    #input$cutoff
    cvs.cut=cvs[cvs$mean<500,]

    plot(cvs.cut$mean, cvs.cut$cv, xlab="Mean",xlim=c(-50,550),
         ylab = 'CV', col=cvs$color, pch=19, cex=0.2)
     lines(smooth.spline(cvs$mean, cvs$cv,df=4), col='blue', lty=1, lwd=2)
     abline(h=50,lty=2)
    # lines(smooth.spline(cvs$mean, cvs$cv,df=4), col='blue', lty=1, lwd=2)
    title("CV vs Mean")
    legend('topright', id, col=id.color, pch=19, cex=0.5)
  })
  
 
  output$hist <- renderPlot({

    data.neg.his=data[data$assay=='Neg Control',]
    x    <- data.neg.his$value
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(data.neg.his$value ,xlab='Spot Counts',ylab='Frequency',main='Histogram of NC Spot Counts', breaks=bins)
    box()
  })
  
  output$Plot <- renderPlot({

    
    temp=data[data$subject==input$sample,]
    temp=temp[!is.na(temp$value),]
    varPlot(value~time.point/assay,temp,Title=list(main=paste0(input$sample)),
            Boxplot=list(jitter=1),VarLab = list(las=2, col='black'),
            Points=list(pch=16, cex=0.5), YLabel = list(text='Count',cex=0.8)
    )

    
    
  })
  
  
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)










