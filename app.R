
# library(shiny)
# data(airquality)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("CDF and pdf of the Gumbel copula"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "parameter",
                  label = "Parameter:",
                  min = -20,
                  max = 20,
                  value = -2,
                  step = .2),
      sliderInput(inputId = "nets",
                  label = "Net",
                  min = 1,
                  max = 200,
                  value = 20),
      sliderInput(inputId = "angel",
                  label = "Angel",
                  min = 1,
                  max = 180,
                  value = 30),
      sliderInput(inputId = "Theta",
                  label = "Theta",
                  min = 1,
                  max = 180,
                  value = 50)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  
  output$distPlot <- renderPlot({
    
    zz=function (x,y,d){
      x1=-log(x)
      x2=-log(y)
      xx1=x1^d
      xx2=x2^d
      Cxy=exp(-(xx1+xx2)^(1/d))
      cxy=Cxy*(1/x)*(1/y)*(xx1+xx2)^((2/d)-2)*(x1*x2)^(d-1)*(d-1+(xx1+xx2)^(1/d))
      return(list(cxy=cxy, Cxy=Cxy))
    }
    
    u1=runif(500)
    u2=runif(500)
    u=u1
    library(VineCopula)
 
    x <- seq(0.01, 1, length= input$nets)
    y<- seq(0.01, 1, length= input$nets)
    z =Z=Z.<- matrix(NA, input$nets,input$nets)
    mm<-function(x,y, d){
      n=length(x)
      b=x[1]
      i=1
      j=1
      #z <- matrix(NA, 10, 10)
      while(i<=n)
      {
        j=1
        while(j<=n){
          
          if ( d>0 )
          {
            u=u
            v=BiCopHinv2(u2, u1, family = 4, par = input$parameter)
            z[i,j]=zz(x[i],y[j],d)$cxy
            Z[i,j]=zz(x[i],y[j],d)$Cxy
          }
          if ( d<0 )
          {
            u=u
            v=1-BiCopHinv2(1-u2, u1, family = 4, par = -input$parameter)
            z[i,j]=zz(1-x[i],y[j],-d)$cxy
            Z[i,j]=y[j]-zz(1-x[i],y[j],-d)$Cxy
          }
          j=j+1
        }
        i=i+1
        
      }
      
      return(list(z, Z,u,v))
    }
    par(mfrow=c(1,3))
    delta=input$parameter #the parameter of Gumbel Copula'
    z=mm(x,y,delta)[[1]]
    Z=mm(x,y,delta)[[2]]
    persp(x, y,z , theta = input$Theta, phi = input$angel, expand = 0.7, col = "lightblue",
          main="pdf of Gumbel Copula")
    persp(x, y, Z, theta = input$Theta, phi = input$angel, expand = 0.7, col = "lightblue",
          main="CDF of Gumbel Copula")
    
    plot(mm(x,y,delta)[[3]],mm(x,y,delta)[[4]], col='red', pch=16)
    
    
    
  })
  
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
