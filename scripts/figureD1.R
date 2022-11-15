############################
### R code for Figure D1 ###
############################

 rm(list=ls())

## The working directory needs to be set up below ##
 setwd("TO BE ADDED")

 data <- readxl::read_excel("results/figureD1input.xlsx",
                            col_names = FALSE)
 
  data <- as.data.frame(data)
     p <- matrix(seq(0.4,0.6,0.01),nrow = 1)
 unity <- data[-1,1]
   OVR <- data[-1,2]

  apr_unity <- data[1,1]
  apr_OVR <- data[1,2]
  
  apr_unity <- rep(apr_unity,length(p))
  apr_OVR <- rep(apr_OVR,length(p))    
  
pdf("results/figureD1.pdf")
par(mfrow=c(1,2))   

plot(p, unity, type = "l", lty = "solid", col = "blue", 
     xlab = "v", ylab = "marginal persuasion rate", 
     main = "Not Vote for Unity")   
lines(p,apr_unity, type = "l", lty = "dashed", col = "black")
legend("topright", c("marginal", "average"), 
       lty =c("solid", "dashed"), col = c("blue", "black"))

plot(p, OVR, type = "l", lty = "solid", col = "blue", 
     xlab = "v", ylab = "marginal persuasion rate", 
     main = "Vote for OVR")   
lines(p,apr_OVR, type = "l", lty = "dashed", col = "black")
legend("topleft", c("marginal", "average"), 
       lty =c("solid", "dashed"), col = c("blue", "black"))
dev.off()

