#!/usr/bin/env Rscript

args<- commandArgs(trailingOnly=TRUE)

#read in the data
data<-read.csv('/home/benjamin/Desktop/PIGv_gvprodigal/viral_screening/viral_screening.csv')
data<- read.csv(args[1], header = TRUE) 
logit_model <- readRDS("./screening_model.rds")
logit_model<-readRDS('/home/benjamin/Desktop/PIGv_gvprodigal/viral_screening/screening_model.rds')
#rename queries alligned to blast
colnames(data)[2]<-'blast'
predicted_prob <- predict(logit_model, newdata = data, type = "response")

#append it to the dataframe
data$Chance_Of_Finding_Genome<-predicted_prob
data
colnames(data)[2]<-'queries_alligned'

data
write.csv(data, file = 'viral_screening.csv', row.names = FALSE)
