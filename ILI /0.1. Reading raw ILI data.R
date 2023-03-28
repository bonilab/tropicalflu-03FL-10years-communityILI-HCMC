
################################################################################

## Reading data ##

###############################################################################
rm(list = ls())
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
A = read.csv("datasets/ILI reports through 191231.csv")

source("functions/daynum.R")

A$date = as.POSIXct(A$date,tz = "UTC",format = "%m/%d/%Y")

A$iliPerc = A$iliNum/A$totalPatients

if(length(which(is.na(A$iliPerc))) > 1){
        
        A = A[-which(is.na(A$iliPerc)),]
        
}

## Specify the timeseries for analysis ###
startdate = as.POSIXct("2010-01-01",tz = "UTC", format = "%Y-%m-%d")

enddate = as.POSIXct("2019-12-31",tz = "UTC")

A = A[which(A$date >= startdate & A$date<= enddate),]

unique_gp_ids = unique( A$gpId )

timeseries = seq(from = startdate,to = enddate,by = "day")

# Build the matrix to contain the ilinum for each clinic on each day
ILINUM = matrix(NA,nrow = length(timeseries),0)

# Calculate ilipercent for each record
PATIENT = matrix(NA,nrow = length(timeseries),0)

clinic_id = c()

for( i in 1:length(unique_gp_ids) ) {
        
        # Extract only this GPs records
        B = A[A$gpId == unique_gp_ids[i],]
        
        # How many dates are there?
        date_index = c(1:dim(B)[1])
        
        # Exclude the clinic with no reports after 2010
        if (dim(B)[1] > 0) {
                
                # Assign the ID the each column
                clinic_id = c(clinic_id, unique_gp_ids[i])
                
                # temporary matrix to be combined to ILINUM
                TEMP = matrix(NA,nrow = length(timeseries),1)
                
                # temporary matrix of PERCENT
                PATIENTTEMP = matrix(NA,nrow = length(timeseries),1)
                
                for( j in date_index ) {
                        # Calculate the daynum for each date
                        daynumber = daynum(B$date[j], startdate) + 1
                        
                        # For this clinic, put the ilinum in corresponding daynum in the ILINUM matrix
                        TEMP[daynumber,1] = B$iliNum[j]
                        
                        # Calculate the ilipercent
                        PATIENTTEMP[daynumber,1] = B$totalPatients[j]
                        
                        #         ## Check if totalpatient has 0 cases, which is missing value ###
                        #         #if (B$totalPatients[j] == 0) {
                        #                 
                        #                 PERCENTTEMP[daynumber,1] = -99
                        #         }
                        #         
                        # }
                        
                }
                ILINUM = cbind(ILINUM, TEMP)
                PATIENT = cbind(PATIENT,PATIENTTEMP)
        }
}


colnames(ILINUM) = clinic_id
colnames(PATIENT) = clinic_id

rm(A,unique_gp_ids,B,daynumber,startdate,enddate,timeseries,TEMP,i,j,date_index,PATIENTTEMP)
rm(list=lsf.str()) 

save.image(file = "Rdata/raw clinical data.Rdata")

