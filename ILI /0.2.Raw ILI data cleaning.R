setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
source('functions/Detrend and Smooth.R')
source('functions/series-util.r')
source('functions/matrix-util.r')
load('Rdata/raw clinical data.Rdata')

################################################################################

### Data Cleaning ###

###############################################################################

# load(file = "./data/raw data.RData")
## Aggregate the data from the clincis below ##

clinicGps <- list( Vic_135    = c(26, 52:56, 73),
                   Vic_079    = c(44 : 51),
                   VietMy     = c(43, 57, 77),
                   HoangKhang = c(33, 71),
                   ThanhCong  = c(35, 85),
                   CongHoa_1  = c(31),
                   CongHoa_2  = c(36),
                   CongHoa_4  = c(34),
                   GiaDinh    = c(88),
                   PhoiViet   = c(91) )


removedcols <- NULL
for (i in 1 : length(clinicGps)) {
        clinicName <- names(clinicGps)[i]
        
        clinicTotalP <- rep( NA, times = nrow(PATIENT) )
        clinicIliNum <- rep( NA, times = nrow(ILINUM) )
        
        for (gp in clinicGps[[i]]) {
                if ( gp %in% colnames(ILINUM) ) {
                        gpcol <- which(colnames(ILINUM) == gp)
                        removedcols <- append(removedcols, gpcol)
                        
                        clinicTotalP <- vectorSum(clinicTotalP, PATIENT[,gpcol])
                        clinicIliNum <- vectorSum(clinicIliNum, ILINUM[,gpcol])
                }
        }
        
        PATIENT <- cbind(PATIENT, clinicTotalP)
        ILINUM <- cbind(ILINUM, clinicIliNum)
        
        colnames(PATIENT)[ncol(PATIENT)] <- clinicName
        colnames(ILINUM)[ncol(ILINUM)] <- clinicName
}
PATIENT <- PATIENT[,-removedcols ]
ILINUM <- ILINUM[,-removedcols ]
dim(PATIENT)

rm(gp, gpcol, i, removedcols, clinicTotalP, clinicIliNum, clinicName, clinicGps)

## Calculate the ILI Percentage

PERCENT = ILINUM/PATIENT

## Remove the clinics on the blacklist ###

blackList <- as.character( c(30, 83, 93 : 104) )
ILINUM = ILINUM[,-which(colnames(ILINUM) %in% blackList)]
PATIENT = PATIENT[,-which(colnames(PATIENT) %in% blackList)]
PERCENT = PERCENT[,-which(colnames(PERCENT) %in% blackList)]

## Only clinics with more than 300 reports are selected to do analysis ##
selected = which(colSums(!is.na(ILINUM))>= 300)

ILINUM = ILINUM[,selected]
PATIENT = PATIENT[,selected]
PERCENT = PERCENT[,selected]
dim(PERCENT)

rm(blackList,selected)

## Remove the clinics that report more than 50 % 0 ilinum  ##
report_threshold = 0.5
total_report = colSums(!is.na(ILINUM))
zero_report = colSums(ILINUM == 0,na.rm = T)
zero_report_percent = zero_report/total_report

removedclinics = which(zero_report_percent > report_threshold)

ILINUM = ILINUM[,-removedclinics]
PATIENT = PATIENT[,-removedclinics]
PERCENT = PERCENT[,-removedclinics]
dim(PERCENT)

clinic_id = colnames(ILINUM)

rm(report_threshold,total_report,zero_report,zero_report_percent,removedclinics)

###############################################################################

##### Data Description ######

###############################################################################

#### How many messeages come from 33 clinics? ###
A = read.csv("datasets/ILI reports through 191231.csv")
sum(A$gpId %in% clinic_id)

### How many visits come from 33 clinics? ###
selected_clinics = A[A$gpId %in% clinic_id,]

sum(selected_clinics$totalPatients)
sum(selected_clinics$iliNum)

summary_num = summary(ILINUM)
summary_percent = summary(PERCENT)

m_patients_daily = rowMeans(PATIENT,na.rm = T)
summary(m_patients_daily)

m_ilinum_daily = rowMeans(ILINUM,na.rm = T)
summary(m_ilinum_daily)

m_iliperc_daily = rowMeans(PERCENT,na.rm = T)
summary(m_iliperc_daily)

boxplot(ILINUM,main = "The number of patients from each clinic")
boxplot(PERCENT,main = "The fraction of ILI patients from each clinic")

rm(clinic_id)
rm(list=lsf.str()) 

save(ILINUM,file = 'Rdata/clean ILINUM.RData')
save(PERCENT,file = 'Rdata/clean ILIPERCENT.RData')
save(PATIENT,file = 'Rdata/clean PATIENTS.RData')
