#Author: K.T.Shreya Parthasarathi
#Date: 12/06/2024
#Purpose: 5 year survival analysis of patients with breast cancer using mRNA expression profiles

setwd('path\\to\\transcriptomic\\data_file')

#Importing libraries
library(RTCGA.clinical)
library(survival)
library(dplyr)
library(survminer)
library(ggplot2)

#Data import and processing
dim(BRCA.clinical)

rna <- read.table('emt_brca.tsv', header=T,sep='\t', fill = TRUE)
dim(rna)
rna = distinct(rna)
dim(rna)
colnames(rna)[1] = 'bcr_patient_barcode' #Rename column name with patient ids as the one similar in clinical data file
head(rna)

clin = survivalTCGA(BRCA.clinical)
head(clin)
table(clin$patient.vital_status)

#Combining gene expression and clinical data
clin = rna %>%
  as_tibble() %>%
  select(bcr_patient_barcode,gene1,gene2,gene3) %>%
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>%
  inner_join(clin, by="bcr_patient_barcode")
names(clin)

dim(clin)
table(clin$patient.vital_status)

#5yr criteria:
all_clin = within(clin, patient.vital_status[patient.vital_status == 1 & times > 1825] <- 0)
all_clin = within(clin, times[times > 1825] <- 1825)
table(all_clin$patient.vital_status)

#Grouping into high and low expression based on median value of expession
gene = cut(all_clin$gene1, breaks = c(0, median(all_clin$gene1), Inf), labels = c("low", "high"))

#Survival anlysis with specific gene
sfit = survfit(Surv(times, patient.vital_status)~gene, data = all_clin)

#Cox regression
fit = coxph(Surv(times, patient.vital_status)~gene1, data = all_clin)

#Plotting
ggsurv = ggsurvplot(sfit, legend.labs=c("Low","High"), legend.title=
                      'Expression', title = 'gene1', pval = TRUE, pval.method = TRUE, xlab = "Time(in days)" )

ggsurv$plot + ggplot2::annotate("text",x = Inf, y = Inf, vjust = 1, hjust=1, 
                                label = "HR = x \n p(HR) = y")







