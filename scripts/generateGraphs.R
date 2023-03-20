#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

quantificationEditingFreq <- Sys.glob(file.path("./CRISPRessoAggregate_on_*","CRISPRessoAggregate_quantification_of_editing_frequency.txt"))
#quantificationEditingFreq is a tsv of with header
# Name	Unmodified%	Modified%	Reads_total	Reads_aligned	Unmodified	Modified	Discarded	Insertions	Deletions	Substitutions	Only Insertions	Only Deletions	Only Substitutions	Insertions and Deletions	Insertions and Substitutions	Deletions and Substitutions	Insertions Deletions and Substitutions

for (tableFile in quantificationEditingFreq){

QEFdf <- read.table(file = tableFile, sep = "\t", header = 1) %>% mutate(Samples = substring(Name, first = 36)) %>% mutate(Samples=factor(Samples,levels=str_sort(Samples, decreasing = TRUE, numeric=TRUE)))
#long format read-classes
QEFdf2 <- QEFdf %>% pivot_longer(c("Unmodified","Modified"),names_to="Reads",values_to="Count")
#plot as stacked barplot
saveName <- substring(tableFile, first = 26, last = 40) %>% paste("reads") %>% paste(".pdf")
modifiedBarPlot <- ggplot(QEFdf2, aes(x=Samples, y = Count, fill=Reads)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  coord_flip() + theme(axis.text=element_text(size=6)) +
  scale_fill_manual(values = c("orange", "darkgreen"))
ggsave(modifiedBarPlot, file = saveName)


#long format read-classes
QEFdf3 <- QEFdf %>% pivot_longer(c("Insertions","Deletions","Substitutions"),names_to="Reads",values_to="Count")
saveName <- substring(tableFile, first = 26, last = 40) %>% paste("variants") %>% paste(".pdf")

#plot as stacked barplot
variantsPlot <- ggplot(QEFdf3, aes(x=Samples, y = Count, fill=Reads)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  coord_flip() + theme(axis.text=element_text(size=6)) +
  scale_fill_manual(values = c("orange", "darkgreen","darkblue"))
ggsave(variantsPlot, file = saveName)

 #convert to frequency (%), add unmodified (6) and indelsub (9,10,11)
QEFdf4 <- QEFdf

for (row in 1:nrow(QEFdf)){
    Reads_total = QEFdf4[row,6] + QEFdf4[row,9] + QEFdf4[row,10] + QEFdf4[row,11]
    QEFdf4[row,6] = QEFdf4[row,6]/Reads_total
    QEFdf4[row,9] = QEFdf4[row,9]/Reads_total
    QEFdf4[row,10] = QEFdf4[row,10]/Reads_total
    QEFdf4[row,11] = QEFdf4[row,11]/Reads_total
}

QEFdf4 <- QEFdf4 %>% pivot_longer(c("Insertions","Deletions","Substitutions","Unmodified"),names_to="Reads",values_to="Frequency")
#plot as stacked barplot
saveName <- substring(tableFile, first = 26, last = 40) %>% paste("frequency") %>% paste(".pdf")
frequencyPlot <- ggplot(QEFdf4, aes(x=Samples, y = Frequency, fill=Reads)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  coord_flip() + theme(axis.text=element_text(size=6)) +
  scale_fill_manual(values = c("orange", "green3","dodgerblue","black"))
ggsave(frequencyPlot, file = saveName)

}


knockoutReport <- Sys.glob(file.path("KnockoutReport.tsv"))
#tsv with header Sample SampleName	Reference	Unmodified	Reads_aligned	Knockout	Reason	Insertions	Deletions	Substitutions
koDF = read.table(file = knockoutReport, sep = "\t", header = 1)

#extract csv values (representing alleles) to long format
decompose <- function(sample,reference,snvtype,indelsub,Reason){

    resultDF <- data.frame(name=character(0),allele = character(0),type=character(0),count=numeric(0),reason=character(0))
    counts <- as.list(strsplit(indelsub, split = ","))[[1]]

    for (i in counts){
        if (i == 0) {
            next
        }
        rowDF <- data.frame(name=sample,allele = reference,type=snvtype,count=as.integer(i),reason=Reason)
        resultDF <- rbind(resultDF, rowDF)
    }

    return(resultDF)
}

indelsubDF <- data.frame(name=character(0),allele = character(0),type=character(0),count=numeric(0),reason=character(0))
#Iterate over koDF and decompose indelsub counts to long format (tsv with csv cell values) # nolint

for (row in 1:nrow(koDF)){

    if( !(is.na(koDF[row, 8]))){
    ins <- decompose(koDF[row, 2], koDF[row, 3], "Insertions", koDF[row, 8], koDF[row, 7])
    indelsubDF <- rbind(indelsubDF, ins)
    }

    if( !(is.na(koDF[row, 8]))){
    dels <- decompose(koDF[row, 2], koDF[row, 3], "Deletions", koDF[row, 9], koDF[row, 7])
    indelsubDF <- rbind(indelsubDF, dels)
    }

    if( !(is.na(koDF[row, 8]))){
    subs <- decompose(koDF[row, 2], koDF[row, 3], "Substitutions", koDF[row, 10], koDF[row, 7])
    indelsubDF <- rbind(indelsubDF, subs)
    }
}
#save just incase
write.csv(indelsubDF,"InDelSub.csv",row.names=FALSE)

#dotplot of indelsubs by allele
saveName <- paste("indelsub_determination") %>% paste(".pdf")
indelsub_determination <- ggplot(indelsubDF, aes(x=reason, y = count, color=type)) +
  geom_jitter(height = 0) +
  scale_fill_manual(values = c("orange", "darkgreen","darkblue")) +
  ylab("Length of mutation") + xlab("Determination")
ggsave(indelsub_determination, file = saveName)

counter <- 1
colors <- c("orange", "darkgreen","darkblue")
for (alleleName in unique(indelsubDF$allele)){

splitDF <- indelsubDF[indelsubDF$allele == alleleName,]

#dotplot of indelsubs by allele
saveName <- paste(sprintf("%s_indelsub_allele", alleleName)) %>% paste(".pdf")
indelsub_allele <- ggplot(splitDF, aes(x=name, y = count, fill=type)) +
  geom_dotplot(binaxis = "y",stackdir = "center", dotsize = 0.3, stackgroups = TRUE) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ylab("Length of mutation") + ggtitle(alleleName) +
  scale_y_log10() +
  theme_bw()
ggsave(indelsub_allele, file = saveName)
counter <- counter + 1

}
