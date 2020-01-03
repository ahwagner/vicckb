###################################################################################################
##################################### Plot disease enrichment #####################################
###################################################################################################

library(ggplot2)
library(plotrix)
library(viridis)
library(reshape2)
library(multtest)
library(stats)

setwd('~/Google Drive/MGI/VICC/Manuscript/Figures/misc_figures/')

###################################################################################################
#### Read in files
###################################################################################################
## Get disease count lists
data <- read.csv(file="~/Google Drive/MGI/VICC/Manuscript/Figures/misc_figures/Data/disease_counts.csv",header=T, stringsAsFactors = F)

###################################################################################################
#### summarize top level data for plotting
###################################################################################################
top_data <- aggregate(. ~ TopNode_disease + TopNode_doid, data=data[,3:ncol(data)], FUN=sum)

## Get a list of all of the database names
group_cols <- colnames(top_data)[!(colnames(top_data) %in% c("TopNode_disease","TopNode_doid"))]

## Get total counts for each cancer across databases
top_data$total <- apply(top_data[,colnames(top_data) %in% group_cols], 1, FUN = sum)
top_data$total_perc <- top_data$total/sum(top_data$total)*100

## Create top-level bins
top_data[(top_data$TopNode_disease %in% c("cancer")),c("TopNode_disease","TopNode_doid")] <- c("other cancers","other cancers")
top_data[(top_data$TopNode_disease %in% c("other")),c("TopNode_disease","TopNode_doid")] <- c("other disease","other disease")

## Figure out the max percent a disease is in any dataset
tmp_top <- top_data
tmp_top[,group_cols] <- apply(tmp_top[,group_cols], 2, function(x) x/sum(x)*100)
tmp_top$max <- apply(tmp_top[,group_cols], 1, FUN = max)
top_data <- merge(top_data, tmp_top[,c("TopNode_disease","TopNode_doid","max")])

# ## Get the top 10 cancers (includes top hit for all 6)
# top <- head(top_data[order(-top_data$total),], n = 5)
## Get diseases for plotting that make up more than 5% of the total
top <- top_data[which(top_data$max > 5),]

## Create a dataframe of only the top 10 + a binned "other cancer" category for all other cancers combined
top_data_scaled <- as.data.frame(top_data[order(-top_data$total),])
top_data_scaled[!(top_data_scaled$TopNode_disease %in% c(top$TopNode_disease, "benign neoplasm", "other disease")),c("TopNode_disease","TopNode_doid")] <- c("other cancers","other cancers")
top_data_scaled <- aggregate(. ~ TopNode_disease + TopNode_doid, data=top_data_scaled, FUN=sum)
top_data_scaled$TopNode_disease <- gsub(" ","\n",top_data_scaled$TopNode_disease)

###################################################################################################
#### Create pie charts of data
###################################################################################################
# pie(top_data_scaled[which(top_data_scaled$cgi > 0),"cgi"], labels=top_data_scaled[which(top_data_scaled$cgi > 0),"TopNode_disease"], main="CGI", col=viridis(length(top_data_scaled[,"cgi"]), option="B"))
# pie3D(top_data_scaled[which(top_data_scaled$cgi > 0),"cgi"], labels=top_data_scaled[which(top_data_scaled$cgi > 0),"TopNode_disease"], explode=0.1, main="CGI")

## Create pie charts for each database
create_a_pie <- function(x,y){
  p <- pie(x[which(x[,y] > 0),y], labels=x[which(x[,y] > 0),"TopNode_disease"], main=y, col=viridis(length(x[,y]), option="B"))
  return(p)
}
# png(file = paste(getwd(),"disease_by_database__piechart_cgi.png", sep = "/"), height=1200, width=1150, res=150)
#   print(create_a_pie(top_data_scaled,"cgi"))
# dev.off()
create_a_pie(top_data_scaled,"civic")

###################################################################################################
#### Create bar charts of data
###################################################################################################
## Melt the df
top_data_scaled_long <- melt(top_data_scaled, id.vars = c("TopNode_disease", "TopNode_doid"))

## Refactor to desired order
top_data_scaled_long$TopNode_disease <- factor(top_data_scaled_long$TopNode_disease, levels = c((unique(top_data_scaled[order(top_data_scaled$total,decreasing=F),"TopNode_disease"]))), exclude = NULL)

## Plot the disease by database
png(file = paste(getwd(),"disease_by_database.png", sep = "/"), height=1200, width=1150, res=150)
  ggplot(top_data_scaled_long[!(top_data_scaled_long$variable %in% c("total","total_perc","max")),], aes(x=TopNode_disease, y=value, fill = variable)) + geom_col(position = 'stack') + xlab("Disease") + ylab("Evidence") + guides(fill=guide_legend(title="Database")) + scale_fill_viridis(discrete = T, direction = -1) 
dev.off()


## Spaces are more desirable for the following plots
top_data_scaled_long$TopNode_disease <- gsub("\n"," ",top_data_scaled_long$TopNode_disease)
## Change the disease order
top_data_scaled_long$TopNode_disease <- factor(top_data_scaled_long$TopNode_disease, levels = as.vector(c("other disease", "benign neoplasm",as.character(unique(top_data_scaled_long[!(top_data_scaled_long$TopNode_disease %in% c("benign neoplasm","other disease")),"TopNode_disease"])))), exclude = NULL)
## Create a custom color palette
plma <- rev(viridis::plasma(n=nrow(top)))
new_palette = c("grey87","grey66", plma)

## Plot the databases by disease proportions
png(file = paste(getwd(),"disease_by_database_proportion.png", sep = "/"), height=1200, width=1150, res=150)
  ggplot(top_data_scaled_long[!(top_data_scaled_long$variable %in% c("total","total_perc")),], aes(x=variable, y=value, fill = TopNode_disease)) + geom_col(position = 'fill') + xlab("Database") + ylab("Proportion of Interpretations") + guides(fill=guide_legend(title="Disease")) + scale_fill_manual(values=new_palette) + theme_bw()
dev.off()

pdf(file = paste(getwd(),"disease_by_database_proportion.pdf", sep = "/"), height=9, width=8.625)
  ggplot(top_data_scaled_long[!(top_data_scaled_long$variable %in% c("total","total_perc")),], aes(x=variable, y=value, fill = TopNode_disease)) + geom_col(position = 'fill') + xlab("Database") + ylab("Proportion of Interpretations") + guides(fill=guide_legend(title="Disease")) + scale_fill_manual(values=new_palette) + theme_bw()
dev.off()

###################################################################################################
#### Look for disease enrichment per database
###################################################################################################
top_data_scaled$TopNode_disease <- gsub("\n"," ",top_data_scaled$TopNode_disease)

## Get a list of cancer types to test
cancer_to_test <- top_data_scaled$TopNode_disease[!(top_data_scaled$TopNode_disease %in% c("other disease","other cancers"))]

## Restrict to desired columns and eliminate generic other disease or cancer rows
disease_vs_dataset <- top_data_scaled[!(top_data_scaled$TopNode_disease %in% c("other disease","other cancers")),c("TopNode_disease",group_cols)]

## Convert to table
rownames(disease_vs_dataset) <- disease_vs_dataset$TopNode_disease
disease_vs_dataset$TopNode_disease <- NULL
chisquared_all <- chisq.test(disease_vs_dataset)

## Restrict to desired columns
all_disease_vs_dataset <- top_data_scaled[,c("TopNode_disease",group_cols)]

## Create a results table
results = matrix(nrow=length(cancer_to_test), ncol=3)
colnames(results) <- c("disease","pvalue", "evidence_count")

## Create contingency tables
for(c in 1:length(cancer_to_test)){
  ctab <- all_disease_vs_dataset
  ctab[!grepl(cancer_to_test[c], ctab$TopNode_disease),"TopNode_disease"] <- paste("not",cancer_to_test[c])
  ctab <- aggregate(. ~ TopNode_disease, data=ctab, FUN=sum)
  rownames(ctab) <- ctab$TopNode_disease
  ctab$TopNode_disease <- NULL
  chisq_results = chisq.test(ctab)
  results[c,"disease"] <- cancer_to_test[c]
  results[c,"pvalue"] <- chisq_results$p.value
  results[c,"evidence_count"] <- sum(ctab[cancer_to_test[c],])
  #write.table(ctab, file=paste0(cancer_to_test[c],"_contingency_table.tsv"), sep="\t", row.names = F, quote = F)
}

#Correct p-values
pvalues=as.numeric(results[,"pvalue"])
pvalues_adj=mt.rawp2adjp(pvalues, proc=c("Bonferroni","BH"))
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
results=cbind(results, pvalues_adj_orig_order[,2:3])
results <- as.data.frame(results)


## Create a results table
db_results = matrix(nrow=length(group_cols), ncol=3)
colnames(db_results) <- c("database","pvalue","evidence_count")

## Create contingency tables
for(d in 1:length(group_cols)){
  ctab <- all_disease_vs_dataset
  others <- group_cols[!(group_cols %in% c(group_cols[d]))]
  ctab$others <- rowSums(ctab[,others])
  ctab <- ctab[,c("TopNode_disease",group_cols[d],"others")]
  rownames(ctab) <- ctab$TopNode_disease
  ctab$TopNode_disease <- NULL
  chisq_db_results = chisq.test(ctab)
  db_results[d,"database"] <- group_cols[d]
  db_results[d,"pvalue"] <- chisq_db_results$p.value
  db_results[d,"evidence_count"] <- sum(ctab[,group_cols[d]])
  #write.table(ctab, file=paste0(group_cols[d],"_contingency_table.tsv"), sep="\t", row.names = F, quote = F)
}

#Correct p-values
pvalues=as.numeric(db_results[,"pvalue"])
pvalues_adj=mt.rawp2adjp(pvalues, proc=c("Bonferroni","BH"))
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
db_results=cbind(db_results, pvalues_adj_orig_order[,2:3])
db_results <- as.data.frame(db_results)


## Create a results table
results = matrix(nrow=length(cancer_to_test)*length(group_cols), ncol=5)
colnames(results) <- c("disease", "database", "pvalue", "interpretation_freq", "interpretation_freq_other_databases")
n=0

## Create contingency tables
for(c in 1:length(cancer_to_test)){
  for(d in 1:length(group_cols)){
    n=n+1
    ctab <- all_disease_vs_dataset
    ctab[!grepl(cancer_to_test[c], ctab$TopNode_disease),"TopNode_disease"] <- paste("not",cancer_to_test[c])
    ctab <- aggregate(. ~ TopNode_disease, data=ctab, FUN=sum)
    others <- group_cols[!(group_cols %in% c(group_cols[d]))]
    ctab$others <- rowSums(ctab[,others])
    ctab <- ctab[,c("TopNode_disease",group_cols[d],"others")]
    rownames(ctab) <- ctab$TopNode_disease
    ctab$TopNode_disease <- NULL
    chisq_results = chisq.test(ctab)
    results[n,"disease"] <- cancer_to_test[c]
    results[n,"database"] <- group_cols[d]
    results[n,"pvalue"] <- chisq_results$p.value
    results[n,"interpretation_freq"] <- ctab[cancer_to_test[c],group_cols[d]]/sum(ctab[,group_cols[d]])
    results[n,"interpretation_freq_other_databases"] <- ctab[cancer_to_test[c],"others"]/sum(ctab[,"others"])
    write.table(ctab, file=paste0("contingency_tables/", group_cols[d], "_", gsub(" ", "_", cancer_to_test[c]), "_contingency_table.tsv"), sep="\t", row.names = F, quote = F)
  }
}

#Correct p-values
pvalues=as.numeric(results[,"pvalue"])
pvalues_adj=mt.rawp2adjp(pvalues, proc=c("Bonferroni","BH"))
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
results=cbind(results, pvalues_adj_orig_order[,2:3])
results <- as.data.frame(results)
results[,3:ncol(results)] <- apply(results[,3:ncol(results)], 2, FUN=as.numeric)
results$representation <- NA
results[which(as.numeric(results$pvalue) >0.05),"representation"] <- "not significant"
results[which(as.numeric(results$pvalue) <0.05 & as.numeric(results$interpretation_freq) > as.numeric(results$interpretation_freq_other_databases)),"representation"] <- "overrepresented"
results[which(as.numeric(results$pvalue) <0.05 & as.numeric(results$interpretation_freq) < as.numeric(results$interpretation_freq_other_databases)),"representation"] <- "underrepresented"

write.table(results, file="database_vs_disease.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



library("gplots")
# 1. convert the data as a table
rownames(all_disease_vs_dataset) <- all_disease_vs_dataset$TopNode_disease
all_disease_vs_dataset$TopNode_disease <- NULL
dt <- as.table(as.matrix(all_disease_vs_dataset))
# 2. Graph
png(file = paste(getwd(),"disease_by_database_balloonplot.png", sep = "/"), height=1200, width=1550, res=150)
  balloonplot(t(dt), main = "Diseases by database",xlab ="", ylab="", label = FALSE, show.margins = FALSE)
dev.off()
png(file = paste(getwd(),"disease_by_database_balloonplot_noother.png", sep = "/"), height=1200, width=1550, res=150)
  balloonplot(t(dt[cancer_to_test[!(cancer_to_test %in% "benign neoplasm")],]), main = "Diseases by database", xlab ="", ylab="", label = FALSE, label.size=0.5)
dev.off()
