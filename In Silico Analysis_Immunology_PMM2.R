
##################################################################################
############# Exploratory analysis of transcriptome profiles #####################
##################################################################################
wd <- "C:/Users/pmgra/OneDrive/Documentos/FCT/FCT/UCIBIO/Glycoimmunology/Fibroblasts/Transcriptomic Code"
setwd(wd)


load("Transcriptomics_data.RData")


# Creating a DGEList object for use in edgeR
library(edgeR)

y <- DGEList(counts=readCounts,group=group)

# Filtering low expressed genes

keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes=FALSE]

# Calculate normalizing factors

y <- calcNormFactors(y)

#Creating a CPMS matrix (normalized values for PCA and clustering)

logCPMs<- edgeR::cpm(y, offset = y$offset, log = TRUE)


#########
#Figure 2 C - Principal component analysis biplots of human skin fibroblasts before and upon TNF-α stimulus
########

library(PCAtools)

#PCA

metadata <- data.frame(group)
rownames(metadata) <- colnames(logCPMs)

#Run PCA

pca.res <- pca(logCPMs, metadata=metadata)

#Plot variance explained by each component

# Plot PC1 versus PC2

biplot(pca.res, colby="group",title = "Biplot PC1 vs PC2") # Biplot with colors by sample type
png("myplot.png", width=900, height=600)
biplot(pca.res, x="PC1", y="PC2",lab="",colby="group",legendPosition = 'right', labSize = 200,legendLabSize = 20,pointSize =6 )
dev.off()




#To Remove PMM2-CDG sample which is more similar to the WT transcriptomic profile

metadata <- data.frame(metadata, "sample" = rownames(metadata))
metadata[,"sample"]<- sub("CTL", "WT", metadata[,"sample"])
metadata[,"sample"] <- sub("S$", "Stim", metadata[,"sample"])
metadata <- metadata[!(metadata[,"sample"] %in% c("4PMM2Stim", "10PMM2NStim")),]

group <- factor(metadata[,"group"], levels=unique(metadata[,"group"]))
colnames(readCounts) <- sub("CTL", "WT", colnames(readCounts))
colnames(readCounts) <- sub("NS$", "NStim", colnames(readCounts))
colnames(readCounts) <- sub("S$", "Stim", colnames(readCounts))

readCounts <- readCounts[, !(colnames(readCounts) %in% c("4PMM2Stim", "10PMM2NStim"))]
colnames(readCounts) <- rownames(metadata)

pca.res <- pca(logCPMs, metadata = metadata)
biplot(pca.res, x="PC1", y="PC2",lab="",colby="group",legendPosition = 'right', labSize = 200,legendLabSize = 15,pointSize =6 ) # Biplot with PC1 and PC2



##########################################################################################
# Differential expression analysis  ####
##########################################################################################


library(edgeR); library(ggplot2)


load("ReadCounts1June2021.RData")
load("MetadataJune2021.RData")

readCounts <- ReadCounts1

#To Remove PMM2-CDG sample which is more similar to the WT transcriptomic profile
Metadata[,"sampletype"][c(1,2,3,7,8)] <- c("WTStim","WTStim","WTStim","WTNStim","WTNStim")
Metadata[,"sample"]<- sub("CTL", "WT", colnames(readCounts))
Metadata[,"sample"]<- sub("NS$", "NStim", Metadata[,"sample"])
Metadata[,"sample"] <- sub("S$", "Stim", Metadata[,"sample"])
Metadata <- Metadata[!(Metadata[,"sample"] %in% c("4PMM2Stim", "10PMM2NStim")),]

group <- factor(Metadata[,"sampletype"], levels=unique(Metadata[,"sampletype"]))


colnames(readCounts) <- sub("CTL", "WT", colnames(readCounts))
colnames(readCounts) <- sub("NS$", "NStim", colnames(readCounts))
colnames(readCounts) <- sub("S$", "Stim", colnames(readCounts))
readCounts <- readCounts[, !(colnames(readCounts) %in% c("4PMM2Stim", "10PMM2NStim"))]




y <-DGEList(counts=readCounts,group=group)


# Filtering low expressed genes
keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes=FALSE]

# Calculate normalizing factors
y <- calcNormFactors(y)

# Estimating the dispersion and Fit model
design_matrix <- model.matrix(~0+group); colnames(design_matrix) <- gsub("group", "", colnames(design_matrix))
contrasts <- makeContrasts("WTStim_vs_WTNStim"=WTStim-WTNStim,"PMM2Stim_vs_PMM2NStim"=PMM2Stim-PMM2NStim,levels=design_matrix)
y <- estimateDisp(y, design_matrix)
fit <- glmQLFit(y,design_matrix)


# Identify DEGs - CTLStim_vs_ControlNS
qlf_WT <- glmQLFTest(fit, contrast = contrasts[,1])
DEG_WT <-topTags(qlf_WT, n=nrow(qlf_WT$table),p.value=0.05)
table(sign(DEG_WT$table$logFC))
#-1   1 
#66 239 

# Identify DEGs - PMM2Stim_vs_PMM2NStim
qlf_PMM2 <- glmQLFTest(fit, contrast = contrasts[,2])
DEG_PMM2 <-topTags(qlf_PMM2, n=nrow(qlf_PMM2$table),p.value=0.05)
table(sign(DEG_PMM2$table$logFC))
#-1   1 
#34 188 

#############################################################################################
###################Writing Output Table with the Differently Expression Genes################
#############################################################################################



######################## Primary Biological Questions##################
#Separated DEGs based on upregulation and downregulation & Sample based
#########################PMM2 vs PMM2WT ############################

#install.packages(openxlsx)
library(openxlsx)

separate_positive_negativeDEG <- function(df){
  # Function that creates variables with DEG upregulated and downregulated
  DEG_positive <- vector()
  DEG_negative <- vector()
  for (i in 1:length(rownames(df))){
    ifelse(sign(df[i,1]) == +1,DEG_positive <- c(DEG_positive, rownames(df)[i]), DEG_negative <- c(DEG_negative,rownames(df)[i]) )
  }
  return (list(pos = as.data.frame(DEG_positive),neg = as.data.frame(DEG_negative)))
}

separetadedDEGs <- separate_positive_negativeDEG(DEG_WT$table)




# write output table for WT VS PMM2-CDG

DEG <- cbind("GeneName"=rownames(logCPMs), logCPMs)
sep_DEGs <- list("WTStim_vs_WTNStim"=DEG,"Background_genes"=unique(rownames(y)),separetadedDEGs["pos"],separetadedDEGs["neg"])
names(sep_DEGs) <- c("DEG","Background_genes","DEGs_upregulated","DEGs_downregulated")
#write.xlsx(sep_DEGs, file = "Expression.xlsx")

#####################################################################################
###############################Expression plot#######################################
#####################################################################################

#Calculate LogFC
LogFC_WT <- topTags(qlf_WT, n=nrow(qlf_WT$table), p.value = 1.1)
LogFC_PMM2 <- topTags(qlf_PMM2, n=nrow(qlf_PMM2$table), p.value = 1.1)

# Common gene counts
allGenes_common <- intersect(rownames(LogFC_WT$table), rownames(LogFC_PMM2$table))
length(allGenes_common)

#DEG Names
CommonDEGs <- intersect(rownames(DEG_WT$table), rownames(DEG_PMM2$table))
Exc_PMM2 <- setdiff(rownames(DEG_PMM2$table), CommonDEGs)
Exc_WT <-setdiff(rownames(DEG_WT$table), CommonDEGs)
  


# Plotting

png("myplot.png", width=900, height=900)
plot(LogFC_WT$table[allGenes_common,1], LogFC_PMM2$table[allGenes_common,1], xlab="Log2FC Stim WT Vs NStim WT ", ylab="Log2FC Stim PMM2-CDG Vs NStim PMM2-CDG", pch=20, col="gray", main="All genes",
          cex.main = 2, cex = 2,cex.xlab = 1.4,cex.ylab=1.4)
points(LogFC_WT$table[Exc_WT,1], LogFC_PMM2$table[Exc_WT,1], col="chartreuse3", pch=20, cex = 2)
points(LogFC_WT$table[Exc_PMM2,1], LogFC_PMM2$table[Exc_PMM2,1], col="springgreen4", pch=20, cex = 2)
points(LogFC_WT$table[CommonDEGs,1], LogFC_PMM2$table[CommonDEGs,1], col="tomato3", pch=20, cex = 2)

legend("bottomright", fill=c("gray","springgreen4","chartreuse3","tomato3"), legend=c("NDEGs", "PMM2-CDG-Exclusive","WT-Exclusive", "Common"), bty="n", cex=1.35)

v_lim <- min(abs(DEG_WT$table[,1]))

abline(v=c(-v_lim,v_lim),h=c(-v_lim ,v_lim),lty=2)

dev.off()

###########################################################################################
###################################Volcano Plot############################################
###########################################################################################

# Produce plots


# Produce plots
library(gplots)
library(gridExtra)
plotData_WT <- as.data.frame(topTags(qlf_WT, n=nrow(qlf_WT$table)))

# Assign colors to upregulated and downregulated genes
plotData_WT$colour <- ifelse(plotData_WT$logFC > 0 & plotData_WT$FDR < 0.05,"springgreen4",
                          ifelse(plotData_WT$logFC < 0 & plotData_WT$FDR < 0.05, "tomato3","gray20"))


plotData_WT$colour <- factor(plotData_WT$colour, levels = c("tomato3", "gray20", "springgreen4"))


#Select DEGs that should be Highlighted in the Plot
DEGs_WT_Higlights <- rownames(DEG_WT$table[DEG_WT$table[,1]> 5 | DEG_WT$table[,1] < -2.2,])


# Create the volcano plot

volcano_WT <- ggplot(plotData_WT, aes(x = logFC, y = -log10(FDR), color = colour)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("tomato3" = "tomato3", "gray20" = "gray20", "springgreen4" = "springgreen4"),
                     labels = c("Downregulated", "Non Significant", "Upregulated"))+
  xlim(c(-6.5, 15)) +
  ylim(c(0, 4.2)) +
  geom_hline(yintercept = -log10(max(DEG_WT$table$FDR)), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(- min(abs(DEG_WT$table$logFC)),0.56), linetype = "dashed", color = "grey50") +
  labs(title = "WTStim Vs WTNStim", x = "Log2(Fold Change)", y = "-log10(FDR)") +
  guides(fill = 'transparent',colour = guide_legend(override.aes = list(size=6),label.theme = element_text(size = 15))) +
  geom_text_repel(data = plotData_WT[DEGs_WT_Higlights , ], 
                  aes(x = logFC, y = -log10(FDR), label = DEGs_WT_Higlights), hjust = 0, fontface= "bold", vjust = 0, color = "grey27")+
  theme( panel.background = element_rect(fill = "white"),
         panel.spacing = unit(1.5, "lines"),
         legend.position = "top",
         legend.key = element_blank(),
         plot.title = element_text(size =20),
         legend.box.background = element_blank(),
         panel.grid.major = element_line(colour = "gray87"), panel.grid.minor = element_line(colour = "gray87"),
         axis.line = element_line(colour = "black"), legend.background = element_rect(fill='transparent')
  )


dev.off()

#PMM2-CDG

plotData_PMM2 <- as.data.frame(topTags(qlf_PMM2, n=nrow(qlf_PMM2$table)))

# Assign colors to upregulated and downregulated genes
plotData_PMM2$colour <- ifelse(plotData_PMM2$logFC > 0 & plotData_PMM2$FDR < 0.05,"springgreen4",
                             ifelse(plotData_PMM2$logFC < 0 & plotData_PMM2$FDR < 0.05, "tomato3","gray20"))


plotData_PMM2$colour <- factor(plotData_PMM2$colour, levels = c("tomato3", "gray20", "springgreen4"))


#Select DEGs that should be Highlighted in the Plot
DEGs_PMM2_Higlights <- rownames(DEG_PMM2$table[DEG_PMM2$table[,1]> 5.2 | DEG_PMM2$table[,1] < -2.2,])
DEG_PMM2$table[,1]

# Create the volcano plot

volcano_PMM2 <- 
  
  ggplot(plotData_PMM2, aes(x = logFC, y = -log10(FDR), color = colour)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("tomato3" = "tomato3", "gray20" = "gray20", "springgreen4" = "springgreen4"),
                     labels = c("Downregulated", "Non Significant", "Upregulated"))+
  xlim(c(-7., 12.5)) +
  ylim(c(0, 3.8)) +
  geom_hline(yintercept = -log10(max(DEG_PMM2$table$FDR)), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-min(abs(DEG_PMM2$table$logFC)), min(abs(DEG_PMM2$table$logFC))), linetype = "dashed", color = "grey50") +
  labs(title = "PMM2Stim Vs PMM2NStim", x = "Log2(Fold Change)", y = "-log10(FDR)") +
  guides(fill = 'transparent',colour = guide_legend(override.aes = list(size=6),label.theme = element_text(size = 15))) +
  geom_text_repel(data = plotData_PMM2[DEGs_PMM2_Higlights , ], 
                  aes(x = logFC, y = -log10(FDR), label = DEGs_PMM2_Higlights), hjust = 0, fontface= "bold", vjust = 0, color = "grey27")+
  theme( panel.background = element_rect(fill = "white"),
         panel.spacing = unit(1.5, "lines"),
         legend.position = "top",
         legend.key = element_blank(),
         plot.title = element_text(size =20),
         legend.box.background = element_blank(),
         panel.grid.major = element_line(colour = "gray87"), panel.grid.minor = element_line(colour = "gray87"),
         axis.line = element_line(colour = "black"), legend.background = element_rect(fill='transparent')
  )


grid.arrange(volcano_WT, volcano_PMM2, ncol = 2)
ggsave("volcanoplots.png", arrangeGrob(volcano_WT, volcano_PMM2, ncol = 2), width = 15.6, height =7)


#################################
# Plot gene ontology enrichment #
#################################

# Requires the package 
############################################################################################
#################################Fold Enrichment Calcullus##################################
############################################################################################

library(readxl)
wd <- "C:/Users/pmgra/OneDrive/Documentos/FCT/UCIBIO/Glycoimmunology/Fibroblasts/OVA/Functional Annotation/Refined Excels"
setwd(wd)


Enr_calc<- function(GO,allDEG,allNDGEs){
  Enr_list <- c("")
  for(i in 1:length(rownames(GO))){
    DEGs_in <- as.numeric(GO[i,1])
    NDEGs_in <- as.numeric(GO[i,2])-as.numeric(GO[i,1])
    DEGs_out <- as.numeric(allDEG - as.numeric(GO[i,1]))
    NDEGs_out <- as.numeric(allNDGEs - NDEGs_in)
    OD <- (DEGs_in*NDEGs_out)/(DEGs_out*NDEGs_in)
    Enr_list <- c(Enr_list,OD)
    return(Enr_list)
    }'ggplot2' (needs to be installed first)
# Load the ggplot2 package
library(ggplot2)

# set the working directory where the tables to use are located
setwd("C:/Users/pmgra/OneDrive/Documentos/CDG_Allies/Carlota's Data/Gene Ontology/New Approach/Refined Excels")


#######################################################
######Plot the representative GO terms by condition####
# or list of DEGs and compare the enrichment with #####
#######################################################
#######################################################

# Prepare dataframe
#------------------
# Import the table containing the enriched GO terms by groups
setwd("C:/Users/pmgra/OneDrive/Documentos/CDG_Allies/Carlota's Data/Gene Ontology/New Approach/Refined Excels")
GO_gp <- read.table("Fold Enrichment Plot.txt",header=T,stringsAsFactors = T)

# List objects and their structure contained in the dataframe 'GO_gp'
ls.str(GO_gp)

# Transform the column 'Gene_number' into a numeric variable
GO_gp$Gene_number <- as.numeric(GO_gp$Gene_number)

# Replace all the "_" by a space in the column containing the GO terms
GO_gp$GO_biological_process <- chartr("_", " ", GO_gp$GO_biological_process)

# Transform the column 'GO_biological_process' into factors
GO_gp$GO_biological_process<-as.factor(GO_gp$GO_biological_process)

# Change factor order
GO_gp$Group<- factor(GO_gp$Group,levels = c("WT_up","PMM2-CDG_up"))
GO_gp$GO_biological_process<-factor(GO_gp$GO_biological_process,levels=rev(levels(GO_gp$GO_biological_process)))

# Create a vector with new names for groups to use in the plot
# Replace the terms by your own (\n allow to start a new line)
group.labs <- c(`WT_up` = "WT upregulated",
                `PMM2-CDG_up` = "PMM2-CDG upregulated")

# Draw the plot in facets by group with ggplot2
# to represent -log10(FDR), Number of genes and 
# Fold enrichment of each GO biological process per group (Figure 3)
#-------------------------------------------------------------------
base_size <- 14  # Increase this value to increase the size of all text in the plot
line_size <- 1.5  # Increase this value to increase the size of lines in the plot

ggplot(GO_gp, aes(x = GO_biological_process, y = Fold_enrichment)) +
  geom_hline(yintercept = 1, linetype="dashed", color = "azure4", size=line_size) +
  geom_point(aes(size = Gene_number, colour = `FDR`)) +  # Increase the size of the points
  geom_density(alpha = 0.2) +
  coord_flip() +
  theme_bw(base_size = base_size) +  # Specify the base size
  theme(
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text.x = element_text(margin = margin(5, 5, 0, 5, "pt"), size = base_size),  # Increase the size of axis text
    axis.text.y = element_text(margin = margin(5, 5, 5, 5, "pt"), size = base_size),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.title.align = 0.5,
    plot.title = element_text(size = base_size * 1.2),  # Increase the size of the plot title
    legend.text = element_text(size = base_size)  # Increase the size of the legend text
  ) +
  xlab("GO biological processes") +
  ylab("Fold enrichment") +
  labs(color = "FDR ≤ 0.05", size = "Number\nof genes") +
  facet_wrap(~Group, ncol = 2, labeller = as_labeller(group.labs)) +
  guides(y = guide_axis(order = 2), colour = guide_legend(override.aes = list(alpha = 1)))

  # Save the plot as a PNG file with specified dimensions and resolution
ggsave(filename = "Enrichment_plot.png", dpi = 300)


#############################################################################################
########################Heatmap Molecular Interactions#######################################
#############################################################################################
read_excel("TNFA_stim.xlsx")
library(gplots)

Interactions <- read_excel("Interactions.xlsx", col_names =TRUE, sheet =2)
View(Interactions[,"Fibroblasts"])
selDEG <- unique(Interactions[,"Fibroblasts"])


###############################################################################
#############################HeatMaps Plot#####################################
###############################################################################


logFC <- cbind(DEG$table[rownames(DEG)[order(rownames(DEG))],1],LogFC_PMM2$table[rownames(DEG)[order(rownames(DEG))],1])
rownames(logFC) <-rownames(DEG)[order(rownames(DEG))]
colnames(logFC) <- c("WT","PMM2-CDG")


#Immunological Response
library(RColorBrewer); library(gplots)
plotCol_exp <- brewer.pal(9, "Greens")
plotColSamples <- c(rep("grey50",3),rep("grey70", 2))
par(mar = rep(2, 4))
Genes_interest <- c("IL1B", "IL6", "IL15", "CXCL1", "CXCL5", "CXCL8", "CCL2", "CCL5", "NFKB1","NFKB2", "NFKBIA")

length(Genes_interest)

logCPMs["IFNG",]

Expression <- logCPMs[rownames(DEG$table),]
heatmap.2(logCPMs[Genes_interest,], Rowv=FALSE,  scale="row", dendrogram = "none", cexCol = 1, cexRow = 0.57, trace="none", density.info="none", ColSideColors = plotColSamples, col=plotCol_exp,main="Immunological Response",margins = c(0.01,25))




###############################################################################3
wd <- "C:/Users/pmgra/OneDrive/Documentos/FCT/FCT/UCIBIO/Glycoimmunology/Fibroblasts/Transcriptomic Code"
setwd(wd)


library(readxl)
library(tidyr)
library(corrplot)
library(RColorBrewer)


Interactions <- as.data.frame(read_excel("PredictedInteraction_WT_PMM2.xlsx", col_names =TRUE))

selDEG <- as.vector(unique(Interactions[,"DEGs"]))
split_elements <-as.vector(unlist(strsplit(selDEG, ",")))
updated_string <- gsub(" ", "", split_elements)
uniqueDEG <- na.omit(unique(updated_string))


LogFC_WT <- topTags(qlf_WT, n=nrow(qlf_WT$table), p.value = 1.1)
LogFC_PMM2 <- topTags(qlf_PMM2, n=nrow(qlf_PMM2$table), p.value = 1.1)

Classification <- na.omit(data.frame(Interactions[,c("DEGs","Condition")]))
Classification$Condition <- factor(Classification$Condition, levels = c("Common", "WT", "PMM2"))
Classification <- Classification %>% arrange(Condition)

expgenes <- data.frame(LogFC_WT[uniqueDEG,c("logFC","FDR")],LogFC_PMM2[uniqueDEG,c("logFC","FDR")])
colnames(expgenes) <- c("LogFC_WT","FDR_WT","LogFC_PMM2","FDR_PMM2")



library(tidyverse)

# Convert 'Classification' column to a factor and specify the levels

# Transform data from wide to long format
expgenes_long <- expgenes %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Variable", value = "Value", -Gene)

# Separate condition and metric
expgenes_long <- expgenes_long %>%
  separate(Variable, into = c("Condition", "Metric"), sep = "_")



# Remove FDR rows and handle NAs
expgenes_heatmap <- expgenes_long %>%
  filter(Condition == "LogFC") %>%
  spread(key = Metric, value = Value)


# Convert the data frame to a matrix
expgenes_matrix <- data.matrix(expgenes_heatmap[,-1])
rownames(expgenes_matrix) <- expgenes_heatmap[,1]



# Get FDR data
expgenes_FDR <- expgenes_long %>%
  filter(Condition == "FDR") %>%
  spread(key = Metric, value = Value)



# Convert the data frame to a matrix
expgenes_FDR_matrix <- data.matrix(expgenes_FDR[,-1])
rownames(expgenes_FDR_matrix) <- expgenes_FDR[,1]



# Function to convert FDR to asterisks
FDR_to_star <- function(fdr) {
  if (fdr < 0.001) {
    return('***')
  } else if (fdr < 0.01) {
    return('**')
  } else if (fdr < 0.05) {
    return('*')
  } else {
    return('')
  }
}

# Apply function to FDR matrix
# Apply function to each value in FDR matrix
expgenes_star_matrix <- apply(expgenes_FDR_matrix, c(1, 2), FDR_to_star)

expgenes_matrix <- expgenes_matrix[,-1]
expgenes_star_matrix <- expgenes_star_matrix[, -1]
  




Interactions_Immune <- Interactions[,8:ncol(Interactions)]


Interactions_Immune <- Interactions_Immune %>%
  filter(`Present in immune cell?` == 1)


# Convert 'Fibroblasts' column to a list column where each cell contains all gene names
Interactions_Immune$Fibroblasts <- lapply(strsplit(Interactions_Immune$Fibroblasts, ",\\s*"), trimws)

# Unnest the data frame so that there is one row per gene
library(tidyr)
Interactions_Immune <- as.data.frame(unnest(Interactions_Immune, Fibroblasts))


Interactions_Immune <- Interactions_Immune[,-(match(c("Other cell","Select interactions","Present in immune cell?",	"Which?"),colnames(Interactions_Immune)))]

Interactions_Immune <- Interactions_Immune <- Interactions_Immune %>%
  group_by(Fibroblasts) %>%
  summarise_all(.funs = sum, na.rm = TRUE)

Interactions_Immune <- as.data.frame(Interactions_Immune %>%
  group_by(Fibroblasts) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)))

rownames(Interactions_Immune) <- Interactions_Immune$Fibroblasts

Interactions_Immune <- Interactions_Immune[,-match("Fibroblasts", colnames(Interactions_Immune))]
dim(Interactions_Immune)

Interactions_Immune <- Interactions_Immune[Classification$DEGs,]

# Remove rows with NA in rownames
Interactions_Immune <- Interactions_Immune[!grepl("^NA", rownames(Interactions_Immune)), ]




# Convert rownames to a column in Interactions_Immune
Interactions_Immune <- Interactions_Immune %>%
  rownames_to_column("DEGs")

# Add the "Condition" column from Classification to Interactions_Immune
Interactions_Immune <- Interactions_Immune %>%
  left_join(Classification[ , c("DEGs", "Condition")], by = "DEGs")

# Define the color mapping
color_mapping <- c("Common" = "#31915A", "WT" = "#C2C2C2", "PMM2" = "#40A068")

# Add a "Color" column based on the "Condition" column

# Initialize color matrix with the same dimensions as Interactions_Immune
color_matrix <- matrix(nrow = nrow(Interactions_Immune), ncol = (ncol(Interactions_Immune)-3))
col <- c("tomato","tomato","tomato","tomato","goldenrod1","#C7E9B4","#C7E9B4","#C7E9B4","#C7E9B4","#C7E9B4","#C7E9B4","#7FCDBB","#7FCDBB","#41B6C4" ,"#1D91C0", "#225EA8","#0C2C84")
length(col)

# Assign colors to matrix positions based on interactions
for (i in 1:nrow(Interactions_Immune)) {
  for (j in 2:(ncol(Interactions_Immune)-2)) {
    if (as.numeric(Interactions_Immune[i, j]) > 0) {  # replace this condition with your own if needed
      color_matrix[i, j-1] <- col[j-1]
    } else {
      color_matrix[i, j-1] <- NA
    }
  }
}


rownames(color_matrix) <- Interactions_Immune$DEGs


expgenes_star_matrix <- expgenes_star_matrix[Interactions_Immune$DEGs,]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

color_gradient <- colorRampPalette(brewer.pal(9, "Greens"))

# Create heatmap

ldg = list(direction = "vertical")

heatmap <- Heatmap(expgenes_matrix[Interactions_Immune$DEGs,],
                   name = "logFC",
                   col = color_gradient(9),
                   show_column_names = TRUE, 
                   show_row_names = TRUE,
                   row_names_side = "left",
                   heatmap_legend_param = ldg,
                   cluster_rows = FALSE,
                   width = rep(0.01, ncol(expgenes_matrix[Interactions_Immune$DEGs,])),
                   cluster_columns = FALSE,
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     # Draw the cell
                     grid.rect(x, y, w, h, gp = gpar(fill = fill))
                     
                     # Add the significance star
                     star <- expgenes_star_matrix[i, j]
                     if (star != "") {
                       grid.text(star, x , y*0.995)
                     }
                   })



draw(heatmap, heatmap_legend_side = "left")


# Set up the plotting device

pin_value <- 5  # Adjust this value
nr <- nrow(color_matrix)
nc <- ncol(color_matrix)
par(pin=c(pin_value, pin_value*nr/nc), oma=c(5, 4, 0, 0), mar=c(0, 0, 0, 0))



# create a blank plot
plot(1, 1, xlim=c(1, nc), ylim=c(1, nr), type="n", xaxt='n', yaxt='n', xlab="", ylab="", asp=1)

# Define the size of the gap
gap <- 0.15

# add colored squares based on the presence or absence of color in the matrix
for (i in 1:nr) {
  for (j in 1:nc) {
    if (!is.na(color_matrix[i, j])) {
      rect(j-1+gap, nr-i+gap, j-gap, nr-i+1-gap, col=color_matrix[i, j], border=NA)
    }
  }
}

box(lty = 0)
