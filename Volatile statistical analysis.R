# Author: Christina Papazlatani
# Date: 29-04-2024


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#         Data Preprocessing   #-#-#-# ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# (I) Extract the peak intensity table from the MZmine in the metaboanalyst form

# Manually curate the Peak Intensity Table
## Remove silica compounds
## Remove artifacts
## Replace missing values with 0
## Transpose the table


## Table example

# Sample	                TSBAControlA.mzdata.xml	TSBAControlB.mzdata.xml	TSBAControlC.mzdata.xml
# Groups	                TSBA_CTL	              TSBA_CTL	              TSBA_CTL
# C1/44.0572mz/1.27min	  4.85E+04	              1.72E+05	              2.98E+04
# C2/40.0339mz/1.37min	  2.84E+05	              4.69E+05	              3.35E+04


# (II) In the project directory, create the folder "Volatilomic analysis" and then prepare the following folders (Copy/Paste to make sure the names are the same)
# "0_Data_preparation"
# "1_Data_norm"
# "2_Diff.Int"
# "3_Burk_Comp_Norm"
# "4_Diff_Int_Media"


# (IΙΙ) Place the scripts in the "Volatilomic analysis" folder and the peak intensity table in the "0_Data_preparation" folder


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#      1st   Data Transformation / Normalization   #-#-#-# ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Set the working directory
library(here)
setwd(paste0(here(), "/Volatilomic analysis/1_Data_norm/") )

# Load the necessary packages
package.v <- list()
package.v[["R Version"]] <- R.version$version.string
package.v[["R Studio Version"]] <- RStudio.Version()$version
library(openxlsx); package.v[["openxlsx"]] <- packageVersion("openxlsx")
library(functional); package.v[["functional"]] <- packageVersion("functional") #Necessary for pretreat command
capture.output(package.v, file = "0_Package versions.txt")

# Load the data
list.files("../0_Data_preparation/")
mix <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "AAmixTransp", startRow = 1, colNames = T, rowNames = T)
dim(mix)
# 12 147
gly <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "GlyTransp", startRow = 2, colNames = T, rowNames = T)
dim(gly)
# 12 116
gln <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "GlnTransp", startRow = 2, colNames = T, rowNames = T)
dim(gln)
# 12 136
asn <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "AspTransp", startRow = 2, colNames = T, rowNames = T)
dim(asn)
# 12 100
arg <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "ArgTransp", startRow = 2, colNames = T, rowNames = T)
dim(arg)
# 12 99

# Set the amino acid treatments
myaminos <- c("mix","gly","gln", "asn", "arg")

for(myaa in myaminos) {
  
  if(myaa == "gly") {
    mydata <- gly
    mytitle <- "Glycine"
    
  } else if( myaa == "gln") {
    mydata <- gln
    mytitle <- "Glutamine"
    
  } else if (myaa == "asn") {
    mydata <- asn
    mytitle <- "Asparagine"
    
  } else if (myaa == "arg") {
    mydata <- arg
    mytitle <- "Arginine"
    
  } else if (myaa == "mix") {
    mydata <- mix
    mytitle <- "Amino acid mixture"
  }
  
  # Get the compounds table
  mydt <- mydata[,-1]
  
  # Get the metadata table
  metadata <- data.frame(mydata[,1])
  rownames(metadata) <- rownames(mydata)
  colnames(metadata) <- "Groups"
  
# Data pretreatment methods
# The normalization procedures are grouped into three categories. You can use one or combine them to achieve better results.

# (1) Centering converts all the concentrations to fluctuations around zero instead of around the mean of the metabolite concentrations. Hereby, it adjusts for differences in the offset between high and low abundant metabolites. It is therefore used to focus on the fluctuating part of the data and leaves only the relevant variation (being the variation between the samples) for analysis. Centering is applied in combination with all the methods described below.

# (2) Transformations are nonlinear conversions of the data. They are generally applied to correct for heteroscedasticity, to convert multiplicative relations into additive relations, and to make skewed distributions (more) symmetric.
#           1) Log2 Transformation, 
#           2) Square Root Transformation,
#           3) Cubic Root Transformation

# (3) Scaling methods are data pretreatment approaches that divide each variable by a factor, the scaling factor, which is different for each variable. They aim to adjust for the differences in fold differences between the different metabolites by converting the data into differences in concentration relative to the scaling factor..
#   Scaling based on data dispersion:
#           1) Autoscaling
#           2) Pareto scaling
#           3) Range scaling
#           4) Vast scaling
#   Scaling based on average value:
#           1) Level scaling


  # pretreat command takes as variables
  # transf = "log2" / "log10 / "sqrt" / "none"
  # center = TRUE / FALSE
  # scale = "none" / auto" / "pareto"
  
  source("../pretreat command V2.1.R")
  pretreat_mydt <- pretreat(mydt, transf = "sqrt", center = T, scale = "pareto")
  mytrasf <- pretreat_mydt$mytrasf
  transf_mydt <- pretreat_mydt$transf_mydt

  ### Visualize the result of data transformation - Calculate the means of each compound/column OR just plot the values (which makes more sense)
  valuesbefore <- unlist(mydt)
  valuesafter <- unlist(transf_mydt)

  # Prepare the plot
  pdf(paste0("Data_Norm_", mytitle,"_", mytrasf, ".pdf"), height = 4, width = 6)
  
  par(mfrow= c(1,2), mar = c(4,4,4,2), cex.main = 1.5)
  plot(density(valuesbefore),
       col = "blue", lwd = 2,
       main = paste0(mytitle," Before Transformation")
       , xlab = "")
  plot(density(valuesafter),
       col = "green", lwd = 2,
       main = paste0(mytitle," After Transformation")
       , xlab = "")
  dev.off()
  
  #  Merge the normalized values table with the metadata file...
  normtbl <- merge(transf_mydt, metadata, by = "row.names")
  rownames(normtbl) <- normtbl$Row.names 
  normtbl <- normtbl[, 2:ncol(normtbl)]
  
  # ...and save it
  write.table(normtbl, paste0("Values_Norm_",mytitle,"_",mytrasf,".txt"), col.names = NA, quote = FALSE, sep = "\t", row.names = T)
  
}


#Clean your environment and continue with the following step
  



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#      Pairwise comparison between abiotic controls and inoculated treatments for each medium #-#-#-# ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
  
# Set the working directory
library(here)
setwd(paste0(here::here(), "/Volatilomic analysis/2_Diff_Int_Inocul/") )

# Load the necessary packages
package.v <- list()
package.v[["R Version"]] <- R.version$version.string
package.v[["R Studio Version"]] <- RStudio.Version()$version
library(effsize); package.v[["effsize"]] <- packageVersion("effsize")
library(plyr); package.v[["plyr"]] <- packageVersion("plyr")
library(vegan); package.v[["vegan"]] <- packageVersion("vegan")
library(RColorBrewer); package.v[["RColorBrewer"]] <- packageVersion("RColorBrewer")
library(pheatmap); package.v[["pheatmap"]] <- packageVersion("pheatmap")
library(ggforce); package.v[["ggforce"]] <- packageVersion("ggforce")
library(ggpubr); package.v[["ggpubr"]] <- packageVersion("ggpubr")
capture.output(package.v, file = "0_Package versions.txt")

# Load the data
list.files("../1_Data_norm/")
mix <- read.table( "../1_Data_norm/Values_Norm_Amino acid mixture_SqRt Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(mix)
# 12 147
gly <- read.table( "../1_Data_norm/Values_Norm_Glycine_SqRt Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(gly)
# 12 116
gln <- read.table( "../1_Data_norm/Values_Norm_Glutamine_SqRt Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(gln)
# 12 136
asn <- read.table( "../1_Data_norm/Values_Norm_Asparagine_SqRt Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(asn)
# 12 100
arg <- read.table( "../1_Data_norm/Values_Norm_Arginine_SqRt Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(arg)
# 12 99

# Set the treatments
myaminos <- c("mix","gly","gln", "asn", "arg")
heatmap_list <- list()

for(myaa in myaminos) {
  
  if(myaa == "gly") {
    analysis_tbl <- gly
    mytitle <- "Glycine"
    
  } else if( myaa == "gln") {
    analysis_tbl <- gln
    mytitle <- "Glutamine"
    
  } else if (myaa == "asn") {
    analysis_tbl <- asn
    mytitle <- "Asparagine"
    
  } else if (myaa == "arg") {
    analysis_tbl <- arg
    mytitle <- "Arginine"
    
  } else if (myaa == "mix") {
    analysis_tbl <- mix
    mytitle <- "Amino acid mixture"
  }

  ## Set the inoculation scheme variable
  analysis_tbl$Inoculation <- rep(c("Burkholderia sp. AD24","Abiotic Control"), each = 3, times = 2 )
  # levels(as.factor(analysis_tbl$Inoculation))
  analysis_tbl$Media <- rep(c("0.1 TSBA + AA","0.1 TSBA"), each = 6)
  # levels(as.factor(analysis_tbl$Media))
  
  ## We will run the analysis for each medium
  mymedia <- levels(as.factor(analysis_tbl$Media))
  
  # Prepare the colors of the barplots
  mycols <- brewer.pal(4, "Set1")
  mycols <- mycols[unique(as.numeric(as.factor(analysis_tbl$Groups)))]
  
  
  ## We will run the analysis for each compound
  mycomps <- colnames(analysis_tbl)[1: c(ncol(analysis_tbl)-3)]
  print(length(mycomps))
  
  # Prepare a stats list to save all the statistic results and the boxplots and a table for the final results
  stats_list <- list()
  groups_tbl <- data.frame()
  identical_intensities <- data.frame()
  plot_list <- vector("list", 2*length(mycomps))
  names(plot_list)<- paste(rep(mymedia, length(mycomps)), mycomps)
  
  for(mymedium in mymedia) {
    tbl_sel <- analysis_tbl[analysis_tbl$Media == mymedium , ]
    
    ### For analysis of variance
    for( mycomp in mycomps) {
      
      if(length(unique(tbl_sel[,mycomp])) == 1 ) {
        
        # Prepare the data table
        my_groups <- plyr::ddply(tbl_sel, c("Inoculation"), summarise,
                                 Means = mean(get(mycomp), na.rm = T),
                                 sd = sd(get(mycomp), na.rm = T),
                                 se = sd/sqrt(length(get(mycomp))),
                                 Compound = mycomp,
                                 Medium = mymedium )
        
        # Calculate the logFC as well
        my_groups$diff <- my_groups$Means[2] - my_groups$Means[1] 
        my_groups$perc.rel.change <- ( (my_groups$diff)/ my_groups$Means[2] ) * 100
        identical_intensities <- rbind(identical_intensities, my_groups)
        
        
      } else {
        myttest <- t.test(tbl_sel[,mycomp]  ~ Inoculation, data = tbl_sel)
        stats_list[[paste(mymedium, mycomp)]][["Student's t-Test"]] <- myttest
        
        # Calculate the effect size
        ef.size <- cohen.d(tbl_sel[,mycomp], tbl_sel$Inoculation, na.rm = T)
        stats_list[[paste(mymedium, mycomp)]][["Effect size: Cohen's d"]] <- ef.size
        
        # Prepare the data table
        my_groups <- plyr::ddply(tbl_sel, c("Inoculation"), summarise,
                                 Means = mean(get(mycomp), na.rm = T),
                                 sd = sd(get(mycomp), na.rm = T),
                                 se = sd/sqrt(length(get(mycomp))),
                                 test = myttest$method,
                                 p.value = myttest$p.value,
                                 eff.size = ef.size$estimate,
                                 CI.lower = ef.size$conf.int[[1]],
                                 CI.upper = ef.size$conf.int[[2]],
                                 label = paste("Effect size:", round(ef.size$estimate, 2) ),
                                 Compound = mycomp,
                                 Medium = mymedium )
        
        # Calculate the logFC as well
        my_groups$diff <- my_groups$Means[2] - my_groups$Means[1] 
        my_groups$perc.rel.change <- ( (my_groups$diff)/ my_groups$Means[2] ) * 100
        groups_tbl <- rbind(groups_tbl, my_groups)
        
        ## Prepare a boxplot with the results
        message(paste(mymedium,mycomp))
        plot_list[[paste(mymedium,mycomp)]] <- local({
          mycomp <- mycomp
          
          dt_sel <- tbl_sel[,c("Groups",mycomp)]
          mytreatorder <- levels(as.factor(tbl_sel$Groups))
          dt_sel$Groups <- factor(dt_sel$Groups, levels = mytreatorder)
          
          myplot <-
            ggplot(dt_sel, aes(x = Groups, y = dt_sel[, mycomp], fill = Groups) ) +
            geom_boxplot(fill = mycols[unique(as.numeric(dt_sel$Groups))]) +
            geom_jitter(shape = 21, size=5, alpha=0.5) +
            scale_fill_manual(values = mycols[unique(as.numeric(dt_sel$Groups))]) +
            labs(title = paste(mymedium,mycomp), x= "Treatment", y = "Normalized Intensities") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust=1),
                  plot.background = element_rect(colour = "black", fill=NA, linewidth =1))
          print(myplot)
        })
        
      }
    }
  }
  
  capture.output(stats_list, file = paste0("1_", toupper(myaa), "_stats.txt"))

  #--#--# Tables with information 
  # Retrieve the significant compounds...
  sign_tbl <-  groups_tbl[groups_tbl$p.value < 0.05,]
  sign_tbl$media.comp <- paste(sign_tbl$Medium, sign_tbl$Compound, sep = ", ")
  # sign_comps <- unique(sign_tbl$media.comp)
  ## Save the compounds table
  write.table(sign_tbl, file = paste0("2_",toupper(myaa),"_SignComp_stat.table.txt"), col.names = NA, quote = FALSE, sep = "\t")
  
  
  # Select the compounds that are significantly present in the INOCULATED treatments...
  Burk_Comps_tbl <- groups_tbl[groups_tbl$diff >= 0 & groups_tbl$p.value < 0.05 ,]
  Burk_comps <- unique(Burk_Comps_tbl$Compound)
  myBurkdata <- analysis_tbl[ ,colnames(analysis_tbl) %in% c(Burk_comps, "Media")]
  ## ...and save it
  write.table(myBurkdata, paste0("3_",toupper(myaa),"_BurkNormValues.txt"), col.names = NA, quote = FALSE, sep = "\t")
  
  
  
  #--#--# Boxplots with the significant compounds
  plot_positions <-   sapply(Burk_comps, function(x) grep(x, names(plot_list)))
  sign_plot <- plot_list[plot_positions]
  # length(sign_plot)
  
  pdf(paste0("4_SignComp_boxplots.pdf"), height = 10, width = 9)
  mybarplots <- ggarrange(plotlist = sign_plot, ncol=2, nrow = 3, legend = "none")
  print(mybarplots)
  dev.off()
  
  dev.off()
  
  
  
  #--#--# Heatmaps with the normalized intensities of the significant compounds
  
  # (I) Retrieve the significant compounds... 
  sign_comps <- unique(sign_tbl$Compound)
  print(length(sign_comps))
  sign_intensities <- analysis_tbl[,colnames(analysis_tbl) %in% c(sign_comps, "Groups")]
  ## ...and prepare the intensities table for the heatmap
  sign_intensities_t <- t(sign_intensities[,-ncol(sign_intensities)]) # Transpose
  sign_intensities_t <- decostand(sign_intensities_t, MARGIN = 1, method = "range") # Transform between 0-1

  
  # (II) Get the p-value annotation column
  ## For TSBA
  groups_tbl_tsba <- groups_tbl[groups_tbl$Medium == "0.1 TSBA", ]
  tsba <-
    data.frame( "Compounds" = groups_tbl_tsba$Compound[seq(from = 1, to = nrow(groups_tbl_tsba), by = 2)],
                "p.value" = groups_tbl_tsba$p.value[seq(from = 1, to = nrow(groups_tbl_tsba), by = 2)]
    )
  tsba$color <- ifelse(tsba$p.value < 0.05, "< 0.05", "> 0.05")
  tsba <- tsba[tsba$Compounds %in% sign_comps, ]
  ## For glutamine only - Include the compounds whose identities were the same between the two inculation treatments on 0.1 TSBA medium
  if(myaa == "gln") {
    tsba <- rbind(tsba, list("X67.61.0362mz.5.96min", "1", "> 0.05") )
    tsba <- rbind(tsba, list("X152.61.0385mz.9.72min", "1", "> 0.05") )
    tsba <- tsba[match(colnames(sign_intensities)[-ncol(sign_intensities)], tsba$Compounds),]
  }
  ## For TSBA + AA
  groups_tbl_tsba_aa <-  groups_tbl[groups_tbl$Medium == "0.1 TSBA + AA", ]
  tsba_aa <-
    data.frame( "Compounds" = groups_tbl_tsba_aa$Compound[seq(from = 1, to = nrow(groups_tbl_tsba_aa), by = 2)],
                "p.value" = groups_tbl_tsba_aa$p.value[seq(from = 1, to = nrow(groups_tbl_tsba_aa), by = 2)]
    )
  tsba_aa$color <- ifelse(tsba_aa$p.value < 0.05, "< 0.05", "> 0.05")
  tsba_aa <- tsba_aa[tsba_aa$Compounds %in% colnames(sign_intensities), ]
  
  
  # (ΙΙΙ) Prepare the experimental design table 
  my_treat_groups <- data.frame("Treatment" = rep(c("Burkholderia AD24","Abiotic Control"), each = 3, times = 2),
                                "Medium" = rep(c("0.1 TSBA + AA", "0.1 TSBA") , each = 6),
                                "row.names" = colnames(sign_intensities_t)  )
  
  
  # (IV) Prepare the colors for the heatmap.... 
  mycolors <- brewer.pal(9, "OrRd")
  
  # ....and for the annotations....
  mytreatcls <- brewer.pal(4, "Set1")
  mypvalcls <- brewer.pal(4, "Dark2")
  
  #....and put them in a list
  myhtmpcls <- list(Treatment = c( "Burkholderia AD24" = mytreatcls[3], "Abiotic Control" = mytreatcls[4] ),
                    Medium = c( "0.1 TSBA" = mytreatcls[1], "0.1 TSBA + AA" = mytreatcls[2]),
                    TSBA_pvals = c( "< 0.05" = mypvalcls[4], "> 0.05" = mypvalcls[1] ) ,
                    TSBA_AA_pvals = c( "< 0.05" = mypvalcls[4], "> 0.05" = mypvalcls[1]  )
  )
  
  # Prepare the row annotation table
  row_annotation <- data.frame("row.names" = rownames(sign_intensities_t),
                               "TSBA_AA_pvals" = tsba_aa$color,
                               "TSBA_pvals" = tsba$color)
  
  # (V) Cluster the rows with the cluster package
  Ks=sapply(2:(ncol(t(sign_intensities_t))-1),
            function(i)
              summary(cluster::silhouette(cluster::pam(sign_intensities_t,k=i)))$avg.width)
  names(Ks) <- 2:(ncol(t(sign_intensities_t))-1)
  mycuttreerows <- as.numeric(names(which(Ks == max(Ks))))
  # mycuttreerows

  
  # (VI) Prepare the heatmap
  heatmap_list[[myaa]] <- 
    pheatmap(sign_intensities_t, 
             color = mycolors,
             fontsize_row = 8,
             main = mytitle,
             cluster_cols = F, gaps_col = 6,
             cluster_rows = T, cutree_rows = mycuttreerows,
             annotation_col = my_treat_groups,
             annotation_row = row_annotation,
             annotation_colors = myhtmpcls,
             border_color = "grey50", cellwidth = 15, cellheight = 8)
}

pdf("Supplementary Figure 2 - Heatmaps Burk vs CTL.pdf", height = 9, width = 9, onefile = T)

ggarrange(heatmap_list[[1]][[4]], 
          heatmap_list[[2]][[4]],
          heatmap_list[[3]][[4]],
          heatmap_list[[4]][[4]],
          heatmap_list[[5]][[4]],
          ncol = 1, nrow = 1)

dev.off()


legend <- c("Heatmap of transformed intensities (square root transformed, pareto scaled and ranged between 0 and 1) of the compounds that were found to be significantly different between inoculated treatments and abiotic controls. Compound names in red indicate that they were found in significantly higher intensities in inoculated media. Left most columns indicate p-value of the pairwise comparisons for each compound in 0.1 TSA (right) and 0.1 TSA + AAmixture (left).")
capture.output(legend, file = "6_Htmp-SqRt,Centering Intens_legend.txt")

#Clean your environment and continue with the following step



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#      Remove control from Burk Compounds and Normalize again    #-#-#-# ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

## Set the working directory
library(here)
setwd(paste0(here::here(), "/Volatilomic analysis/3_Burk_Comp_Norm/") )

## load the necessary packages
package.v <- list()
package.v[["R Version"]] <- R.version$version.string
package.v[["R Studio Version"]] <- RStudio.Version()$version
capture.output(package.v, file = "0_Package versions.txt")


# Load the data
list.files("../0_Data_preparation/")
mix <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "AAmixTransp", startRow = 1, colNames = T, rowNames = T)
dim(mix)
# 12 147
gly <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "GlyTransp", startRow = 2, colNames = T, rowNames = T)
dim(gly)
# 12 116
gln <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "GlnTransp", startRow = 2, colNames = T, rowNames = T)
dim(gln)
# 12 136
asn <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "AspTransp", startRow = 2, colNames = T, rowNames = T)
dim(asn)
# 12 100
arg <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "ArgTransp", startRow = 2, colNames = T, rowNames = T)
dim(arg)
# 12 99

# Set the amino acid treatments
myaminos <- c("mix","gly","gln", "asn", "arg")

for(myaa in myaminos) {
  
  if(myaa == "gly") {
    mydata <- gly
    mytitle <- "Glycine"
    Burknormvalues <- read.table(file = "../2_Diff_Int_Inocul/3_GLY_BurkNormValues.txt", header = T, row.names = 1, sep = "\t", check.names = T)
    
  } else if( myaa == "gln") {
    mydata <- gln
    mytitle <- "Glutamine"
    Burknormvalues <- read.table(file = "../2_Diff_Int_Inocul/3_GLN_BurkNormValues.txt", header = T, row.names = 1, sep = "\t", check.names = T)

  } else if (myaa == "asn") {
    mydata <- asn
    mytitle <- "Asparagine"
    Burknormvalues <- read.table(file = "../2_Diff_Int_Inocul/3_ASN_BurkNormValues.txt", header = T, row.names = 1, sep = "\t", check.names = T)

  } else if (myaa == "arg") {
    mydata <- arg
    mytitle <- "Arginine"
    Burknormvalues <- read.table(file = "../2_Diff_Int_Inocul/3_ARG_BurkNormValues.txt", header = T, row.names = 1, sep = "\t", check.names = T)

  } else if (myaa == "mix") {
    mydata <- mix
    mytitle <- "Amino acid mixture"
    Burknormvalues <- read.table(file = "../2_Diff_Int_Inocul/3_MIX_BurkNormValues.txt", header = T, row.names = 1, sep = "\t", check.names = T)
  }
  
  # Get the names of the compounds found in significantly higher intensitites in the inoculated treatments
  BurkComps <- colnames(Burknormvalues)[-ncol(Burknormvalues)]

  # Fix the column names of the raw data
  colnames(mydata) <- gsub("/",".", colnames(mydata)) ; 
  colnames(mydata)  <- gsub("^","X", colnames(mydata), perl = T)
  
  # Select the bacterial compounds
  myBurkdata <- mydata[, colnames(mydata) %in% c(BurkComps, "XGroups")]
  print(dim(myBurkdata))
  # 12 23
  
  myBurkdata$Media <- gsub("CTL ", "" ,myBurkdata$XGroups)
  
  # Set the media variable
  mymedia <- levels(as.factor(myBurkdata$Media))
  
  # Remove the controls
  mydatalist <- list()
  for (mymedium in mymedia) {
    
    cleanData <- data.frame("Medium" = rep(mymedium, 3))
    
    for( BurkComp in BurkComps ) {
      
      myogBurkComps_sel <- myBurkdata[ myBurkdata$Media == mymedium , c("XGroups" ,"Media", BurkComp)]
      
      # Calculate the control mean
      CntrlMean <- mean(myogBurkComps_sel[ grep("CTL", myogBurkComps_sel$XGroups),3])
      
      #Substract the CntrlMean from the Burk intensities
      substrInten <- myogBurkComps_sel[ -grep("CTL", myogBurkComps_sel$XGroups),3] - CntrlMean
      
      cleanData <- cbind(cleanData, substrInten)
      colnames(cleanData)[ncol(cleanData)] <- BurkComp
    }
    
    mydatalist[[mymedium]] <- cleanData
  }
  
  
  # Combine the tables by column name
  mydatafin <- merge(mydatalist[[1]], mydatalist[[2]], all = T)
  rownames(mydatafin) <-  paste0(mydatafin$Medium,"_", rep(c("A","B","C"), times = 2))
  
  write.table(mydatafin, file = paste0("1_",toupper(myaa) ,"_Intensities without Control.txt"), quote = F, col.names = NA, sep = "\t")
  
  
  # Intensities Transformation
  mydt <- mydatafin[,-1]
  metadata <- data.frame(mydatafin[,1])
  rownames(metadata) <- rownames(mydatafin)
  colnames(metadata) <- "Groups"
  
  source("../pretreat command V2.1.R")
  pretreat_mydt <- pretreat(mydt, transf = "none", center = T, scale = "pareto")
  mytrasf <- pretreat_mydt$mytrasf
  transf_mydt <- pretreat_mydt$transf_mydt
  
  ### Visualize the result of data transformation - Calculate the means of each compound/column OR just plot the values (which makes more sense)
  valuesbefore <- unlist(mydt)
  shapiro.test(valuesbefore)
  
  valuesafter <- unlist(transf_mydt)
  shapiro.test(valuesafter)
  
  # Prepare the plot
  pdf(paste0("2_DataNorm_",toupper(myaa), "_", mytrasf, ".pdf"), height = 4, width = 6)
  
  par(mfrow= c(1,2), mar = c(4,4,4,2), cex.main = 1.5)
  plot(density(valuesbefore),
       col = "blue", lwd = 2,
       main = paste0(mytitle,"Before Transformation")
       , xlab = "")
  plot(density(valuesafter),
       col = "green", lwd = 2,
       main = paste0(mytitle, "After Transformation")
       , xlab = "")
  dev.off()
  
  #  Merge the normalized values table with the metadata file...
  normtbl <- merge(transf_mydt, metadata, by = "row.names")
  rownames(normtbl) <- normtbl$Row.names 
  normtbl <- normtbl[, 2:ncol(normtbl)]
  
  # ...and save it
  write.table(normtbl, paste0("3_NormInten_without_control_",toupper(myaa),"_", mytrasf,".txt"), col.names = NA, quote = FALSE, sep = "\t", row.names = T)
}



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#      Focus on Burk compounds: Paiwise comparison between media    #-#-#-# ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

## Set the working directory
library(here)
setwd(paste0(here::here(), "/Volatilomic analysis/4_Diff_Int_Media/") )

## load the necessary packages
package.v <- list()
package.v[["R Version"]] <- R.version$version.string
package.v[["R Studio Version"]] <- RStudio.Version()$version
library(effsize); package.v[["effsize"]] <- packageVersion("effsize")
library(plyr); package.v[["plyr"]] <- packageVersion("plyr")
library(vegan); package.v[["vegan"]] <- packageVersion("vegan")
library(RColorBrewer); package.v[["RColorBrewer"]] <- packageVersion("RColorBrewer")
library(pheatmap); package.v[["pheatmap"]] <- packageVersion("pheatmap")
library(ggforce); package.v[["ggforce"]] <- packageVersion("ggforce")
library(ggpubr); package.v[["ggpubr"]] <- packageVersion("ggpubr")
capture.output(package.v, file = "0_Package versions.txt")


# Load the data
list.files("../3_Burk_Comp_Norm/")
mix <- read.table( "../3_Burk_Comp_Norm/3_NormInten_without_control_MIX_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(mix)
# 6 23
gly <- read.table( "../3_Burk_Comp_Norm/3_NormInten_without_control_GLY_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(gly)
# 6 4
gln <- read.table( "../3_Burk_Comp_Norm/3_NormInten_without_control_GLN_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(gln)
# 6 7
asn <- read.table( "../3_Burk_Comp_Norm/3_NormInten_without_control_ASN_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(asn)
# 6 3
arg <- read.table( "../3_Burk_Comp_Norm/3_NormInten_without_control_ARG_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(arg)
# 6 3

# Set the treatments
myaminos <- c("mix","gly","gln", "asn", "arg")
heatmap_list <- list()

for (myaa in myaminos) {
  
  if(myaa == "gly") {
    mydata <- gly
    mytitle <- "Glycine"

  } else if( myaa == "gln") {
    mydata <- gln
    mytitle <- "Glutamine"

  } else if (myaa == "asn") {
    mydata <- asn
    mytitle <- "Asparagine"

  } else if (myaa == "arg") {
    mydata <- arg
    mytitle <- "Arginine"

  } else if (myaa == "mix") {
    mydata <- mix
    mytitle <- "Amino acid mixture"

  }
  
  # Get the compounds
  mycomps <- colnames(mydata[,-ncol(mydata)])
  
  stats_list <- list()
  groups_tbl <- data.frame()
  
  for(mycomp in mycomps) {
    
    #Select the compound
    mydata_comp <- mydata[,colnames(mydata) %in% c( mycomp, "Groups") ]
    
    # Run Student's t-test
    myform <- as.formula(paste(mycomp, " ~ Groups", sep = ""))
    myttest <- t.test(myform, data = mydata_comp)
    stats_list[[paste(mycomp)]][["Student's t-Test"]] <- myttest
    
    # Calculate the effect size
    ef.size <- cohen.d(mydata_comp[,1], mydata_comp$Groups)
    stats_list[[paste(mycomp)]][["Effect size: Cohen's d"]] <- ef.size
    
    # Prepare the data table
    my_groups <- plyr::ddply(mydata_comp, c("Groups"), summarise,
                             Means = mean(get(mycomp)),
                             sd = sd(get(mycomp)),
                             se = sd/sqrt(length(get(mycomp)))  )
    
    # Calculate the logFC 
    my_groups$logFC <-  my_groups$Means[2] - my_groups$Means[1] / my_groups$Means[1]
    my_groups$Compound <- mycomp
    my_groups$test <- myttest$method
    my_groups$p.value <- myttest$p.value
    my_groups$eff.size <- ef.size$estimate
    my_groups$CI.lower = ef.size$conf.int[[1]]
    my_groups$CI.upper = ef.size$conf.int[[2]]
    my_groups$label = paste("Effect size:", round(ef.size$estimate, 2) )
    
    groups_tbl <- rbind(groups_tbl, my_groups)
    
  }
  
  capture.output(stats_list, file = paste0("1_Stats_list",toupper(myaa),".txt"))
  write.table(groups_tbl, file = paste0("2_", toupper(myaa), "_Means,sd,se,pval.txt"), sep = "\t", quote = FALSE, col.names = NA)
  
  #--#--# Heatmaps with the normalized intensities 
  # (I) Prepare the table for the heatmap
  mydata_clean <- mydata[, -ncol(mydata)]
  mydata_clean_t <- t(mydata_clean)
  mydata_clean_t_01 <- decostand(mydata_clean_t,MARGIN = 1, method = "range")
  
  # (II) Prepare the metadata table (col_annotation)
  metadata <- data.frame( mydata$Groups )
  rownames(metadata) <- rownames(mydata)
  colnames(metadata) <- "Medium"
  
  # (III) Prepare the colors for the heatmap...
  mycolors <- viridisLite::viridis(9)[9:1]
  # ...and for the experimental set up and put it in a list
  mytreatcls <- brewer.pal(4, "Set1")
  myhtmpcls <- list(Medium = c( "0.1 TSBA + AA" = mytreatcls[1], "0.1 TSBA" = mytreatcls[2] ) )


  # Draw the heatmap 
  myhtmp <-   
    pheatmap(mydata_clean_t_01, 
             color = mycolors,
             main = mytitle,
             fontsize_row = 10,
             # cluster_cols = T, cutree_cols = mycuttreecols,
             cluster_cols = F, gaps_col = seq(0,6, by = 3),
             # cluster_rows = T, cutree_rows = mycuttreerows,
             # cluster_rows = F,
             annotation_col = metadata,
             # annotation_row = logFC.fdr_fin,
             annotation_colors = myhtmpcls,
             border_color = "grey70", cellwidth = 25, cellheight = 12
             ,show_rownames = F
    )
  # dev.off()
  heatmap_list[[paste(myaa)]][["Heatmap"]] <- myhtmp
  
  
  
  #### Prepare the effect size plot ###
  # (Ι) Prepare the table to plot
  tbl_to_plot <- groups_tbl[,c(6,8:12)]
  tbl_to_plot <- tbl_to_plot[-c(seq(2,nrow(tbl_to_plot), by = 2)) ,]
  
  # Prepare the point color
  tbl_to_plot$point.color.label <- ifelse(tbl_to_plot$p.value < 0.05, "< 0.05", "> 0.05" )
  allpointcols <- brewer.pal(4,"Dark2")[c(4,1)]
  if( any(tbl_to_plot$p.value < 0.05) ) { 
    mypointcols <-  allpointcols 
  } else { 
    mypointcols <-  allpointcols[2]
  }

  # Order the compounds according to the heatmap
  htmpcomps <- rownames(mydata_clean_t)[ myhtmp$tree_row$order]
  orderedcomps <- htmpcomps[length(htmpcomps) : 1]
  tbl_to_plot$Compound <- factor(tbl_to_plot$Compound, levels = orderedcomps)
  
  # Draw the effect size plot
  myeffsizeplot <-
    ggplot(tbl_to_plot, aes(x=Compound, y= eff.size, color = factor(point.color.label) )) +
    geom_point(size = 3) +
    scale_color_manual(name = "Significance",
                       values =  mypointcols  ) +
    coord_flip() +
    geom_errorbar(aes(ymin = CI.lower, ymax = CI.upper), width = 0.2, color = "black"  ) +
    labs(y = "Effect Size" ) +
    theme_minimal() + 
    theme(legend.position = "top",
          legend.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank()    ) +
    geom_vline(xintercept = seq(0.5,nrow(tbl_to_plot)+0.5, by = 1), color = "gray90") +
    geom_hline(yintercept = -1 , color = "forestgreen", linewidth = 1) +
    scale_x_discrete(position = "top")
  
  
  if (myaa == "mix") {
    myeffsizeplot <- myeffsizeplot +
      theme(plot.margin = margin(0.8, 2, 2.7, 3, "cm"))
  } else if(myaa == "gly") {
    myeffsizeplot <- myeffsizeplot +
      theme(plot.margin = margin(4.6, 2, 7, 3, "cm"))
  } else if (myaa == "gln") {
    myeffsizeplot <- myeffsizeplot +
      # theme(plot.margin = margin(3.7, 2, 6.5, 3, "cm"))
      theme(plot.margin = margin(3.8, 2, 6.6, 3, "cm"))
  } else {
    myeffsizeplot <- myeffsizeplot +
      theme(plot.margin = margin(4.6, 2, 7.4, 3, "cm"))
  }
  # print(myeffsizeplot)
  # dev.off()
  heatmap_list[[myaa]][["Eff.size plot"]] <- myeffsizeplot
  
  
  
}


pdf("Figure 6 - Heatmaps .pdf", height = 6, width = 11, onefile = T)
ggarrange(heatmap_list$mix$Heatmap[[4]], heatmap_list$mix$`Eff.size plot`,
          heatmap_list$gly$Heatmap[[4]], heatmap_list$gly$`Eff.size plot`,
          heatmap_list$gln$Heatmap[[4]], heatmap_list$gln$`Eff.size plot`,
          heatmap_list$asn$Heatmap[[4]], heatmap_list$asn$`Eff.size plot`,
          heatmap_list$arg$Heatmap[[4]], heatmap_list$arg$`Eff.size plot`,
          ncol = 2, nrow = 1, widths = c(1,2))

dev.off()

# Prepare the legend
figure_legend <- 
  c("Figure 6. Left: Heatmap of the pareto-scaled intensities ranged between 0 and 1 for each compound that was emitted by Burkholderia AD24 when cultivated in 0.1 TSBA + AA (red, left columns) and 0.1 TSBA (blue, right columns). Middle: Effect size graph of the difference in the intensities of each compound that was produced when Burkholderia AD24 was cultivated in the two different media. Points indicate the p value of the Welch Two Sample t-testperformed for each compound (pink: p < 0.05, green: p > 0.05. Green line represents effect size of -1 standard deviation which is portraying a threshold to strong effect. Left: Annotation and retention time of the compounds emitted by Burkholderia AD24")
capture.output(figure_legend, file = "4_Htmp_effsize_figure viridis ranged_legend.txt")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#      3rd   Data Transformation / Normalization   #-#-#-# ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#      1st   Data Transformation / Normalization   #-#-#-# ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Set the working directory
library(here)
setwd(paste0(here::here(), "/Volatilomic analysis/5_DataNorm_fr_multivariate") )

# Load the necessary packages
package.v <- list()
package.v[["R Version"]] <- R.version$version.string
package.v[["R Studio Version"]] <- RStudio.Version()$version
library(openxlsx); package.v[["openxlsx"]] <- packageVersion("openxlsx")
library(functional); package.v[["functional"]] <- packageVersion("functional") #Necessary for pretreat command
capture.output(package.v, file = "0_Package versions.txt")

# Load the data
list.files("../0_Data_preparation/")
mix <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "AAmixTransp", startRow = 1, colNames = T, rowNames = T)
dim(mix)
# 12 147
gly <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "GlyTransp", startRow = 2, colNames = T, rowNames = T)
dim(gly)
# 12 116
gln <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "GlnTransp", startRow = 2, colNames = T, rowNames = T)
dim(gln)
# 12 136
asn <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "AspTransp", startRow = 2, colNames = T, rowNames = T)
dim(asn)
# 12 100
arg <- read.xlsx("../0_Data_preparation/MZmine Tables.xlsx", sheet = "ArgTransp", startRow = 2, colNames = T, rowNames = T)
dim(arg)
# 12 99

# Set the amino acid treatments
myaminos <- c("mix","gly","gln", "asn", "arg")

for(myaa in myaminos) {
  
  if(myaa == "gly") {
    mydata <- gly
    mytitle <- "Glycine"
    
  } else if( myaa == "gln") {
    mydata <- gln
    mytitle <- "Glutamine"
    
  } else if (myaa == "asn") {
    mydata <- asn
    mytitle <- "Asparagine"
    
  } else if (myaa == "arg") {
    mydata <- arg
    mytitle <- "Arginine"
    
  } else if (myaa == "mix") {
    mydata <- mix
    mytitle <- "Amino acid mixture"
  }
  
  # Get the compounds table
  mydt <- mydata[,-1]
  
  # Get the metadata table
  metadata <- data.frame(mydata[,1])
  rownames(metadata) <- rownames(mydata)
  colnames(metadata) <- "Groups"
  
  # Data pretreatment methods
  # The normalization procedures are grouped into three categories. You can use one or combine them to achieve better results.
  
  # (1) Centering converts all the concentrations to fluctuations around zero instead of around the mean of the metabolite concentrations. Hereby, it adjusts for differences in the offset between high and low abundant metabolites. It is therefore used to focus on the fluctuating part of the data and leaves only the relevant variation (being the variation between the samples) for analysis. Centering is applied in combination with all the methods described below.
  
  # (2) Transformations are nonlinear conversions of the data. They are generally applied to correct for heteroscedasticity, to convert multiplicative relations into additive relations, and to make skewed distributions (more) symmetric.
  #           1) Log2 Transformation, 
  #           2) Square Root Transformation,
  #           3) Cubic Root Transformation
  
  # (3) Scaling methods are data pretreatment approaches that divide each variable by a factor, the scaling factor, which is different for each variable. They aim to adjust for the differences in fold differences between the different metabolites by converting the data into differences in concentration relative to the scaling factor..
  #   Scaling based on data dispersion:
  #           1) Autoscaling
  #           2) Pareto scaling
  #           3) Range scaling
  #           4) Vast scaling
  #   Scaling based on average value:
  #           1) Level scaling
  
  
  # pretreat command takes as variables
  # transf = "log2" / "log10 / "sqrt" / "none"
  # center = TRUE / FALSE
  # scale = "none" / auto" / "pareto"
  
  source("../pretreat command V2.1.R")
  pretreat_mydt <- pretreat(mydt, transf = "none", center = T, scale = "pareto")
  mytrasf <- pretreat_mydt$mytrasf
  transf_mydt <- pretreat_mydt$transf_mydt
  
  ### Visualize the result of data transformation - Calculate the means of each compound/column OR just plot the values (which makes more sense)
  valuesbefore <- unlist(mydt)
  valuesafter <- unlist(transf_mydt)
  
  # Prepare the plot
  pdf(paste0("Data_Norm_", mytitle,"_", mytrasf, ".pdf"), height = 4, width = 6)
  
  par(mfrow= c(1,2), mar = c(4,4,4,2), cex.main = 1.5)
  plot(density(valuesbefore),
       col = "blue", lwd = 2,
       main = paste0(mytitle," Before Transformation")
       , xlab = "")
  plot(density(valuesafter),
       col = "green", lwd = 2,
       main = paste0(mytitle," After Transformation")
       , xlab = "")
  dev.off()
  
  #  Merge the normalized values table with the metadata file...
  normtbl <- merge(transf_mydt, metadata, by = "row.names")
  rownames(normtbl) <- normtbl$Row.names 
  normtbl <- normtbl[, 2:ncol(normtbl)]
  
  # ...and save it
  write.table(normtbl, paste0("Values_Norm_",mytitle,"_",mytrasf,".txt"), col.names = NA, quote = FALSE, sep = "\t", row.names = T)
  
}


#Clean your environment and continue with the following step



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#       Multivariate Analysis / Sample Ordination   #-#-#-# ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

## Set the working directory
library(here)
setwd(paste0(here::here(), "/Volatilomic analysis/6_PLSDA") )


## load the necessary packages
package.v <- list()
package.v[["R Version"]] <- R.version$version.string
package.v[["R Studio Version"]] <- RStudio.Version()$version
library(mixOmics); package.v[["mixOmics"]] <- packageVersion("mixOmics")
library(RColorBrewer); package.v[["RColorBrewer"]] <- packageVersion("RColorBrewer")
library(ggforce); package.v[["ggforce"]] <- packageVersion("ggforce")
library(vegan); package.v[["vegan"]] <- packageVersion("vegan")
capture.output(package.v, file = "0_Package versions.txt")


# Load the data
list.files("../5_DataNorm_fr_multivariate/")
mix <- read.table( "../5_DataNorm_fr_multivariate/Values_Norm_Amino acid mixture_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(mix)
# 12 147
gly <- read.table( "../5_DataNorm_fr_multivariate/Values_Norm_Glycine_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(gly)
# 12 116
gln <- read.table( "../5_DataNorm_fr_multivariate/Values_Norm_Glutamine_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(gln)
# 12 136
asn <- read.table( "../5_DataNorm_fr_multivariate/Values_Norm_Asparagine_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(asn)
# 12 100
arg <- read.table( "../5_DataNorm_fr_multivariate/Values_Norm_Arginine_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t", check.names = T)
dim(arg)
# 12 99

# Set the treatments
myaminos <- c("mix","gly","gln", "asn", "arg")
plot_list <- list()
# ### Run Partial Least Square Discriminant Analysis (PLS-DA)
for(myaa in myaminos) {
  
  if(myaa == "gly") {
    mydata <- gly
    mytitle <- "Glycine"
    
  } else if( myaa == "gln") {
    mydata <- gln
    mytitle <- "Glutamine"

  } else if (myaa == "asn") {
    mydata <- asn
    mytitle <- "Asparagine"

  } else if (myaa == "arg") {
    mydata <- arg
    mytitle <- "Arginine"

  } else if (myaa == "mix") {
    mydata <- mix
    mytitle <- "Amino acid mixture"
  }
  
  # Get the compounds table
  mydt <- mydata[,-ncol(mydata)]
  
  # Get the metadata table
  metadata <- data.frame(mydata[,ncol(mydata)])
  rownames(metadata) <- rownames(mydata)
  colnames(metadata) <- "Groups"
  
  # Set the order you want the treatments to be presented in the graphs
  metadata$Groups <- factor(metadata$Groups, levels = levels(as.factor(metadata$Groups)))
  
  ##  run PLS-DA without scaling since our samples are already scaled. 
  set.seed(18)
  myplsda <- plsda(mydt, metadata$Groups, ncomp = 4, scale = F) # With mixOmics
  # plotIndiv(myplsda)
  
  # Focus on Compartments 1 and 2
  xaxis <- 1
  yaxis <- 2
  
  # Get the loading of the samples
  mysampleloadings <- data.frame(myplsda$variates$X[,c(xaxis,yaxis)])
  plot_tbl <- merge(metadata, mysampleloadings, by = "row.names")
  
  # Prepare the colors
  myplotcls <- brewer.pal(4, "Set1")
  # myplotcls <- myplotcls[unique(as.numeric(as.factor(metadata$Groups)))]
  
  # Set the main title 
  mymain <- paste0(mytitle)
  
  # Prepare the axes titles with the percentage of explained variance
  myexplvar <- myplsda$prop_expl_var$X
  # Prepare the axes titles with the percentage of explained variability
  xlab <- paste0("Component ", xaxis, ": ", round( 100 * myexplvar[xaxis], 2), "%")
  ylab <- paste0("Component ", yaxis, ": ", round( 100 * myexplvar[yaxis], 2), "%")
  
  #Prepare the plot
  myplsdaplot <-
    ggplot(plot_tbl, aes(x= comp1, y = comp2, color = Groups)) +
    geom_mark_ellipse(aes(x= comp1, y = comp2, group = Groups),color = "gray20", expand = 0, show.legend = F, inherit.aes = F) +
    geom_point(size = 3) + 
    scale_color_manual(values = myplotcls) +
    labs(title = mymain, x = xlab, y = ylab ) +
    theme_classic() +
    theme(title = element_text(face = "bold", size = 14),
          legend.position = "bottom",
          legend.title = element_blank())
  
  plot_list[[myaa]] <- myplsdaplot
  
}

pdf("Figure 5 - PLSDA plots.pdf", height = 9, width = 9)
ggarrange(plotlist = plot_list, ncol = 2, nrow = 3, common.legend = T, legend = "right")
dev.off()










































## load the necessary packages
package.v <- list()
package.v[["R Version"]] <- R.version$version.string
package.v[["R Studio Version"]] <- RStudio.Version()$version
library(vegan); package.v[["vegan"]] <- packageVersion("vegan")
library(mixOmics); package.v[["mixOmics"]] <- packageVersion("mixOmics")
library(Polychrome); package.v[["Polychrome"]] <- packageVersion("Polychrome")
library(RColorBrewer); package.v[["RColorBrewer"]] <- packageVersion("RColorBrewer")
capture.output(package.v, file = "0_Package versions.txt")
# Load the data
list.files("../4_2nd_data_norm//")

mydata <- read.table(file = "../4_2nd_data_norm/1st_Normalized_Values_No Transf,Centering,Paretoscaling.txt", header = T, row.names = 1, sep = "\t")
# str(mydata)
dim(mydata)
# 12 284
# 12 147

# Get the compounds table
mydt <- mydata[,-ncol(mydata)]

# Prepare the metadata file
metadata <- data.frame(mydata[,ncol(mydata)])
rownames(metadata) <- rownames(mydata)
colnames(metadata) <- "Groups"


##  run PLS-DA without scaling since our samples are already scaled. Study PCA scree plot and scatterplot in order to use a meaningful number of components to include in the model
myplsda <- plsda(mydt, metadata$Groups, ncomp = 4, scale = F) # With mixOmics
# plotIndiv(myplsda)

# Focus on Compartments 1 and 2
xaxis <- 1
yaxis <- 2

# Get the loading of the samples
mysampleloadings <- myplsda$variates$X[,c(xaxis,yaxis)]

# Prepare the colors
myplotcls <- createPalette(12,  brewer.pal(6, "Set1"))
myplotcls <- myplotcls[unique(as.numeric(as.factor(metadata$Groups)))]

# Set the main title 
mymain <- paste0("PLS-DA scores plot")

# Set the order you want the treatments to be presented in the graphs
# mytreatorder <- c("TSBA_CTL", "TSBA_BURK","ARG_CTL", "ARG_BURK", "ASP_CTL", "ASP_BURK", "GLU_CTL", "GLU_BURK", "GLY_CTL", "GLY_BURK")
# metadata$Groups <- factor(metadata$Groups, levels = mytreatorder)

# Prepare the axes titles with the percentage of explained variance
myexplvar <- myplsda$prop_expl_var$X
# Prepare the axes titles with the percentage of explained variability
xlab <- paste0("Component ", xaxis, ": ", round( 100 * myexplvar[xaxis], 2), "%")
ylab <- paste0("Component ", yaxis, ": ", round( 100 * myexplvar[yaxis], 2), "%")

# Prepare & save the plot
pdf(paste0("PLSDA_plot_good.pdf"), width = 5, height = 4 )
par(bty = "n", mar = c(5,5,3,2), xpd =T)

par(adj = 0.5)
plot(mysampleloadings[,1], mysampleloadings[,2], type = "n", frame.plot = F, 
     xlab = xlab, ylab = ylab, 
     cex.axis = 1.1, cex.lab = 1.1, cex.main = 2)
# title(main = mymain, line = 3, cex.main = 2, adj=0)
mtext(mymain, side = 3, line = 1, adj = 0.5, cex = 1.5, font = 2)
box(col = "gray80", which  = "figure")

# (I) First add the ellipses in order to keep them in the background and prevent them to obstruct information
par(xpd = NA)
ordiellipse(mysampleloadings[,1:2], 
            groups = metadata$Groups, 
            col = myplotcls[1:length(levels(as.factor(metadata$Groups)))], 
            label = F, #label = T,
            kind = "ehull",
            lty = 1, lwd=1)

# (II) Add the sample points and give them colors one for each treatment
points(mysampleloadings[,1], mysampleloadings[,2], 
       bg = myplotcls[as.numeric(as.factor(metadata$Groups))], 
       pch = 21, cex = 1.5, col = "grey20")
# Add the text to the points if needed
# text(mysampleloadings[,1], mysampleloadings[,2], labels = rownames(mysampleloadings), cex = 0.6, adj = c(-0.1, -0.4), srt = 10, xpd = T)

# (III) Add the legend
par(xpd = T)
legend("topleft", inset = c(-0.23,-0.22), ncol = 1,
       legend = levels(as.factor(metadata$Groups)), cex = 0.6,
       # x.intersp = 0.5, 
       # text.width = mean(strwidth(levels(as.factor(metadata$Groups)), units = "user", cex = 1.2))+0.2,
       pt.bg = myplotcls[1:length(levels(as.factor(metadata$Groups)))], pch = 21, pt.cex = 0.8, 
       bty = "o", box.col = "gray80")

dev.off()

























