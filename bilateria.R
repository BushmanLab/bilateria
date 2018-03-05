# Libraries----

library(pander)
library(ape)
library(phytools)
library(phylobase)
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(dnar)
library(Rtsne)
library(lattice)
library(grid)
library(gridExtra)
library(ggimage)
library(phyloseq)
library(eclectic)
library(ade4)

# Setup----

work_dir <- "/home/kevin/projects/bilateria/"

mapping_file_fp <- file.path(work_dir, "map/map.tsv")
otu_table_fp <- file.path(work_dir, "otu/otu_table.txt")
uu_fp <- file.path(work_dir, "beta_diversity/unweighted_unifrac_otu_table.original.txt")
wu_fp <- file.path(work_dir, "beta_diversity/weighted_unifrac_otu_table.original.txt")

# Functions----

E <- new.env()
data_dir=paste0(work_dir, "validation/anteater/")

resetValidationData <- function(E, data_dir, filterbyweight=T) {
  library(qiimer)
  E$otuTableFile <- paste0(data_dir, "otu/otu_table.txt")
  E$mapFile <- paste0(data_dir, "map/map.tsv")
  
  E$o <- read_qiime_otu_table(E$otuTableFile)
  E$s <- read_qiime_mapping_file(E$mapFile)
  E$s <- E$s[E$s$toKeep=="keep",]
  if(filterbyweight) {
    E$s <- E$s[!is.na(E$s$weight),]
  }

  alignQiimeData(E)
  E$o$metadata <- E$o$metadata[!grepl("*Chloroplast*",E$o$metadata)]
  E$o$counts <- E$o$counts[names(E$o$metadata),]
  E$o$otu_ids <- names(E$o$metadata)
  
  E$sampleOrdering <- with(E$s,order(ave(1:nrow(E$s),E$s$phylum,FUN=length),phylum,
                                         ave(1:nrow(E$s),E$s$class,FUN=length),class,
                                         ave(1:nrow(E$s),E$s$order,FUN=length),order,
                                         ave(1:nrow(E$s),E$s$family,FUN=length),family,
                                         ave(1:nrow(E$s),E$s$genus,FUN=length),genus,
                                         ave(1:nrow(E$s),E$s$species,FUN=length),species,
                                         ave(1:nrow(E$s),E$s$common,FUN=length),common))
  
  E$speciesTax <- unique(E$s[,c("phylum","class","order","family","genus","common")])
  E$speciesTax <- unique(E$speciesTax)
  E$speciesOrdering <- with(E$speciesTax,
                              order(ave(1:nrow(E$speciesTax),E$speciesTax$phylum,FUN=length),phylum,
                                    ave(1:nrow(E$speciesTax),E$speciesTax$class,FUN=length),class,
                                    ave(1:nrow(E$speciesTax),E$speciesTax$order,FUN=length),order,
                                    ave(1:nrow(E$speciesTax),E$speciesTax$family,FUN=length),family,
                                    ave(1:nrow(E$speciesTax),E$speciesTax$genus,FUN=length),genus))
  
  #  E$figure1SpeciesOrdering <- match(E$,E$speciesTax$common)
  E$s$weight <- as.numeric(E$s$weight)
  E$s$log.weight <- log10(E$s$weight)
  E$s$numOtus <- colSums(E$o$counts>0)
  E$s$filteredReadCount <- colSums(E$o$counts)
  
  md <- sub("(; [kpcofgs]__)+$", "", E$o$metadata, perl=T)
  library(tibble)
  
  adf <- split_assignments(md) %>% 
    mutate_all(funs(sub("[kpcofgs]__", "", .))) %>%
    mutate_all(funs(ifelse(. %in% "", NA, .))) %>%
    mutate(Phylum = ifelse(Genus %in% "Enterococcus", "Firmicutes", Phylum)) %>%
    mutate(Class = ifelse(Genus %in% "Enterococcus", "Bacilli", Class)) %>%
    mutate(Order = ifelse(Genus %in% "Enterococcus", "Lactobacillales", Order)) %>%
    mutate(Family = ifelse(Genus %in% "Enterococcus", "Enterococcaceae", Family))
  E$a <- simplify_assignments(adf)
  names(E$a) <- rownames(E$o$counts)
  
  adf <- adf %>% rownames_to_column("OtuID") %>%
    mutate(OtuID = rownames(E$o$counts))
  
  library(dplyr)
  library(readr)
  aero <- read_tsv(paste0(work_dir, "oxygen/bilaterian_aerobic_status_data.tsv")) %>%
    mutate(TaxonRank = factor(TaxonRank, levels=rev(taxonomic_ranks))) %>%
    arrange(TaxonRank, TaxonName)
  
  E$adf_aero <- adf %>%
    mutate(aerobic_status = NA, explanation = NA)
  for (row_idx in seq_len(nrow(aero))) {
    taxon_rank <- as.character(aero$TaxonRank[row_idx])
    taxon_name <- aero$TaxonName[row_idx]
    taxon_status <- aero$AerobicStatus[row_idx]
    explanation <- aero$Explanation[row_idx]
    is_taxon <- E$adf_aero[[taxon_rank]] %in% taxon_name
    is_unlabeled <- is.na(E$adf_aero$aerobic_status)
    E$adf_aero$aerobic_status[is_taxon & is_unlabeled] <- taxon_status
    E$adf_aero$explanation[is_taxon & is_unlabeled] <- explanation
  }
  
  E$cts_df <- E$o$counts %>%
    reshape2:::melt.array(
      varnames = c("OtuID", "SampleID"), 
      value.name = "Counts",
      as.is=TRUE) %>%
    filter(Counts > 0) %>%
    group_by(SampleID) %>%
    mutate(Proportion = Counts / sum(Counts)) %>%
    ungroup() %>%
    left_join(E$adf_aero, by="OtuID") %>%
    mutate(., Assignment = simplify_assignments(.[,taxonomic_ranks])) %>%
    left_join(E$s[,c("SampleID","common")], by="SampleID")
  
  E$samples_aero <- E$cts_df %>%
    filter(!aerobic_status %in% c("NA", "no taxonomic assignment", "Streptophyta")) %>%
    group_by(SampleID, aerobic_status, common) %>%
    summarize(Proportion = sum(Proportion), Counts = sum(Counts))

  E$s <- E$s %>% left_join(E$samples_aero %>%
                             filter(aerobic_status=="obligate anaerobe") %>%
                             ungroup() %>%
                             select(SampleID, Proportion), by="SampleID") %>%
    rename(obligateAnaerobeProp = Proportion) %>%
    mutate_at(.vars=vars(obligateAnaerobeProp),
              funs(replace(obligateAnaerobeProp, which(is.na(.)), 0)))
    
  # Load phylogenetic diversity
  E$pdwt_fp <- file.path(data_dir, "alpha_diversity.tsv")
  E$pdwt <- read.table(E$pdwt_fp)
  E$s$phyDiv <- E$pdwt[E$s$SampleID,"PD_whole_tree"]
}

resetData <- function(status = "all") {
  library(qiimer)
  library(phylobase)
  library(adephylo)
  library(ape)
  
  G <<- new.env()
  G$s <- read_qiime_mapping_file(mapping_file_fp)
  if(status=="toUse") {
    G$s <- G$s[G$s$toKeep=="keep" & !is.na(G$s$Weight_to_use),]
  }
  
  if(status=="toUsePlusControls") {
    G$s <- G$s[G$s$toKeep=="keep" | G$s$sample_type %in% c("GeneBlock","Extraction_blank",
                                                           "Dissection_blank", "PCR_H2O"),]
  }

  # Code to read in om.csv and create G$o object  
  # om <- read.csv(paste0(work_dir,"om.csv"))
  # o <- dcast(data = om,formula = OTU~SampleID,fun.aggregate = sum,value.var = "count")
  # rownames(o) <- o$OTU
  # G$o$counts <- as.matrix(o)
  # uniq_metadata <- unique(om[,c("OTU","metadata")])
  # G$o$sample_ids <- colnames(G$o$counts)
  # G$o$otu_ids <- rownames(G$o$counts)
  # G$o$metadata <- as.character(uniq_metadata$metadata)
  # names(G$o$metadata) <- as.character(uniq_metadata$OTU)
  # rm(uniq_metadata)
  # rm(o)
  # rm(om)
  
  # Use this code to create OTU table to save time on future loads
  # write.table(cbind(G$o$counts,G$o$metadata),
  #            paste0(work_dir,"otu/otu_table.txt"), quote=FALSE, sep = "\t")
  # system(paste0("sed -i '1s/^/# Constructed from biom file\\n#OTU ID\t/' ",work_dir,"otu/otu_table.txt"))
  # system(paste0("sed -i '2s/$/Consensus Lineage/' ",work_dir,"otu/otu_table.txt"))
  
  # Use this line in place of the above code block defining G$o if you have already created OTU table
  G$o <- read_qiime_otu_table(otu_table_fp)
  
  alignQiimeData(G)
  rownames(G$s) <- G$s$SampleID
  colnames(G$o$counts) <- G$s$SampleID
  G$o$sample_ids <- G$s$SampleID
  if(!status=="all") {
    G$o$metadata <- G$o$metadata[!grepl("*Chloroplast*",G$o$metadata)]
  }
  G$o$counts <- G$o$counts[names(G$o$metadata),]
  G$o$otu_ids <- names(G$o$metadata)
 
  # Code used to create om.csv (om = OTU melted)
  #om <- melt(G$o$counts)
  #om <- om[om$value>0,]
  #om$metadata <- G$o$metadata[om$Var1]
  #colnames(om) <- c("OTU","SampleID","count","metadata")
  #write.csv(om, paste0(work_dir,"om.csv"), row.names = F)
  
  G$sampleOrdering <- with(G$s,order(ave(1:nrow(G$s),G$s$phylum,FUN=length),phylum,
                                     ave(1:nrow(G$s),G$s$supersuperclass,FUN=length),supersuperclass,
                                     ave(1:nrow(G$s),G$s$superclass,FUN=length),superclass,
                                     ave(1:nrow(G$s),G$s$class,FUN=length),class,
                                     ave(1:nrow(G$s),G$s$superorder,FUN=length),superorder,
                                     ave(1:nrow(G$s),G$s$clade,FUN=length),clade,
                                     ave(1:nrow(G$s),G$s$order,FUN=length),order,
                                     ave(1:nrow(G$s),G$s$family,FUN=length),family,
                                     ave(1:nrow(G$s),G$s$genus,FUN=length),genus,
                                     ave(1:nrow(G$s),G$s$species,FUN=length),species,
                                     ave(1:nrow(G$s),G$s$common,FUN=length),common))
  
  G$speciesTax <- unique(G$s[,c("phylum","supersuperclass","superclass","class","superorder","clade",
                                "order","family","genus","common","timetree_proxy")])
  G$speciesTax$genus <- as.character(G$speciesTax$genus)
  G$speciesTax$genus[G$speciesTax$common=="Bee"] <- "Bee genus"
  G$speciesTax <- unique(G$speciesTax)
  G$speciesOrdering <- with(G$speciesTax,
                            order(ave(1:nrow(G$speciesTax),G$speciesTax$phylum,FUN=length),phylum,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$supersuperclass,FUN=length),supersuperclass,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$superclass,FUN=length),superclass,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$class,FUN=length),class,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$superorder,FUN=length),superorder,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$clade,FUN=length),clade,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$order,FUN=length),order,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$family,FUN=length),family,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$genus,FUN=length),genus))

  G$s$speciesOrder <- match(G$s$common,G$speciesTax$common[G$figure1SpeciesOrdering])
    
  md <- sub("(; [kpcofgs]__)+$", "", G$o$metadata, perl=T)
  adf <- split_assignments(md)
  G$a <- simplify_assignments(adf)

  # Load alpha diversity, calculated using PD whole tree (PWT)
  G$pdwt_fp <- file.path(work_dir, "alpha_diversity.tsv")
  G$pdwt <- read.table(G$pdwt_fp)
  G$s$phyDiv <- G$pdwt[G$s$SampleID,"PD_whole_tree"]
  
  library(ape)  
  G$newick <- read.tree(paste0(work_dir, "speciesPhylogeny.nwk"))

  G$timetreeLabels <-
    tolower(gsub("_"," ",
                 G$newick$tip.label[G$newick$edge[,2][G$newick$edge[,2]<=length(G$newick$tip.label)]]))
  
  G$figure1SpeciesOrdering <- match(G$timetreeLabels,G$speciesTax$timetree_proxy)
  G$s$Weight_to_use <- as.numeric(G$s$Weight_to_use)
  G$s$log.weight <- log10(G$s$Weight_to_use)
  
  G$species <- G$s %>%
    select(common,phylum,class,order,family,diet,specificDiet,timetree_proxy) %>%
    filter(common %in% G$speciesTax$common) %>%
    slice(match(G$speciesTax$common[G$figure1SpeciesOrdering],common)) %>%
    as.data.frame
  rownames(G$species) <- G$species$common
}

resetUnifracData <- function(){
  library(qiimer)
  G$uur_fp <- file.path(work_dir, "bd_uur", "unweighted_unifrac_otu_rare.txt")
  G$uur <- read_qiime_distmat(G$uur_fp)
  G$uur <- dist_subset(G$uur, G$s$SampleID)

  G$wnu_fp <- file.path(work_dir, "bd_wnu", "weighted_normalized_unifrac_otu_table.txt")
  G$wnu <- read_qiime_distmat(G$wnu_fp)
  G$wnu <- dist_subset(G$wnu, G$s$SampleID)
}

alignQiimeData <- function(E) {
  o <- E$o
  s <- E$s
  if(!is.list(o)) stop("Object o is not of type list")
  if(!is.data.frame(s)) stop("Object s is not of type data frame")
  if(!"sample_ids" %in% names(o)) stop("Object o needs element named 'sample_ids'")
  if(!"otu_ids" %in% names(o)) stop("Object o needs element named 'otu_ids'")
  if(!"metadata" %in% names(o)) stop("Object o needs element named 'metadata'")
  if(!"counts" %in% names(o)) stop("Object o needs element named 'counts'")
  if(!"SampleID" %in% colnames(s)) stop("Object s needs column named 'SampleID'")
  if(is.null(colnames(o$counts))) stop("o$counts needs named columns")
  if(is.null(rownames(o$counts))) stop("o$counts needs named rows")
  if(is.null(names(o$metadata))) stop("o$metadata must be a named vector")
  
  s <- s[s$SampleID %in% colnames(o$counts),]
  s <- droplevels(s)
  
  o$counts <- o$counts[,colnames(o$counts) %in% s$SampleID]
  o$counts <- o$counts[rowSums(o$counts)>0,]
  o$counts <- o$counts[,match(s$SampleID, colnames(o$counts))]
  o$sample_ids <- colnames(o$counts)
  o$metadata <- o$metadata[rownames(o$counts)]
  o$otu_ids <- o$otu_ids[o$otu_ids %in% rownames(o$counts)]
  E$s <- s
  E$o <- o
}

expandTree <- function(tree,replacingLabel,newLabel,newEdgeLength){
  replacingTip <- which(tree$tip.label==replacingLabel)
  replacingEdge <- which(tree$edge[,2]==replacingTip)
  
  numNodes <- max(tree$edge)
  newNode <- numNodes + 2
  newTip <- replacingTip + 1
  numTips <- length(tree$tip.label)
  
  tree$edge[tree$edge > replacingTip] <- tree$edge[tree$edge > replacingTip] + 1
  tree$edge[tree$edge==replacingTip] <- newNode
  tree$edge <- cbind(append(append(tree$edge[,1],newNode,replacingEdge),newNode,replacingEdge),
                     append(append(tree$edge[,2],newTip,replacingEdge),replacingTip,replacingEdge))
  
  tree$edge.length[replacingEdge] <- tree$edge.length[replacingEdge] - newEdgeLength
  tree$edge.length <- append(append(tree$edge.length,newEdgeLength,replacingEdge),newEdgeLength,replacingEdge)
  
  tree$tip.label <- append(tree$tip.label,newLabel,replacingTip)
  tree$Nnode <- tree$Nnode + 1
  tree$node.label <- c(tree$node.label, paste0("'",max(as.numeric(gsub("'","",tree$node.label[-1])))+1,"'"))
  
  return(tree)
}

expandTreeInternal <- function(tree,replacingNode,newLabel,newTipEdgeLength){
  replacingEdge <- which(tree$edge[,2]==replacingNode)
  replacingNodeLeaves <- as(subset(as(tree,"phylo4"),node.subtree=replacingNode),"phylo")$tip.label
  newTip <- max(match(replacingNodeLeaves,tree$tip.label)) + 1
  
  numNodes <- max(tree$edge)
  newNode <- numNodes + 2
  numTips <- length(tree$tip.label)
  replacingNodeHeight <- max(node.depth.edgelength(as(subset(as(tree,"phylo4"),
                                                             node.subtree=replacingNode),"phylo")))
  newInternalEdgeLength <- newTipEdgeLength - replacingNodeHeight
  newReplacingEdgeLength <- tree$edge.length[replacingEdge] - newInternalEdgeLength
  
  
  tree$edge[tree$edge >= newTip] <- tree$edge[tree$edge >= newTip] + 1
  replacingNode <- replacingNode + 1
  tree$edge[tree$edge[,1]==replacingNode,1] <- newNode
  tree$edge <- cbind(append(append(tree$edge[,1],replacingNode,replacingEdge),replacingNode,replacingEdge),
                     append(append(tree$edge[,2],newNode,replacingEdge),newTip,replacingEdge))
  
  tree$edge.length <- append(append(tree$edge.length,newInternalEdgeLength,replacingEdge),
                             newTipEdgeLength,replacingEdge)
  tree$edge.length[replacingEdge] <- newReplacingEdgeLength
  
  tree$tip.label <- append(tree$tip.label,newLabel,newTip - 1)
  tree$Nnode <- tree$Nnode + 1
  tree$node.label <- c(tree$node.label, paste0("'",max(as.numeric(gsub("'","",tree$node.label[-1])))+1,"'"))
  
  return(tree)
}

shannonDiv <- function(data){
  H <- numeric(ncol(data))
  mult <- numeric(nrow(data))
  for(i in 1:ncol(data)){
    prop <- data[,i] / sum(data[,i])
    for(j in 1:nrow(data)){
      mult[j] <- prop[j] * log(prop[j])
    }
    H[i] <-- sum(mult, na.rm=T)
  }
  return(cbind(colnames(data),H))
}

plotFigure1 <- function(E, study) {
  
  library(reshape2)
  library(dplyr)
  library(ggplot2)
  library(dnar)
  library(grid)
  
  E$aPhylum <- sub("; c__(.)*$", "", E$o$metadata, perl=T)
  E$oByPhylum <- rowsum(E$o$counts,E$aPhylum)
  
  E$obpProp <- apply(E$oByPhylum,2,function(x) ifelse (x,x / sum(x),0))
  E$obpProp <- data.frame(t(E$obpProp[rowSums(E$obpProp)>0,]))
  
  E$obpProp$common <- as.character(E$s$common[match(rownames(E$obpProp),E$s$SampleID)])
  E$obpPropBySpecies <- aggregate(.~common,data=E$obpProp,mean)
  
  E$otuPhylaToUse <-
    as.character(read.table(paste0(work_dir,"otu_phyla.txt"),header = F,stringsAsFactors = F)$V1)
  
  E$obpPropBySpeciesToUse <- E$obpPropBySpecies[,names(E$obpPropBySpecies) == 'common' |
                                                  names(E$obpPropBySpecies) %in% E$otuPhylaToUse]
  
  colnames(E$obpPropBySpeciesToUse)[-1] <-
    colnames(E$obpPropBySpeciesToUse)[rev(order(colSums(E$obpPropBySpeciesToUse[-1])))+1]
  E$obpPropBySpeciesToUse[,-1] <-
    E$obpPropBySpeciesToUse[,rev(order(colSums(E$obpPropBySpeciesToUse[-1])))+1]
  
  E$obpPropBySpeciesToUse$Other <- 1 - rowSums(E$obpPropBySpeciesToUse[-1])
  
  E$obpHeatmapData <- melt(cbind(E$obpPropBySpeciesToUse), id.vars=c('common'))
  
  plotHeatmap <- ggplot(E$obpHeatmapData, aes(x=variable, y=common, fill=value)) +
    scale_y_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    scale_x_discrete(limits=c(G$otuPhylaToUse,"Other")) +
    geom_tile(color="grey80", size=0.4) + theme_grey() +
    theme(
      legend.position = "none",
      strip.text.y = element_text(angle=0, vjust=0),
      strip.text.x = element_text(angle=90, vjust=0),
      panel.border = element_blank(),
      axis.text.y = element_text(size=12),
      axis.text.x = element_text(angle = 45, hjust = 1, size=10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()) +
    eclectic::saturated_rainbow_pct()
  
  
  # Bar chart of max OTU
  
  E$o$prop <- apply(E$o$counts, 2, function(x) ifelse (x, x / sum(x),0))
  
  E$o$maxProp <- apply(E$o$prop, 2, max)
  E$o$maxProp <- E$o$maxProp[match(E$s$SampleID[E$sampleOrdering],
                                   names(E$o$maxProp))]
  E$maxPropBySpecies <-
    tapply(E$o$maxProp[match(E$s$SampleID, names(E$o$maxProp))],
           E$s$common,mean)
  E$maxPropBySpecies <-
    E$maxPropBySpecies[match(E$speciesTax$common[E$speciesOrdering],
                             names(E$maxPropBySpecies))]
  
  E$mpbsDf <- data.frame(common=names(E$maxPropBySpecies),value=unname(E$maxPropBySpecies))
  E$mpDf <- data.frame(common=E$s$common[match(names(E$o$maxProp),E$s$SampleID)],
                       value=E$o$maxProp)
  plotMax <- ggplot(E$mpbsDf,aes(x=common,y=value,fill="cyan")) +
    geom_bar(stat="identity",color="cyan") +
    scale_fill_manual(values=c("cyan")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="[title]") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$mpDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)
  
  
  # Bar chart of shannon diversity
  
  E$otuShannon <- shannonDiv(E$o$counts)
  E$otuShannonDf <- data.frame(SampleID=as.character(E$otuShannon[,1]),
                               shannon=as.numeric(E$otuShannon[,2]))
  E$otuShannonDf$common <- E$s$common[match(E$s$SampleID,E$otuShannonDf$SampleID)]
  E$otuShannonBySpecies <- aggregate(shannon ~ common, data=E$otuShannonDf, mean)
  E$otuShannonBySpecies <-
    E$otuShannonBySpecies[match(E$speciesTax$common[E$speciesOrdering],
                                E$otuShannonBySpecies$common),]
  
  plotShannon <- ggplot(E$otuShannonBySpecies,aes(x=common,y=shannon, fill="purple")) +
    geom_bar(stat="identity",color="purple") +
    scale_fill_manual(values=c("purple")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="[title]") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$otuShannonDf,aes(x=common,y=shannon),color="black",shape="+",size=4,show.legend = FALSE)
  
  
  # Bar chart of percent reads unassigned
  
  E$unassignedOtuIds <- names(E$a)[E$a=="Unassigned"]
  E$percentTotalUnassignedBySample <-
    apply(E$o$counts,2,function(X) sum(X[names(X) %in% E$unassignedOtuIds]) / sum(X))
  E$percentTotalUnassignedBySample <-
    E$percentTotalUnassignedBySample[match(E$s$SampleID[E$sampleOrdering],
                                           names(E$percentTotalUnassignedBySample))]
  E$s$percentTotalUnassigned <-
    E$percentTotalUnassignedBySample[match(E$s$SampleID, names(E$percentTotalUnassignedBySample))]
  
  
  E$percentTotalUnassignedBySpecies <-
    tapply(E$percentTotalUnassignedBySample[match(E$s$SampleID, names(E$percentTotalUnassignedBySample))],
           E$s$common,mean)
  E$percentTotalUnassignedBySpecies <-
    E$percentTotalUnassignedBySpecies[match(E$speciesTax$common[E$speciesOrdering],
                                            names(E$percentTotalUnassignedBySpecies))]
  
  E$ptubsDf <- data.frame(common=names(E$percentTotalUnassignedBySpecies),
                          value=unname(E$percentTotalUnassignedBySpecies))
  E$ptuDf <- data.frame(common=E$s$common[match(names(E$percentTotalUnassignedBySample),E$s$SampleID)],
                        value=E$percentTotalUnassignedBySample)
  
  plotUnassigned <- ggplot(E$ptubsDf,aes(x=common,y=value,fill="red")) +
    geom_bar(stat="identity",color="red") +
    scale_fill_manual(values=c("red")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="[title]") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$ptuDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)
  
  
  
  # Barchart of Wolbachia OTUs
  
  E$wolbachiaOtus <- names(E$a[grep("Wolbachia",E$a)])
  E$percentWolbachiaBySample <-
    apply(E$o$counts,2,function(X) sum(X[names(X) %in% E$wolbachiaOtus]) / sum(X))
  E$percentWolbachiaBySample <-
    E$percentWolbachiaBySample[match(E$s$SampleID[E$sampleOrdering],
                                     names(E$percentWolbachiaBySample))]
  
  E$percentWolbachiaBySpecies <-
    tapply(E$percentWolbachiaBySample[match(E$s$SampleID, names(E$percentWolbachiaBySample))],
           E$s$common,mean)
  E$percentWolbachiaBySpecies <-
    E$percentWolbachiaBySpecies[match(E$speciesTax$common[E$speciesOrdering],
                                      names(E$percentWolbachiaBySpecies))]
  
  E$pwbsDf <- data.frame(common=names(E$percentWolbachiaBySpecies),
                         value=unname(E$percentWolbachiaBySpecies))
  E$pwDf <- data.frame(common=E$s$common[match(names(E$percentWolbachiaBySample),E$s$SampleID)],
                       value=E$percentWolbachiaBySample)
  plotWolbachia <- ggplot(E$pwbsDf,aes(x=common,y=value, fill="orange")) +
    geom_bar(stat="identity",color="orange") +
    scale_fill_manual(values=c("orange")) + guides(fill=FALSE) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="[title]") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$pwDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


  # Barchart of phylogenetic diversity
  
  E$phyDivBySpecies <- tapply(E$s$phyDiv, E$s$common,mean)
  E$phyDivBySpecies <-
    E$phyDivBySpecies[match(E$speciesTax$common[E$speciesOrdering],
                            names(E$phyDivBySpecies))]
  
  E$pdbsDf <- data.frame(common=names(E$phyDivBySpecies),
                         value=unname(E$phyDivBySpecies))
  E$pdDf <- data.frame(common=E$s$common,
                       value=E$s$phyDiv)
  
  plotPhyDiv <- ggplot(E$pdbsDf,aes(x=common,y=value,fill="gray")) +
    geom_bar(stat="identity",color="gray") +
    scale_fill_manual(values=c("gray")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="[title]") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$pdDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)
  

  # Bar chart of obligate anaerobes
  
  E$obligateAnaerobeSpecies <- tapply(E$s$obligateAnaerobeProp,E$s$common,mean)
  E$obligateAnaerobeSpecies <-
    E$obligateAnaerobeSpecies[match(E$speciesTax$common[E$speciesOrdering],
                                    names(E$obligateAnaerobeSpecies))]
  
  E$oapbsDf <- data.frame(common=names(E$obligateAnaerobeSpecies),
                          value=unname(E$obligateAnaerobeSpecies))
  E$oapDf <- data.frame(common=E$s$common,
                        value=E$s$obligateAnaerobeProp)
  
  plotObligateAnaerobes <- ggplot(E$oapbsDf,aes(x=common,y=value,fill="pink")) +
    geom_bar(stat="identity",color="pink") +
    scale_fill_manual(values=c("pink")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="[title]") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$oapDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)
  
  
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

  library(grid)
  pdf(paste0(paste0(work_dir, "R_plots/figure1.",study,".pdf")),
      height=15, width=20, onefile = FALSE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1000, 14)))
  print(plotHeatmap, vp = vplayout(35:1000, 5:8))
  print(plotMax, vp = vplayout(33:921, 9))
  print(plotShannon, vp = vplayout(33:921, 10))
  print(plotPhyDiv, vp = vplayout(33:921, 11))
  print(plotObligateAnaerobes, vp = vplayout(33:921, 12))
  print(plotWolbachia, vp = vplayout(33:921, 13))
  print(plotUnassigned, vp = vplayout(33:921, 14))
  
  dev.off()
  
}

plotTsneLegends <- function() {
  legend(
    title="Class",
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4], #y pos in middle of plot
    names(G$classColors), col=G$classColors, pch=15, # labels and annotations here
    xjust=0,xpd=NA #left justify and plot outside region
  )
  
  G$mammalOrderShapeLegend <-
    G$classOrderNumber[names(G$classOrderNumber) %in% G$s$order[G$s$class=='Mammalia']]
  G$mammalOrderShapeLegend <- G$mammalOrderShapeLegend[order(G$mammalOrderShapeLegend)]
  legend(
    title="Mammalia",
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4]-0.4*diff(par('usr')[3:4]), #y pos in middle of plot
    names(G$mammalOrderShapeLegend), pch=as.numeric(G$mammalOrderShapeLegend),
    col=G$classColors['Mammalia'], pt.lwd=2,
    xjust=0,xpd=NA #left justify and plot outside region
  )
  
  G$insectOrderShapeLegend <-
    G$classOrderNumber[names(G$classOrderNumber) %in% G$s$order[G$s$class=='Insecta']]
  G$insectOrderShapeLegend <- G$insectOrderShapeLegend[order(G$insectOrderShapeLegend)]
  legend(
    title="Insecta",
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4]-0.7*diff(par('usr')[3:4]), #y pos in middle of plot
    names(G$insectOrderShapeLegend), pch=as.numeric(G$insectOrderShapeLegend),
    col=G$classColors['Insecta'], pt.lwd=2,
    xjust=0,xpd=NA #left justify and plot outside region
  )
}

plotExtraTsneLegends <- function() {
  G$actinOrderShapeLegend <-
    G$classOrderNumber[names(G$classOrderNumber) %in% G$s$order[G$s$class=='Actinopteri']]
  G$actinOrderShapeLegend <- G$actinOrderShapeLegend[order(G$actinOrderShapeLegend)]
  legend(
    title="Actinopteri",
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4]-0*diff(par('usr')[3:4]), #y pos in middle of plot
    names(G$actinOrderShapeLegend), pch=as.numeric(G$actinOrderShapeLegend),
    col=G$classColors['Actinopteri'],
    xjust=0,xpd=NA #left justify and plot outside region
  )
  
  G$avesOrderShapeLegend <-
    G$classOrderNumber[names(G$classOrderNumber) %in% G$s$order[G$s$class=='Aves']]
  G$avesOrderShapeLegend <- G$avesOrderShapeLegend[order(G$avesOrderShapeLegend)]
  legend(
    title="Aves",
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4]-0.25*diff(par('usr')[3:4]), #y pos in middle of plot
    names(G$avesOrderShapeLegend), pch=as.numeric(G$avesOrderShapeLegend),
    col=G$classColors['Aves'],
    xjust=0,xpd=NA #left justify and plot outside region
  )
  
  G$chonOrderShapeLegend <-
    G$classOrderNumber[names(G$classOrderNumber) %in% G$s$order[G$s$class=='Chondrichthyes']]
  G$chonOrderShapeLegend <- G$chonOrderShapeLegend[order(G$chonOrderShapeLegend)]
  legend(
    title="Chondrichthyes",
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4]-0.5*diff(par('usr')[3:4]), #y pos in middle of plot
    names(G$chonOrderShapeLegend), pch=as.numeric(G$chonOrderShapeLegend),
    col=G$classColors['Chondrichthyes'],
    xjust=0,xpd=NA #left justify and plot outside region
  )
  
  G$malaOrderShapeLegend <-
    G$classOrderNumber[names(G$classOrderNumber) %in% G$s$order[G$s$class=='Malacostraca']]
  G$malaOrderShapeLegend <- G$malaOrderShapeLegend[order(G$malaOrderShapeLegend)]
  legend(
    title="Malacostraca",
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4]-0.75*diff(par('usr')[3:4]), #y pos in middle of plot
    names(G$malaOrderShapeLegend), pch=as.numeric(G$malaOrderShapeLegend),
    col=G$classColors['Malacostraca'],
    xjust=0,xpd=NA #left justify and plot outside region
  )
}

countsToProp <- function(counts) {
  return(apply(counts, 2, function(x)   ifelse (x, x / sum(x),0)))
}

ss <- function(d) {
  sum(d[lower.tri(d)]^2)
}

cd <- function(m,xlabs,ylabs) {
  x <- as.matrix(m[xlabs,xlabs])
  n.x <- nrow(x)
  y <- as.matrix(m[ylabs,ylabs])
  n.y <- nrow(y)
  xy <- m[c(xlabs,ylabs),c(xlabs,ylabs)]
  return(sqrt( (ss(xy) - (n.x + n.y) * (ss(x)/n.x + ss(y)/n.y)) / (n.x*n.y)))
}

getCdCommon <- function(m,s,common1,common2) {
  return(cd(m,
            rownames(m)[rownames(m) %in% s$SampleID[s$common==common1]],
            rownames(m)[rownames(m) %in% s$SampleID[s$common==common2]]))
}

makeLegendColors <- function(vector) {
  legendColors <- rainbow(length(unique(vector)))
  names(legendColors) <- sort(unique(vector))
  return(legendColors)
}

makeLegendShapes <- function(vector) {
  if(length(unique(vector)) > 8) { stop("Can't make order legend for more than 8 unique values.")}
  vectorOrder <- as.matrix(unique(cbind(as.character(vector),as.numeric(factor(vector)))))
  orderShapeList <- c(15,17,18,19,21,22,23,24,25)
  newVectorOrder <- as.numeric(orderShapeList[as.numeric(vectorOrder[,2])])
  names(newVectorOrder) <- vectorOrder[,1]
  return(newVectorOrder)
}

makeHeatmapColors <- function (..., na.value = "gray", zero.value = "white",
                               low.value = "dark blue",threshold = 0.4, scale = 1) {
  breaks = c(0,0.000001,0.0001,seq(0.01,1,by=0.01))
  rainbow_colors <- rev(rainbow(100 * threshold, start = 0, 
                                end = 0.6))
  last_color <- tail(rainbow_colors, n = 1)
  lastColorReps <- length(breaks) * (1 - threshold)
  colors <- c(zero.value, low.value, rainbow_colors, rep(last_color, lastColorReps))
  ggplot2::scale_fill_gradientn(na.value = na.value, colors = colors, 
                                limits = c(0, 1) * scale,
                                values = breaks,
                                breaks = breaks,
                                ...)
}

makeHeatmapColorsLog <- function (..., na.value = "white", scale = 1, threshold = 0.4, 
                                  low.value = "dark blue", very.low.value = "black",
                                  nrainbowcolors = 8) {
  breaks = log10(c(1e-8,1e-6,1e-4,seq(0.01,1,by=0.05)))
  labels = c(expression("10"^{-8}),expression("10"^{-6}),expression("10"^{-4}),"0.01")
  rainbow_colors <- rev(rainbow(nrainbowcolors, start = 0, end = (1-threshold)))
  last_color <- tail(rainbow_colors, n = 1)
  lastColorReps <- length(breaks) - nrainbowcolors - 2
  colors <- c(very.low.value, low.value, rainbow_colors,
              rep(last_color, lastColorReps))
  legendBreaks = breaks[1:length(labels)]
  ggplot2::scale_fill_gradientn(na.value = na.value, colors = colors, 
                                limits = log10(c(1e-8, 1) * scale), labels = labels,
                                values = 1 + breaks / 8, breaks = legendBreaks,
                                ...)
}

isoTax <- function(metadata, tax, counts) {
  taxOtus <- names(metadata[grep(tax,metadata)])
  percentBySample <- apply(counts, 2, function(X) sum(X[names(X) %in% taxOtus]) / sum(X))
  return(percentBySample)
}



# Export OTU tables----
# Top 10000 OTUs by total abundance for used samples + controls (311 total)

resetData("toUsePlusControls")
G$o$prop <- apply(G$o$counts, 2, function(x) ifelse (x, x / sum(x),0))
G$maxOtuPropAll <- apply(G$o$prop, 1, max)
G$finalOtuTable <- G$o$counts[names(sort(G$maxOtuPropAll, decreasing=T))[1:10000],]
write.table(cbind(G$o$metadata[rownames(G$finalOtuTable)],G$finalOtuTable),
            paste0(work_dir, "otuTable.final.csv"),
            row.names=T,sep=",",quote = FALSE)




# Run Unifrac----
resetData("toUse")
system(paste0("beta_diversity.py -i ",work_dir,"otu/otu_table.biom -o ",work_dir,
              "bd_wnu -t ",work_dir,"otu/rep_set.tre -m weighted_normalized_unifrac"))

G$o$rare <- apply(G$o$counts,2,function(x) rarefyCounts(x,1000))
rownames(G$o$rare) <- rownames(G$o$counts)
write.table(cbind(G$o$rare,G$o$metadata),paste0(work_dir,"otu/otu_rare.txt"), quote = F, sep = "\t")
system(paste0("sed -i '1s/^/# Constructed from biom file\\n#OTU ID\t/' ",work_dir,"otu/otu_rare.txt"))
system(paste0("sed -i '2s/$/Consensus Lineage/' ",work_dir,"otu/otu_rare.txt"))
system(paste0("biom convert -i ",work_dir,"otu/otu_rare.txt -o ",work_dir,
              "otu/otu_rare.biom --table-type=\"OTU table\" --to-hdf5"))
system(paste0("beta_diversity.py -i ",work_dir,"otu/otu_rare.biom -o ",work_dir,
              "bd_uur -t ",work_dir,"otu/rep_set.tre -m weighted_unifrac"))

system(paste0("alpha_diversity.py -i ",work_dir,"otu/otu_table.biom -o ",work_dir,
              "alpha_diversity.tsv -t ",work_dir,"otu/rep_set.tre -m PD_whole_tree"))

# Figure 1----

resetData("toUse")

library(reshape2)
library(ggplot2)
library(phyloseq)
#plot(newick)
#plotTree(newick,node.numbers=TRUE)
#plot(newick,show.node.label=TRUE)
#edgelabels(round(newick$edge.length,0), font=1)
#plot_tree(newick,"treeonly",nodeplotblank,base.spacing=0.55)

newickPlot <- plot_tree(G$newick,"treeonly",nodeplotblank,base.spacing=0.55)

# Heatmap of OTUs by phylum grouped by species

G$aPhylum <- sub("; c__(.)*$", "", G$o$metadata, perl=T)
G$oByPhylum <- rowsum(G$o$counts,G$aPhylum)

G$obpProp <- apply(G$oByPhylum,2,function(x) ifelse (x,x / sum(x),0))
G$obpProp <- data.frame(t(G$obpProp[rowSums(G$obpProp)>0,]))

G$obpProp$common <- as.character(G$s$common[match(rownames(G$obpProp),G$s$SampleID)])
G$obpPropBySpecies <- aggregate(.~common,data=G$obpProp,mean)

G$obpMaxProps <- apply(G$obpPropBySpecies[,-1],2,max)

numOtuPhylaToUse <- 10
G$otuPhylaToUse <- names(sort(G$obpMaxProps,decreasing=TRUE))[1:(numOtuPhylaToUse-1)]

# Write to file so validation plots can use same set of phyla
write.table(G$otuPhylaToUse, paste0(work_dir, "otu_phyla.txt"), quote = F, row.names = F, col.names = F)

G$obpPropBySpeciesToUse <- G$obpPropBySpecies[,names(G$obpPropBySpecies) == 'common' |
                                                names(G$obpPropBySpecies) %in% G$otuPhylaToUse]

colnames(G$obpPropBySpeciesToUse)[-1] <-
  colnames(G$obpPropBySpeciesToUse)[rev(order(colSums(G$obpPropBySpeciesToUse[-1])))+1]
G$obpPropBySpeciesToUse[,-1] <-
  G$obpPropBySpeciesToUse[,rev(order(colSums(G$obpPropBySpeciesToUse[-1])))+1]

G$obpPropBySpeciesToUse$Other <- 1 - rowSums(G$obpPropBySpeciesToUse[-1])

G$obpHeatmapData <- melt(cbind(G$obpPropBySpeciesToUse), id.vars=c('common'))

fill.scale <- makeHeatmapColors()
plotHeatmap <- ggplot(G$obpHeatmapData, aes(x=variable, y=common, fill=value)) +
  scale_y_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  geom_tile(color="grey80", size=0.4) + theme_grey() +
  theme(
    legend.position = "none",
    strip.text.y = element_text(angle=0, vjust=0),
    strip.text.x = element_text(angle=90, vjust=0),
    panel.border = element_blank(),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(angle = 45, hjust = 1, size=10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) +
  fill.scale


# Bar chart of max OTU

G$o$prop <- countsToProp(G$o$counts)

G$o$maxProp <- apply(G$o$prop, 2, max)
G$o$maxProp <- G$o$maxProp[match(G$s$SampleID[G$sampleOrdering],
                                 names(G$o$maxProp))]
G$maxPropBySpecies <-
  tapply(G$o$maxProp[match(G$s$SampleID, names(G$o$maxProp))],
         G$s$common,mean)
G$maxPropBySpecies <-
  G$maxPropBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                           names(G$maxPropBySpecies))]

G$mpbsDf <- data.frame(common=names(G$maxPropBySpecies),value=unname(G$maxPropBySpecies))
G$mpDf <- data.frame(common=G$s$common[match(names(G$o$maxProp),G$s$SampleID)],
                     value=G$o$maxProp)
plotMax <- ggplot(G$mpbsDf,aes(x=common,y=value,fill="cyan")) +
  geom_bar(stat="identity",color="cyan") +
  scale_fill_manual(values=c("cyan")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="[title]") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$mpDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


# Bar chart of shannon diversity

G$otuShannon <- shannonDiv(G$o$counts)
G$otuShannonDf <- data.frame(SampleID=as.character(G$otuShannon[,1]),
                             shannon=as.numeric(G$otuShannon[,2]))
G$otuShannonDf$common <- G$s$common[match(G$s$SampleID,G$otuShannonDf$SampleID)]
G$otuShannonBySpecies <- aggregate(shannon ~ common, data=G$otuShannonDf, mean)
G$otuShannonBySpecies <-
  G$otuShannonBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                              G$otuShannonBySpecies$common),]

plotShannon <- ggplot(G$otuShannonBySpecies,aes(x=common,y=shannon, fill="purple")) +
  geom_bar(stat="identity",color="purple") +
  scale_fill_manual(values=c("purple")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="[title]") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$otuShannonDf,aes(x=common,y=shannon),color="black",shape="+",size=4,show.legend = FALSE)


# Bar chart of percent reads unassigned

G$unassignedOtuIds <- names(G$a)[G$a=="Unassigned"]
G$percentTotalUnassignedBySample <-
  apply(G$o$counts,2,function(X) sum(X[names(X) %in% G$unassignedOtuIds]) / sum(X))
G$percentTotalUnassignedBySample <-
  G$percentTotalUnassignedBySample[match(G$s$SampleID[G$sampleOrdering],
                                         names(G$percentTotalUnassignedBySample))]
G$s$percentTotalUnassigned <-
  G$percentTotalUnassignedBySample[match(G$s$SampleID, names(G$percentTotalUnassignedBySample))]


G$percentTotalUnassignedBySpecies <-
  tapply(G$percentTotalUnassignedBySample[match(G$s$SampleID, names(G$percentTotalUnassignedBySample))],
         G$s$common,mean)
G$percentTotalUnassignedBySpecies <-
  G$percentTotalUnassignedBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                          names(G$percentTotalUnassignedBySpecies))]

G$ptubsDf <- data.frame(common=names(G$percentTotalUnassignedBySpecies),
                        value=unname(G$percentTotalUnassignedBySpecies))
G$ptuDf <- data.frame(common=G$s$common[match(names(G$percentTotalUnassignedBySample),G$s$SampleID)],
                      value=G$percentTotalUnassignedBySample)

plotUnassigned <- ggplot(G$ptubsDf,aes(x=common,y=value,fill="red")) +
  geom_bar(stat="identity",color="red") +
  scale_fill_manual(values=c("red")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="[title]") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$ptuDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


# Barchart of phylogenetic diversity

G$phyDivBySpecies <- tapply(G$s$phyDiv, G$s$common,mean)
G$phyDivBySpecies <-
  G$phyDivBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                          names(G$phyDivBySpecies))]

G$pdbsDf <- data.frame(common=names(G$phyDivBySpecies),
                        value=unname(G$phyDivBySpecies))
G$pdDf <- data.frame(common=G$s$common,
                      value=G$s$phyDiv)

plotPhyDiv <- ggplot(G$pdbsDf,aes(x=common,y=value,fill="gray")) +
  geom_bar(stat="identity",color="gray") +
  scale_fill_manual(values=c("gray")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="[title]") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$pdDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


# Barchart of Wolbachia OTUs

G$wolbachiaOtus <- names(G$a[grep("Wolbachia",G$a)])
G$percentWolbachiaBySample <-
  apply(G$o$counts,2,function(X) sum(X[names(X) %in% G$wolbachiaOtus]) / sum(X))
G$percentWolbachiaBySample <-
  G$percentWolbachiaBySample[match(G$s$SampleID[G$sampleOrdering],
                                   names(G$percentWolbachiaBySample))]

G$percentWolbachiaBySpecies <-
  tapply(G$percentWolbachiaBySample[match(G$s$SampleID, names(G$percentWolbachiaBySample))],
         G$s$common,mean)
G$percentWolbachiaBySpecies <-
  G$percentWolbachiaBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                          names(G$percentWolbachiaBySpecies))]

G$pwbsDf <- data.frame(common=names(G$percentWolbachiaBySpecies),
                       value=unname(G$percentWolbachiaBySpecies))
G$pwDf <- data.frame(common=G$s$common[match(names(G$percentWolbachiaBySample),G$s$SampleID)],
                     value=G$percentWolbachiaBySample)
plotWolbachia <- ggplot(G$pwbsDf,aes(x=common,y=value, fill="orange")) +
  geom_bar(stat="identity",color="orange") +
  scale_fill_manual(values=c("orange")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="[title]") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$pwDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)




# Barchart of obligate anaerobes

G$anaerobeTable <-read.csv(paste0(work_dir,"bilaterian_aerotolerance.csv"),header=T,
                           na.strings="<NA>",stringsAsFactors = F)
G$anaerobeTable$SampleID[G$anaerobeTable$SampleID=="Bonobo.BI0331.Feces.1.1"] <- "Pig.BI0331.Feces.1.1"
G$anaerobeTable$speciesOrder <- G$s[G$anaerobeTable$SampleID,"speciesOrder"]
G$anaerobeTable <-
  G$anaerobeTable[!G$anaerobeTable$aerobic_status_corrected=="NA",]
G$anaerobeTable <-
  G$anaerobeTable[!G$anaerobeTable$aerobic_status_corrected=="no taxonomic assignment",]
G$anaerobeTable <-
  G$anaerobeTable[!G$anaerobeTable$aerobic_status_corrected=="Streptophyta",]
G$anaerobeTable$Proportion <-
  ave(G$anaerobeTable$Counts,G$anaerobeTable$SampleID,
      FUN=function(X) X / sum(X))


G$obligateAnaerobeTable <-
  G$anaerobeTable[G$anaerobeTable$aerobic_status_corrected=="obligate anaerobe",]
G$s <- G$s %>% left_join(G$obligateAnaerobeTable[,c("SampleID","Proportion")] %>%
                           rename(obligateAnaerobeProp = Proportion)) %>%
  mutate_at(.vars=vars(obligateAnaerobeProp),
            funs(replace(obligateAnaerobeProp, which(is.na(.)), 0)))

G$obligateAnaerobeSpecies <- tapply(G$s$obligateAnaerobeProp,G$s$common,mean)
G$obligateAnaerobeSpecies <-
  G$obligateAnaerobeSpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                  names(G$obligateAnaerobeSpecies))]

G$oapbsDf <- data.frame(common=names(G$obligateAnaerobeSpecies),
                        value=unname(G$obligateAnaerobeSpecies))
G$oapDf <- data.frame(common=G$s$common,
                      value=G$s$obligateAnaerobeProp)

plotObligateAnaerobes <- ggplot(G$oapbsDf,aes(x=common,y=value,fill="pink")) +
  geom_bar(stat="identity",color="pink") +
  scale_fill_manual(values=c("pink")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="[title]") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$oapDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


# Print out figure

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf(paste0(work_dir, "R_plots/figure1.final.new.pdf"),
    height=15, width=20, onefile = FALSE)
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1000, 14)))
print(newickPlot, vp = vplayout(1:935, 1:2))
print(plotHeatmap, vp = vplayout(35:1000, 5:8))
print(plotMax, vp = vplayout(35:923, 9))
print(plotShannon, vp = vplayout(35:923, 10))
print(plotPhyDiv, vp = vplayout(35:923, 11))
print(plotObligateAnaerobes, vp = vplayout(35:923, 12))
print(plotWolbachia, vp = vplayout(35:923, 13))
print(plotUnassigned, vp = vplayout(35:923, 14))
dev.off()


# Figure 2 (tSNE)----
# PCoA, tsne, heatmaps for final figures

resetData("toUse")
resetUnifracData()

set.seed(1001)
library(Rtsne)
G$uur_tsne <- Rtsne(G$uur,is_distance=TRUE,perplexity=10,max_iter=3000)

#G$classColors <- makeLegendColors(G$s$class)
G$classColors <- structure(
  c('#8dd3c7','#ffed6f','#bc80bd','#fccde5','#80b1d3','#fdb462','#b3de69','#fb8072','#bebada','#ccebc5'),
  .Names = c("Actinopteri", "Aves", "Chondrichthyes", "Diplopoda", "Holothuroidea", "Insecta", "Malacostraca", "Mammalia", "Polychaeta", "Reptilia")
)
G$s$classColor <- G$classColors[G$s$class]

G$classOrderNumber <-
  unlist(sapply(unique(G$s$class),function(X) makeLegendShapes(G$s$order[G$s$class==X])))

G$s$orderShape <- as.numeric(G$classOrderNumber[as.character(G$s$order)])

pdf(paste0(work_dir, "R_plots/figure2.tsne.uur.pdf"),
    width=11,height=8,useDingbats = FALSE)
par(mar=c(4,4,1,12)) # need to leave margin room
plot(G$uur_tsne$Y, col=G$s$classColor, pch=as.numeric(G$s$orderShape),
     xlab="t-SNE 1",ylab="t-SNE 2",lwd=2)
plotTsneLegends()
dev.off()

pdf(paste0(work_dir, "R_plots/figure2.tsne.uur.extra.pdf"),
    width=11,height=8,useDingbats = FALSE)
par(mar=c(4,4,1,12)) # need to leave margin room
plot(G$uur_tsne$Y, col=G$s$classColor, pch=as.numeric(G$s$orderShape),
     xlab="t-SNE 1",ylab="t-SNE 2",lwd=2)
plotExtraTsneLegends()
dev.off()




# Figure 1G/H----
# plotDNA for unassigned vs assigned reads, both 5' and 3'

library(dnaplotr)
library(dnar)
library(ape)

GH <- new.env()

GH$otuTree <- read.tree(paste0(work_dir, "rep_set.tre"))

# This file is not included in our repository, as it is quite large (>4GB)
# You will need to run Qiime (instructions above) and locate this file under
# the "otu/pynast_aligned_seqs/" directory
GH$repSetAligned <-
  read.fa(paste0(work_dir, "rep_set_aligned.fasta"))
rownames(GH$repSetAligned) <- sapply(GH$repSetAligned$name,function(X) strsplit(X," ")[[1]][1])

GH$repSetAA <- GH$repSetAligned[!rownames(GH$repSetAligned) %in% names(GH$a)[GH$a=="Unassigned"],c("name","seq")]

#GH$repSetAAordered <-
#  GH$repSetAA[match(GH$otuTree$tip.label[GH$otuTree$tip.label %in% rownames(GH$repSetAA)],
#                         rownames(GH$repSetAA)),]

set.seed(1001)
GH$repSetAAsample <- GH$repSetAA[sample(nrow(GH$repSetAA),10000),]

GH$repSetAAsample1 <- GH$repSetAA[sample(nrow(GH$repSetAA),50),]

# code to identify alignment positions with < 10% gaps
# For each position (1-~7000), count # of gaps
# Return True if < 10% of positions are gaps, False otherwise

maxGaps <- 0.9 * nrow(GH$repSetAAsample1)
positionsToKeep <-
  sapply(1:nchar(GH$repSetAAsample1$seq[1]),
         function(X) { return(
           sum(sapply(GH$repSetAAsample1$seq,
                      function(Y) { unlist(strsplit(Y, ""))[X] }) == "-") < maxGaps) })

#GH$repSetAAordered$seq2 <-
#  sapply(GH$repSetAAordered$seq,
#         function(X) { paste0(unlist(strsplit(X,""))[positionsToKeep],collapse="") })
#plotDNA(GH$repSetAAordered$seq2, res=4000)


GH$repSetAAsample$seq2 <-
  sapply(GH$repSetAAsample$seq,
         function(X) { paste0(unlist(strsplit(X,""))[positionsToKeep],collapse="") })

#GH$repSetAAsampleOrdered <-
#  GH$repSetAAsample[match(GH$otuTree$tip.label[GH$otuTree$tip.label %in% rownames(GH$repSetAAsample)],
#                         rownames(GH$repSetAAsample)),]

#plotDNA(GH$repSetAAsample$seq2, res=4000)
#plotDNA(GH$repSetAAsampleOrdered$seq2, res=4000)
#plotDNA(GH$repSetAAsampleOrdered$seq2[c(700:1100,5050:5600)], res=4000)

#sampleSeqNames <- sapply(GH$repSetAAsample$name, function(X) {paste0(">",paste(strsplit(X," ")[[1]][1:2],collapse=" "))})
#write(sampleSeqNames,file = paste0(work_dir, "repSetSample.txt")

GH$repSetAU <- GH$repSetAligned[rownames(GH$repSetAligned) %in% names(GH$a)[GH$a=="Unassigned"],c("name","seq")]
GH$repSetAUsample <- GH$repSetAU[sample(nrow(GH$repSetAU),10000),]

GH$repSetAUsample$seq2 <-
  sapply(GH$repSetAUsample$seq,
         function(X) { paste0(unlist(strsplit(X,""))[positionsToKeep],collapse="") })

#GH$repSetAUsampleOrdered <-
#  GH$repSetAUsample[match(GH$otuTree$tip.label[GH$otuTree$tip.label %in% rownames(GH$repSetAUsample)],
#                         rownames(GH$repSetAUsample)),]


#save.image(paste0(work_dir, "R//good.20170830.RData")
load(paste0(work_dir, "R//good.20170830.RData"))

pdf(paste0(work_dir, "R_plots/figure1.panelG.pdf"),
    height=4, width=4, onefile = FALSE)
plotDNA(GH$repSetAAsample$seq2, res=4000)
#plotDNA(GH$repSetAAsampleOrdered$seq2, res=4000)

dev.off()
pdf(paste0(work_dir, "R_plots/figure1.panelH.pdf"),
    height=4, width=4, onefile = FALSE)
plotDNA(GH$repSetAUsample$seq, res=4000)
dev.off()





# Figure A1----

resetData("toUse")

G$sampleMetadataDf <- do.call(rbind,
                              list(data.frame(level="Phylum",
                                              count=table(G$s$phylum[G$sampleOrdering])),
                                   data.frame(level="Class",
                                              count=table(G$s$class[G$sampleOrdering])),
                                   data.frame(level="Order",
                                              count=table(G$s$order[G$sampleOrdering])),
                                   data.frame(level="Family",
                                              count=table(G$s$family[G$sampleOrdering])),
                                   data.frame(level="Species",
                                              count=table(G$s$common[G$sampleOrdering]))))

G$sampleMetadataDf$count.Var1 <-
  factor(G$sampleMetadataDf$count.Var1,
         levels = c(unique(as.character(G$speciesTax$phylum[G$figure1SpeciesOrdering])),
                    unique(as.character(G$speciesTax$class[G$figure1SpeciesOrdering])),
                    unique(as.character(G$speciesTax$order[G$figure1SpeciesOrdering])),
                    unique(as.character(G$speciesTax$family[G$figure1SpeciesOrdering])),
                    unique(as.character(G$speciesTax$common[G$figure1SpeciesOrdering]))))
G$sampleMetadataDf <- G$sampleMetadataDf[match(levels(G$sampleMetadataDf$count.Var1),
                                               G$sampleMetadataDf$count.Var1),]
G$sampleMetadataDf <- G$sampleMetadataDf %>% group_by(level) %>%
  mutate(cum.freq = cumsum(count.Freq) - 0.5*count.Freq)
pal <- sapply(table(G$sampleMetadataDf$level),
              function(X) colorRampPalette(c('red','blue','dark green'))(X))
pal <- unname(unlist(pal))

palSet <- colorRampPalette(c('blue','dark green'))(nrow(G$s))
pal <- sapply(unique(G$sampleMetadataDf$level),
              function(X) palSet[G$sampleMetadataDf$cum.freq[G$sampleMetadataDf$level==X]+0.5])
pal <- unname(unlist(pal))

G$sampleMetadataDf$labelSize <- sapply(G$sampleMetadataDf$count.Freq, function(X) 2 * min(X,3))

pdf(paste0(work_dir, "R_plots/A1.final.pdf"),
    height=30, width=15, onefile = FALSE)
ggplot(G$sampleMetadataDf, aes(x = level, y = count.Freq, fill = count.Var1)) + 
  geom_bar(stat = "identity", show.legend=FALSE) + 
  scale_x_discrete(limits = c("Phylum","Class","Order","Family","Species"))  + 
  scale_fill_manual(values=pal) +
  theme_bw() + labs(y="Percent of OTUs in Wolbachia") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(), axis.title.x=element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank()) +
  geom_text(aes(label=paste0(count.Var1,": ", count.Freq),
                y=max(cum.freq)-cum.freq+0.5), colour="white",size=G$sampleMetadataDf$labelSize) +
  scale_y_continuous(labels=waiver(), trans="reverse")
dev.off()


# Figure A2----

resetData("all")

pdf(paste0(work_dir, "R_plots/A2.boxplot.readCounts.pdf"),
    width=11,height=8,useDingbat=FALSE)
par(mar=c(5,7,6,6))
par(cex.lab=2)
par(cex.main=2)
boxplot(readCount ~ sample_type2, data = G$s, log = "y", yaxt="n",
        xlab="Sample type", ylab="Number of reads")
logAxis(2,las=1,axisMin = 1)
dev.off()


# Figure A3----
# Reads by OTU Phylum

resetData("toUse")

#load(paste0(work_dir, "R//good.20170830.RData")

G$aPhylum <- sub("; c__(.)*$", "", G$o$metadata, perl=T)
G$oByPhylum <- rowsum(G$o$counts,G$aPhylum)

dataToPlot <- sort(rowSums(G$oByPhylum), decreasing = T)[1:9]
dataToPlot <- data.frame(value=dataToPlot,phylum=names(dataToPlot))

dataToPlot$phylum <- with(dataToPlot, factor(phylum, levels=c("Other",as.character(phylum[order(value)]))))
dataToPlot[10,] <- cbind(sum(sort(rowSums(G$oByPhylum), decreasing = T)[-c(1:9)]), "Other")

dataToPlot$value <- as.numeric(dataToPlot$value) / 1e6

pdf(paste0(work_dir, "R_plots/supp3.pdf"), width = 10, height = 7)
ggplot(dataToPlot, aes(x=phylum, y=value)) + geom_bar(stat="identity") +
  coord_flip() + theme_bw() +
  labs(y = "Number of reads (in millions)", x = "OTU Phylum") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
dev.off()





### Figures A4-A13 (Val.)----

VA <- new.env()
resetValidationData(VA, paste0(work_dir, "validation/anteater/"))
plotFigure1(VA, "anteater")

VBi <- new.env()
resetValidationData(VBi, paste0(work_dir, "validation/bird/"))
plotFigure1(VBi, "bird")

VBu <- new.env()
resetValidationData(VBu, paste0(work_dir, "validation/bug/"))
plotFigure1(VBu, "bug")

VBR <- new.env()
resetValidationData(VBR, paste0(work_dir, "validation/bugRev/"))
plotFigure1(VBR, "bugRev")

VF <- new.env()
resetValidationData(VF, paste0(work_dir, "validation/fish/"))
plotFigure1(VF, "fish")

VM <- new.env()
resetValidationData(VM, paste0(work_dir, "validation/muegge/"))
plotFigure1(VM, "muegge454")

VMI <- new.env()
resetValidationData(VMI, paste0(work_dir, "validation/mueggeIllumina/"))
plotFigure1(VMI, "mueggeIllumina")

VP <- new.env()
resetValidationData(VP, paste0(work_dir, "validation/primates/"))
plotFigure1(VP, "primates")

VW <- new.env()
resetValidationData(VW, paste0(work_dir, "validation/whale454/"))
VW$s$phyDiv <- VW$pdwt[VW$s$old_id, "PD_whole_tree"]
plotFigure1(VW, "whale454")

VWI <- new.env()
resetValidationData(VWI, paste0(work_dir, "validation/whaleIllumina/"))
VWI$s$phyDiv <- VWI$pdwt[VWI$s$old_id, "PD_whole_tree"]
plotFigure1(VWI, "whaleIllumina")







# Figure A14-A15----
# Aerobic status of bacteria

G$anaerobeTable <- read.csv(paste0(work_dir,"bilaterian_aerotolerance.csv"),header=T,
                           na.strings="<NA>",stringsAsFactors = F)
G$anaerobeTable <- G$anaerobeTable %>%
  left_join(G$s[,c("SampleID","speciesOrder")]) %>%
  filter(!aerobic_status_corrected %in% c("NA","no taxonomic assignment","Streptophyta")) %>%
  mutate(Proportion=ave(Counts,SampleID,
                        FUN=function(X) X / sum(X)))


G$aerobeTable <-
  G$anaerobeTable[G$anaerobeTable$aerobic_status_corrected=="aerobe",]
G$s <- G$s %>% left_join(G$aerobeTable[,c("SampleID","Proportion")] %>%
                           rename(aerobeProp = Proportion)) %>%
  mutate_at(.vars=vars(aerobeProp),
            funs(replace(aerobeProp, which(is.na(.)), 0)))

G$s$aerobeToObligateAnaerobeRatio <- log10((G$s$aerobeProp + 0.00001) / (G$s$obligateAnaerobeProp + 0.00001))

G$classColors <- makeLegendColors(G$s$class)
G$s$classColor <- G$classColors[G$s$class]

G$s <- G$s %>%
  left_join(G$anaerobeTable %>%
              filter(aerobic_status_corrected=="obligate anaerobe") %>%
              select(c(SampleID,Proportion)) %>%
              rename(obligateAnaerobeProp = Proportion)) %>%
  mutate_at(.vars=vars(obligateAnaerobeProp),
            funs(replace(obligateAnaerobeProp, which(is.na(.)), 0)))

G$s$aerobeToObligateAnaerobeRatio <-
  log10((G$s$aerobeProp + 0.00001) / (G$s$obligateAnaerobeProp + 0.00001))

G$anaerobeTableSpecies <- aggregate(Counts ~ common + aerobic_status_corrected + speciesOrder,
                                    G$anaerobeTable,sum)
G$anaerobeTableSpecies$Proportion <-
  ave(G$anaerobeTableSpecies$Counts,G$anaerobeTableSpecies$common,
      FUN=function(X) X / sum(X))

pdf(paste0(work_dir, "R_plots/A14_aerobic_status_species.pdf"),
    height=11, width=8, onefile = FALSE)
G$anaerobeTableSpecies %>%
  mutate(ProportionNA = ifelse(is.na(aerobic_status_corrected), Proportion, 0)) %>%
  droplevels() %>%
  mutate(common = reorder(common, speciesOrder)) %>%
  ggplot() +
  geom_bar(aes(x=common, y=Proportion, fill=aerobic_status_corrected), stat="identity", position = "stack") +
  coord_flip() +
  scale_fill_brewer(palette="Paired", na.value="#CCCCCC") +
  theme_bw()
dev.off()


pdf(paste0(work_dir, "R_plots/A15_aerobe_ratio_samples.pdf"),
    height=8, width=11, onefile = FALSE)
plot(G$s$log.weight,G$s$aerobeToObligateAnaerobeRatio,
     col=G$s$classColor,xlab="Weight (grams, log10)",
     ylab="Ratio of aerobic reads to anaerobic reads",
     main="Ratio of aerobic/anaerobic reads vs sample weight.\nColored by host class.")
dev.off()

G$species$aerobeToObligateAnaerobeRatio <-
  with(aggregate(aerobeToObligateAnaerobeRatio ~ common,G$s, mean),
       aerobeToObligateAnaerobeRatio[match(G$species$common,common)])
G$species$log.weight <-
  with(aggregate(log.weight ~ common,G$s, mean),
       log.weight[match(G$species$common,common)])
G$classColors<-structure(
  c('#8dd3c7','#ffed6f','#bc80bd','#fccde5','#80b1d3','#fdb462','#b3de69','#fb8072','#bebada','#ccebc5'),
  .Names = c("Actinopteri", "Aves", "Chondrichthyes", "Diplopoda", "Holothuroidea", "Insecta", "Malacostraca", "Mammalia", "Polychaeta", "Reptilia")
)
G$species$classColor <- G$classColors[G$species$class]



library(dnar)
pdf(paste0(work_dir, "R_plots/A15_aerobe_ratio_species.pdf"),
    height=8, width=11, onefile = FALSE,useDingbats = FALSE)
plot(10^G$species$log.weight,10^G$species$aerobeToObligateAnaerobeRatio,
     log='xy',xaxt='n',yaxt='n', bg=G$species$classColor, pch=21, cex=1.2,
     xlab='Animal weight (g)',
     ylab='Average ratio of aerobic/anaerobic reads')
logAxis(1)
logAxis(2,las=1)
text(10^G$species$log.weight,0.7*10^G$species$aerobeToObligateAnaerobeRatio,
     rownames(G$species),cex=.5)
thisLm <- lm(aerobeToObligateAnaerobeRatio~log.weight,G$species)
fakeWeights <- seq(min(G$species$log.weight,na.rm=TRUE)-1,
                   max(G$species$log.weight)+1, length.out=1000)
preds <- predict(thisLm,data.frame('log.weight'=fakeWeights),interval='conf')
lines(10^fakeWeights,10^preds[,'fit'])
polygon(10^c(fakeWeights, rev(fakeWeights)),
        10^c(preds[,'lwr'],rev(preds[,'upr'])), border=NA, col='#00000022')
dev.off()


# Figure A16----
# Phylogenetic Diversity

pdf("/home/kevin/projects/bilateria/R_plots/diet_vs_phyDiv.mammals.pdf",
    height=8, width=11, onefile = FALSE)
G$s %>% filter(class=="Mammalia") %>%
  ggplot(aes(x=diet, y=phyDiv)) + geom_boxplot(aes(fill=diet)) +
  scale_y_log10() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Phylogenetic diversity of mammals by diet")
dev.off()

pdf("/home/kevin/projects/bilateria/R_plots/diet_vs_phyDiv.arthropods.pdf",
    height=8, width=11, onefile = FALSE)
G$s %>% filter(phylum=="Arthropoda") %>%
  ggplot(aes(x=diet, y=phyDiv)) + geom_boxplot(aes(fill=diet)) +
  scale_y_log10() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Phylogenetic diversity of arthropods by diet")
dev.off()


library(ggplot2)
library(scales)

G$host.el <- data.frame(el=G$newick$edge.length)

ggplot(G$el, aes(x=el)) +
  stat_density(aes(y=..count..), color="black", fill="blue", alpha=0.3) +
  scale_x_continuous(breaks=c(0,1,10,30,100,300,1000), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks=c(0,5,10,25,50), expand=c(0,0)) +
  theme_bw()

library(phyloseq)
library(dnaplotr)
library(dnar)
library(ape)

G$otuTree <- read.tree(paste0(work_dir, "otu/rep_set.tre"))

G$otu.el <- data.frame(el=G$otuTree$edge.length)

ggplot(G$otu.el, aes(x=el)) +
  stat_density(aes(y=..count..), color="black", fill="blue", alpha=0.3) +
  scale_x_continuous(breaks=c(0,0.001,0.01,0.1,0.5), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks=c(0,100,1000,10000,100000), trans="log",expand=c(0,0)) +
  theme_bw()

ggplot(G$otu.el,aes(x=el)) + geom_histogram(colour="darkblue", size=1, fill="blue") +
  scale_x_log10(breaks = c(0,0.001,0.01,0.1,0.5,1))

#sapply(seq(1,length(tempTree$tip.label),100000), function(X) maxEdgeLength(G$otuTree, X))

#library(parallel)

#cluster <- makeCluster(24)
#tempTree <- G$otuTree
#clusterExport(cl=cluster, varlist=c("maxEdgeLength", "tempTree"))
#maxOtuEdgeLengths <- parLapply(cluster, seq(1,length(tempTree$tip.label),100000),
#                               function(X) maxEdgeLength(tempTree,X))
#stopCluster(cluster)

G$remoteOtus <- unlist(sapply(G$otuTree$edge[G$otuTree$edge.length > 0.1,2],
                              function(X) getEdgeLeaves(G$otuTree,X), simplify=TRUE))
G$RONames <- G$otuTree$tip.label[G$remoteOtus]
G$RONames <- unique(G$RONames[G$RONames %in% G$o$otu_ids])
G$ROTaxa <- G$a[G$RONames]

G$ROUnassigned <- G$RONames[G$ROTaxa=="Unassigned"]

# Heatmap of remote OTUs

# Do we want "% OTUS by Phylum for all REMOTE Otus"...
#G$ROPhyla <- sub("; c__(.)*$", "", G$o$metadata[G$RONames], perl=T)
#G$ROCounts <- G$o$counts[G$RONames,]
#G$ROCbyPhylum <- rowsum(G$ROCounts,G$ROPhyla)

# or "% OTUs by Phylum for ALL Otus" ?
G$ROmetadata <- G$a
G$ROmetadata <- sapply(G$ROmetadata, function(X) "non-remote")
G$ROmetadata[G$RONames] <- G$a[G$RONames]
G$ROTaxa <- sub("; c__(.)*$", "", G$ROmetadata, perl=T)
G$ROTaxa[G$ROTaxa=="Unassigned"] <- make.unique(G$ROTaxa[G$ROTaxa=="Unassigned"])
G$ROCbyTaxa <- rowsum(G$o$counts,G$ROTaxa)

G$rocbtProp <- apply(G$ROCbyTaxa,2,function(x) ifelse (x,x / sum(x),0))
G$rocbtProp <- data.frame(t(G$rocbtProp[rowSums(G$rocbtProp)>0,]))
G$rocbtProp$non.remote <- NULL

G$rocbtProp$common <- as.character(G$s$common[match(rownames(G$rocbtProp),G$s$SampleID)])
G$rocbtPropBySpecies <- aggregate(.~common,data=G$rocbtProp,mean)

G$rocbtSumProps <- apply(G$rocbtPropBySpecies[,-1],2,sum)
G$rocbtCountProps <- apply(G$rocbtPropBySpecies[,-1],2,function(X) sum(X>0))

numOtuTaxaToUse <- 25
#G$ROTaxaToUse <- names(sort(G$rocbtSumProps,decreasing=TRUE))[1:(numOtuTaxaToUse-1)]
G$ROTaxaToUse <- names(sort(G$rocbtCountProps,decreasing=TRUE))[1:(numOtuTaxaToUse-1)]

G$rocbtPropBySpeciesToUse <-
  G$rocbtPropBySpecies[,names(G$rocbtPropBySpecies) == 'common' |
                         names(G$rocbtPropBySpecies) %in% G$ROTaxaToUse]

#G$rocbtPropBySpeciesToUse <-
#  G$rocbtPropBySpeciesToUse[,c(1,rev(order(colSums(G$rocbtPropBySpeciesToUse[-1])))+1)]
G$rocbtPropBySpeciesToUse <-
  G$rocbtPropBySpeciesToUse[,c(1,rev(order(colSums(G$rocbtPropBySpeciesToUse[-1]>0)))+1)]

G$rocbtPropBySpeciesToUse$Other <-
  rowSums(G$rocbtPropBySpecies[,-1]) - rowSums(G$rocbtPropBySpeciesToUse[,-1])

library(reshape2)
G$rocbtHeatmapData <- melt(cbind(G$rocbtPropBySpeciesToUse), id.vars=c('common'))
G$rocbtHeatmapData$value[G$rocbtHeatmapData$value==0] <- NA
G$rocbtHeatmapData$value <- log10(G$rocbtHeatmapData$value)

fill.scale.log <- makeHeatmapColorsLog()
library(ggplot2)
pdf(paste0(work_dir,"R_plots/heatmap.divergentOtus.pdf"),
    height=8, width=11, onefile = FALSE)
ggplot(G$rocbtHeatmapData, aes(x=variable, y=common, fill=value)) +
  scale_y_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  scale_x_discrete(limits=c(G$ROTaxaToUse,"Other")) +
  geom_tile(color="grey80", size=0.4) + theme_grey() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text.y = element_text(angle=0, vjust=0),
    strip.text.x = element_text(angle=90, vjust=0),
    panel.border = element_blank(),
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(angle = 45, hjust = 1, size=10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) +
  fill.scale.log + ggtitle("Heatmap of divergent OTUs across host species")
dev.off()



# Run blast on unassigned reads
library(dnar)
G$repSet <- read.fa(paste0(work_dir,"/otu/rep_set/seqs_rep_set.fasta"))
rownames(G$repSet) <- sapply(G$repSet$name,function(X) strsplit(X," ")[[1]][1])
G$repSet$id <- paste(">",rownames(G$repSet),sep="")

# Rep Set Remote Unassigned
G$rsru <- G$repSet[rownames(G$repSet) %in% ROUnassigned,c("id","seq")]

write.table(G$rsru,row.names = FALSE,col.names = FALSE,
            paste0(work_dir,"unassigned/rsru.fasta"),sep="\n",quote = FALSE)

system(paste0("blastn -db /home/common/blastdbs/nt ",
              "-query ", work_dir,"unassigned/rsru.fasta -outfmt 6 ",
              "-num_threads 30 -culling_limit 10 | gzip > ",
              work_dir,"/unassigned/rsru.blast_results.txt.gz"))


#Load Blast Results
G$br <- read.blast(paste0(work_dir,"/unassigned/rsru.blast_results.txt.gz"))
#blastResults$accession <- sapply(strsplit(blastResults$tName,'\\|'),'[[',3)
G$br$maxBit <- ave(G$br$score,G$br$qName,FUN=max)
# could also do e.g. >.95*blastResults$maxBit to get close hits
G$brf <- G$br[G$br$score>0.95*G$br$maxBit,]
library(taxonomizr)
G$brf$taxa <- accessionToTaxa(G$brf$tName,"/home/common/taxonomizr/accessionTaxa.sql")

taxaNodes<-read.nodes('/home/common/taxonomizr/nodes.dmp')
taxaNames<-read.names('/home/common/taxonomizr/names.dmp')

G$brf$taxonomy <-
  getTaxonomy(G$brf$taxa,taxaNodes,taxaNames,mc.cores=5)

write.csv(G$brf, row.names = FALSE, "/home/kevin/projects/placenta/blastResultsFiltered.csv")

G$brf <- read.csv("/home/kevin/projects/placenta/blastResultsFiltered.csv")
#humanOtus <-G$brf$qName[grep("Homo sapiens", G$brf$taxonomy[,"species"])]
G$humanOtus <- as.character(G$brf$qName[grep("Homo sapiens", G$brf$taxonomy.species)])


# Figure A17----
# UniFrac analysis of sample sets
resetData("toUse")
resetUnifracData()

G$s$common <- as.character(G$s$common)
G$s$common[G$s$species=="stephensi"] <- "Mosquito"
G$s$common[G$s$species=="molitor"] <- "Mealworm"
G$s$common[G$s$species=="leucas"] <- "Bull Shark"
G$s$common[G$s$species=="plumbeus"] <- "Sandbar Shark"
G$s$common[G$s$species=="brevirostris"] <- "Lemon Shark"
G$s$common[G$s$species=="cuvier"] <- "Tiger Shark"
G$s$common[G$s$species=="americana"] <- "Dagger Moth"
G$s$common[G$s$species=="trivittata"] <- "Boxelder Bug"

#G$classColors <- makeLegendColors(G$s$class)
G$classColors <- structure(
  c('#8dd3c7','#ffed6f','#bc80bd','#fccde5','#80b1d3','#fdb462','#b3de69','#fb8072','#bebada','#ccebc5'),
  .Names = c("Actinopteri", "Aves", "Chondrichthyes", "Diplopoda", "Holothuroidea", "Insecta", "Malacostraca", "Mammalia", "Polychaeta", "Reptilia")
)
G$species$classColor <- G$classColors[G$species$class]

G$uurMatrix <- as.matrix(G$uur)

G$uurDist <- matrix(nrow=length(unique(G$s$common)),
                    ncol=length(unique(G$s$common)))
G$uurDist[,] <- 0
rownames(G$uurDist) <- unique(G$s$common)
colnames(G$uurDist) <- unique(G$s$common)
for(i in rownames(G$uurDist)) {
  for(j in colnames(G$uurDist)) {
    if(i==j) {
      G$uurDist[i,j] <- 0
    } else {
      G$uurDist[i,j] <- getCdCommon(G$uurMatrix,G$s,i,j)
    }
  }
}

G$uurDist <- G$uurDist[G$figure1SpeciesOrdering,G$figure1SpeciesOrdering]

set.seed(1001)
library(Rtsne)
G$uurDist_tsne <- Rtsne(G$uurDist,is_distance=TRUE,perplexity=10,max_iter=3000)

pdf(paste0(work_dir,"R_plots/A17.uur.speciesCentroid.tsne.pdf"),
    height=8, width=11, onefile = FALSE, useDingbats = FALSE)
par(mar=c(4,4,1,12)) # need to leave margin room
plot(G$uurDist_tsne$Y,bg=G$species$classColor,
     xlab="t-SNE 1",ylab="t-SNE 2", cex=1.2, pch=21)
text(G$uurDist_tsne$Y[,1],G$uurDist_tsne$Y[,2]-0.3,
     rownames(G$uurDist), cex=0.4, xpd=NA)
dev.off()


G$wnuMatrix <- as.matrix(G$wnu)

G$wnuDist <- matrix(nrow=length(unique(G$s$common)), ncol=length(unique(G$s$common)))
G$wnuDist[,] <- 0
rownames(G$wnuDist) <- unique(G$s$common)
colnames(G$wnuDist) <- unique(G$s$common)
for(i in rownames(G$wnuDist)) {
  for(j in colnames(G$wnuDist)) {
    if(i==j) {
      G$wnuDist[i,j] <- 0
    } else {
      G$wnuDist[i,j] <- getCdCommon(G$wnuMatrix,G$s,i,j)
    }
  }
}

G$wnuDist <- G$wnuDist[G$figure1SpeciesOrdering,G$figure1SpeciesOrdering]

set.seed(1001)
library(Rtsne)
G$wnuDist_tsne <- Rtsne(G$wnuDist,is_distance=TRUE,perplexity=10,max_iter=3000)

pdf(paste0(work_dir,"R_plots/A17.wnu.speciesCentroid.tsne.pdf"),
    height=8, width=11, onefile = FALSE, useDingbats = FALSE)
par(mar=c(4,4,1,12)) # need to leave margin room
plot(G$wnuDist_tsne$Y,bg=G$species$classColor,
     xlab="t-SNE 1",ylab="t-SNE 2", cex=1.2, pch=21)
text(G$wnuDist_tsne$Y[,1],G$wnuDist_tsne$Y[,2]-0.3,
     rownames(G$wnuDist), cex=0.4, xpd=NA)
dev.off()





# DADA2 analysis----

# Before running this analysis, sample files will need to be downloaded from SRA
# and added to "SRA" directory

library(dada2)
packageVersion("dada2")

D2 <- new.env()
D2$path <- file.path(paste0(work_dir,"/demux/"))
list.files(D2$path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
D2$fnFs <- sort(list.files(D2$path, pattern="_1.fastq", full.names = TRUE))
D2$fnRs <- sort(list.files(D2$path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
D2$sample.names <- sapply(strsplit(basename(D2$fnFs), "_"), `[`, 1)

plotQualityProfile(D2$fnFs[c(1,7,13)])
plotQualityProfile(D2$fnRs[c(1,7,13)])

# Based on above plots, will trim R1 at 200 and R2 at 160

D2$filt_path <- file.path(D2$path, "filtered") # Place filtered files in filtered/ subdirectory
D2$filtFs <- file.path(D2$filt_path, paste0(D2$sample.names, "_F_filt.fastq.gz"))
D2$filtRs <- file.path(D2$filt_path, paste0(D2$sample.names, "_R_filt.fastq.gz"))

D2$out <- filterAndTrim(D2$fnFs, D2$filtFs, D2$fnRs, D2$filtRs, truncLen=c(200,160),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(D2$out)

D2$filtered.sample.names <- D2$sample.names[sapply(D2$filtFs, file.exists)]

D2$filtFs <- D2$filtFs[sapply(D2$filtFs, file.exists)]
D2$filtRs <- D2$filtRs[sapply(D2$filtRs, file.exists)]

D2$derepFs <- derepFastq(D2$filtFs, verbose=TRUE)
D2$derepRs <- derepFastq(D2$filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(D2$derepFs) <- D2$filtered.sample.names
names(D2$derepRs) <- D2$filtered.sample.names

D2$errF <- learnErrors(D2$derepFs, multithread=TRUE)
D2$errR <- learnErrors(D2$derepRs, multithread=TRUE)

plotErrors(D2$errF, nominalQ=TRUE)
plotErrors(D2$errR, nominalQ=TRUE)

D2$dadaFs <- dada(D2$derepFs, err=D2$errF, multithread=TRUE)
D2$dadaRs <- dada(D2$derepRs, err=D2$errR, multithread=TRUE)

D2$dadaFs[[1]]

D2$mergers <- mergePairs(D2$dadaFs, D2$derepFs, D2$dadaRs, D2$derepRs, verbose=TRUE)

D2$seqtab <- makeSequenceTable(D2$mergers)

table(nchar(getSequences(D2$seqtab)))

D2$seqtab.nochim <-removeBimeraDenovo(D2$seqtab, method="consensus",
                                      multithread=TRUE, verbose=TRUE)

write.table(t(D2$seqtab.nochim),file="/home/kevin/projects/islandGut/toShare/dada2.otu.tsv",
            sep="\t", quote=F, row.names = F)















# Marine bacteria----

resetData("toUse")

marineBacteriaStrings <- c("Synechococcales","Pelagibac","hodobacter")

G$percentSynechocoBySample <- isoTax(G$o$metadata,"Synechococcales",G$o$counts)
G$percentSynechocoBySample <-
  G$percentSynechocoBySample[match(G$s$SampleID[G$sampleOrdering],
                                   names(G$percentSynechocoBySample))]


G$percentSynechocoBySpecies <-
  tapply(G$percentSynechocoBySample[match(G$s$SampleID, names(G$percentSynechocoBySample))],
         G$s$common,mean)
G$percentSynechocoBySpecies <-
  G$percentSynechocoBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                    names(G$percentSynechocoBySpecies))]

G$percentPelagibacBySample <- isoTax(G$o$metadata,"Pelagibac",G$o$counts)
G$percentPelagibacBySample <-
  G$percentPelagibacBySample[match(G$s$SampleID[G$sampleOrdering],
                                   names(G$percentPelagibacBySample))]
G$percentPelagibacBySpecies <-
  tapply(G$percentPelagibacBySample[match(G$s$SampleID, names(G$percentPelagibacBySample))],
         G$s$common,mean)
G$percentPelagibacBySpecies <-
  G$percentPelagibacBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                    names(G$percentPelagibacBySpecies))]

G$percentRoseobacterBySample <- isoTax(G$o$metadata,"hodobacter",G$o$counts)
G$percentRoseobacterBySample <-
  G$percentRoseobacterBySample[match(G$s$SampleID[G$sampleOrdering],
                                     names(G$percentRoseobacterBySample))]
G$percentRoseobacterBySpecies <-
  tapply(G$percentRoseobacterBySample[match(G$s$SampleID, names(G$percentRoseobacterBySample))],
         G$s$common,mean)
G$percentRoseobacterBySpecies <-
  G$percentRoseobacterBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                      names(G$percentRoseobacterBySpecies))]
















# Adonis Unifrac stats----

resetData("toUse")
resetUnifracData()
nperm <- 1e6


adonis(G$uur ~ phylum + class + order + family + genus + common,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ diet,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ phylum + diet,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ diet + phylum,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)



G$uurMammal <- as.dist(as.matrix(G$uur)[match(G$s$SampleID[G$s$class=="Mammalia"],attr(G$uur,"Labels")),
                                        match(G$s$SampleID[G$s$class=="Mammalia"],attr(G$uur,"Labels"))])


G$sUurMammal <- G$s[match(attr(G$uurMammal,"Labels"),G$s$SampleID),]
G$mammalOrders <- G$s %>% filter(class=="Mammalia") %>%
  distinct(order) %>% pull(order)



adonis(G$uurMammal ~ order == "Carnivora",G$s[match(attr(G$uurMammal,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ order == "Carnivora",G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ class == "Chondrichthyes",G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)


library(qiimer)
VW$uur_fp <- file.path(work_dir, "validation/whale454/beta_diversity", "unweighted_unifrac_otu_table.txt")
VW$uur <- read_qiime_distmat(VW$uur_fp)
attr(VW$uur,"Labels") <- VW$s$SampleID[match(attr(VW$uur,"Labels"),VW$s$old_id)]
VW$uur <- dist_subset(VW$uur, VW$s$SampleID)
adonis(VW$uur ~ (order == "Cetacea") + (family=="Delphinidae"|family=="Monodontidae"),
       VW$s[match(attr(VW$uur,"Labels"),VW$s$SampleID),], permutations = nperm)

VWI$uur_fp <- file.path(work_dir, "validation/whaleIllumina/beta_diversity", "unweighted_unifrac_otu_table.txt")
VWI$uur <- read_qiime_distmat(VWI$uur_fp)
attr(VWI$uur,"Labels") <- VWI$s$SampleID[match(attr(VWI$uur,"Labels"),VWI$s$old_id)]
VWI$uur <- dist_subset(VWI$uur, VWI$s$SampleID)
adonis(VWI$uur ~ (order == "Cetacea") + (family=="Delphinidae"|family=="Monodontidae"),
       VWI$s[match(attr(VWI$uur,"Labels"),VWI$s$SampleID),], permutations = nperm)









### Adonis on species phylogeny

resetData("toUse")
resetUnifracData()

G$phyMatrix <- cophenetic.phylo(G$newick)
rownames(G$phyMatrix) <-
  G$speciesTax$common[match(tolower(gsub("_"," ",rownames(G$phyMatrix))),
                            G$speciesTax$timetree_proxy)]
colnames(G$phyMatrix) <- rownames(G$phyMatrix)
G$phyMatrix <-
  G$phyMatrix[match(G$speciesTax$common[G$figure1SpeciesOrdering],rownames(G$phyMatrix)),
              match(G$speciesTax$common[G$figure1SpeciesOrdering],colnames(G$phyMatrix))]

G$phyDist <- as.dist(G$phyMatrix)

adonis(G$phyDist ~ diet, G$species, permutations = nperm)
adonis(G$phyDist ~ specificDiet, G$species, permutations = nperm)
adonis(G$phyDist ~ phylum, G$species, permutations = nperm)
adonis(G$phyDist ~ class, G$species, permutations = nperm)












# Fisher's tests----
# Fisher's test for unassigned OTUs vs usearch chimeras
fisher.test(as.matrix(rbind(c(145842,8814225),c(630689,22588000))))
# Fisher's test for singletons vs usearch chimeras
fisher.test(as.matrix(rbind(c(326406,16058371),c(450125,15343854))))
# Fisher's test for unassigned OTUs vs chimera_slayer chimeras
fisher.test(as.matrix(rbind(c(241,150944),c(1024,627523))))
# Fisher's test for singletons vs chimera_slayer chimeras
fisher.test(as.matrix(rbind(c(525,334908),c(740,443559))))
