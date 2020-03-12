#' @title 16s pipeline plot helper --- plot volcano
#'
#' @description helps generating the plots
#'
#' @param phyloseq object
#'
#' @return ggplots
#'
#' @examples plotBeta(ps0,distance = "unifrac,category = "Extraction.Method")
#'
#' @export
#'
#'
#'

plotVolcano <- function(ps0,ps0.subset,factorOfInterest = NULL,truthFacet = NULL,alpha = NULL,baseMean = NULL,title = NULL,output = NULL){

  ##count table replacement function
  replace_counts = function(ps0, dds) {

    dds_counts = DESeq2::counts(dds, normalized = TRUE)
    if (!identical(phyloseq::taxa_names(ps0), rownames(dds_counts))) {
      stop("OTU ids don't match")
    }


    phyloseq::otu_table(ps0) = phyloseq::otu_table(dds_counts, taxa_are_rows = TRUE)
    return(ps0)

  }




  # Make deseq ready object
  formula <- as.formula(paste0("~",factorOfInterest))
  ds.all <- phyloseq::phyloseq_to_deseq2(ps0, formula) # Can be any variable


  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }

  # get the table with all the columns
  geoMeans <- apply(DESeq2::counts(ds.all), 1, gm_mean)
  ds.all <- DESeq2::estimateSizeFactors(ds.all, geoMeans = geoMeans)
  dds.all <- DESeq2::DESeq(ds.all, fitType = "local")
  rlog.all <- replace_counts(ps0, dds.all)
  rlog.all <- phyloseq::psmelt(rlog.all)
  rlog.all <- dplyr::rename(rlog.all, "ASV" = "OTU")


  # convert phyloseq to deseq2 table
  ds.subset <- phyloseq::phyloseq_to_deseq2(ps0.subset,formula)
  #compute geomeans
  geoMeans <- apply(DESeq2::counts(ds.subset),1,gm_mean)
  #estimate size factors
  ds.subset <- DESeq2::estimateSizeFactors(ds.subset,geoMeans = geoMeans)
  # run deseq
  dds.subset <- DESeq2::DESeq(ds.subset)
  #get deseq result
  res.dds.subset <- results(dds.subset,cooksCutoff = FALSE)

  res.dds.subset <- as.data.frame(res.dds.subset[ which(res.dds.subset$baseMean > baseMean), ])


  #tabulate the result table
  restable.dds.subset <- data.table(as(res.dds.subset,"data.frame"),keep.rownames = TRUE)
  setnames(restable.dds.subset,"rn","OTU")
  setkeyv(restable.dds.subset,"OTU")



  #extract taxa tables
  taxa.subset <- data.table(data.frame(as(phyloseq::tax_table(ps0.subset),"matrix")),keep.rownames = TRUE)
  setnames(taxa.subset,"rn","OTU")
  setkeyv(taxa.subset,"OTU")

  restable.dds.subset <- as.data.frame(restable.dds.subset)
  taxa.subset <- as.data.frame(taxa.subset)

  rownames(restable.dds.subset) <- restable.dds.subset$OTU
  rownames(taxa.subset) <- taxa.subset$OTU

  restable_taxa.subset <- left_join(restable.dds.subset,taxa.subset)


  # set colors
  l <- unique(as.factor(cleaned.table.subset$Phylum))
  color <- colorRampPalette(brewer.pal( 9 , "Set1" ))(length(l))
  names(color) <- l



  # clean the table
  cleaned.table.subset <- restable_taxa.subset %>%
    filter(padj != "NA") %>%
    mutate(significant = padj < alpha)


  outdf <- subset(cleaned.table.subset,padj < 0.05)
  write.csv(outdf, output)

  print("start to plot volcano plot:")
  # valcano plot
  p.volcano <- ggplot(data = cleaned.table.subset,
                      aes(x = log2FoldChange,
                          y = -log10(padj),
                          color = Phylum,
                          label1 = Family,
                          label2 = Genus,
                          label3 = Phylum)) +
    scale_fill_manual(values = color)+
    geom_point(data = subset(cleaned.table.subset,cleaned.table.subset$significant == FALSE),
               color = "grey") +
    geom_point(data = subset(cleaned.table.subset,cleaned.table.subset$significant == TRUE),
               aes(color = Phylum,
                   size = baseMean))+
    geom_vline(xintercept = 0,lty = 2)+
    geom_hline(yintercept = -log10(0.05))+
    geom_text_repel(data = subset(cleaned.table.subset,cleaned.table.subset$significant == TRUE),aes(label=Family)) +
    labs(title = "",subtitle = paste0("alpha = ",alpha,",baseMean = ",baseMean)) +
    theme_bw() +
    theme(axis.title = element_text(size=12),
          axis.text = element_text(size=12),
          axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
          legend.text = element_text(size=8),
          legend.position = "bottom")




  print("start to prepare data for ground truth plot:")
  # ground truth
  # Create merge taxa table
  tax <- as.data.frame(tax_table(ps0))
  tax <-  as.tibble(rownames_to_column(tax, var = "ASV"))


  df.subset <- as.data.frame(res.dds.subset[which(res.dds.subset$padj < alpha & res.dds.subset$baseMean >baseMean), ])
  df.subset <- rownames_to_column(df.subset, var = "ASV")
  df.subset <- left_join(df.subset, tax, var = "ASV")

  #plot
  theme_set(theme_bw(base_size =8,
                     base_family = "Arial"))

  df.rlog.subset <- inner_join(df.subset, rlog.all, by = c("ASV"))
  formula.truth <- as.formula(paste0(truthFacet,"~ASV"))


  print("start to plot ground truth plot:")
  # Boxplots
  p.box <- ggboxplot(df.rlog.subset, x = factorOfInterest, y="Abundance", color = "Phylum.y",outlier.shape = NA) +
    scale_fill_manual("Phylum",values = color)+
    geom_jitter(alpha = 0.7, width = 0.2, size = 1)+
    scale_y_log10() +
    facet_wrap(formula.truth)+
    labs(title = "",x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    stat_compare_means(label = "p.signif", hide.ns = TRUE) +
    labs(color='Phylum')


  annotate_figure(ggarrange(p.volcano,p.box,ncol = 2,
                            top = text_grob(title)))


}

