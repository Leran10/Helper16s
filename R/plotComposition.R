#' @title
#'
#' @description
#'
#' @param
#'
#' @return
#'
#' @examples
#'
#' @export
#'


#source("makeColor.R")
plotComposition <- function(ps0,category = NULL, facet = NULL){

  print("before remove")
  print(ps0)


  cat("remove mitochondria in Family level\nremove Chloroplast from Class level\nremove Cyanobacteria/Chloroplast from Phylum level")
  ps0 <- ps0 %>%
    subset_taxa(
      Family  != "mitochondria" &
        Class   != "Chloroplast" &
        Phylum != "Cyanobacteria/Chloroplast"
    )

  print("after remove")
  print(ps0)

  sample_data(ps0)[[category]] <- as.factor(sample_data(ps0)[[category]])


  ps0.phylum <- ps0 %>% tax_glom(taxrank = "Phylum")
  topnames <- names(sort(taxa_sums(ps0.phylum), TRUE)[1:8])
  ps0.top.phylum <- prune_taxa(topnames, ps0.phylum)


  df.phylum <-  ps0.top.phylum %>%
    transform_sample_counts(function(x) {x/sum(x)}) %>%
    psmelt() %>%
    arrange(Phylum)


  p <- ggplot(df.phylum, aes(df.phylum[[category]], y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", width = 0.9, position = "fill") +
    scale_color_brewer(palette = "Dark2")+
    labs(title = "phylum level compocition plot", y="Relative Abundance", x="") +
    theme(legend.position = "right",
          legend.title = element_text(size=5, face="bold"),
          axis.text.x = element_text(angle = 90, hjust = 1))

  if (!is.null(facet)){
    facetFormula <- as.formula(paste0("~",facet))
    p <- p + facet_wrap(facetFormula,scales= "free_x", nrow=1)

  }


  return(p)

}

