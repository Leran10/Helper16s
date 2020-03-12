#' @title 16s pipeline plot helper --- plot prevelance
#'
#' @description helps generating the plots
#'
#' @param phyloseq object ps0, colorby,facetBy
#'
#' @return ggplots
#'
#' @examples plotAlpha(ps0,colorby = "Family",facetby = "Phylum")
#'
#' @export
#'
#'


plotPrevelance <- function(ps0,colorby = NULL,facetby = NULL){

  prevdf <- apply(X = phyloseq::otu_table(ps0),MARGIN = ifelse(phyloseq::taxa_are_rows(ps0), yes = 1, no = 2),FUN = function(x){sum(x > 0)})

  prevdf <- data.frame(Prevalence = prevdf,TotalAbundance = phyloseq::taxa_sums(ps0),phyloseq::tax_table(ps0))

  p.prev <- ggplot2::ggplot(prevdf,ggplot2::aes(TotalAbundance,Prevalence/phyloseq::nsamples(ps0),color = prevdf[[colorby]]))+
    ggplot2::geom_point(size = 3,alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0.05,alpha = 0.5,linetype = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::theme(legend.position = "none")+
    ggplot2::labs(title = paste0("Phylum Prevalence in All Samples\nColored by ",colorby),
                  caption = "Horizontal Line Placed at 3% Prevalence",
                  x = "Total Abundance",
                  y = "Prevalence")

  if (!is.null(facetby)){
    facet <- as.formula(paste0("~",facetby))
    p.prev <- p.prev + ggplot2::facet_wrap(facet)

  }


  return(p.prev)


}


