#' @title 16s pipeline plot helper --- plot beta
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

# plot beta

plotBeta <- function(ps0,distance = NULL,factor = NULL,binaryGroup = FALSE,title = NULL){

  if(any(is.na(sample_data(ps0)[[factor]]))){
    stop("please clear up all the NA values in the category you are interested in.")
  }

  if(length(distance) == 0){
    stop("please specify the name of distance metric you want to use.")
  }

  if(is.na(factor)){
    stop("please specify the factor you are interested in.")
  }

  phyloseq::sample_data(ps0)[[factor]] <- as.factor(phyloseq::sample_data(ps0)[[factor]])

  if (binaryGroup){
    # Function to run adonis test on a physeq object and a variable from metadata
    doadonis<- function(ps0, distance = distance, factor = factor) {
      bdist <- phyloseq::distance(ps0, distance)
      col <- as(sample_data(ps0), "data.frame")[ ,factor]

      # Adonis test
      adonis.bdist <- vegan::adonis(bdist ~ col)
      p <- round(as.data.frame(adonis.bdist[1])$aov.tab.Pr..F.[1],digits = 3)
      r2 <- round(as.data.frame(adonis.bdist[1])$aov.tab.R2[1],digits = 2)

      res <- c(p,r2)
      return(res)

    }

    plist <- list()
    for (i in (1:length(distance))){


      res <- doadonis(ps0,distance = distance[i],factor = factor)

      ps <- phyloseq::ordinate(ps0,method = "PCoA",distance = distance[i])

      p <- phyloseq::plot_ordination(ps0,ps,color = factor) +
        ggplot2::geom_point(alpha = 0.8,size = 3) +
        ggplot2::stat_ellipse(type = "norm") +
        ggplot2::scale_color_brewer(palette = "Dark2") +
        ggplot2::ggtitle(paste0(distance[i]),subtitle = paste0("Adonis: p = ",as.character(res[1]),",R2 = ", " ",as.character(res[2])))

      plist[[i]] <- p

    }


  }else if(!binaryGroup){

    # Function to run adonis test on a physeq object and a variable from metadata
    doanosim<- function(ps0, distance = distance, factor = factor) {

      bdist <- phyloseq::distance(ps0, distance)
      print("test1")
      col <- as(sample_data(ps0), "data.frame")[ ,factor]
      print("test2")
      attach(sample_data(ps0))


      # Anosim test
      anosim.dist <- vegan::anosim(bdist, col, permutations = 999)
      print("test3")
      p <- round(anosim.dist$signif,digits = 3)
      r2 <- round(anosim.dist$statistic,,digits = 3)

      res <- c(p,r2)
      return(res)

    }

    plist <- list()
    for (i in (1:length(distance))){


      res <- doanosim(ps0,distance = distance[i],factor = factor)

      ps <- phyloseq::ordinate(ps0,method = "PCoA",distance = distance[i])

      p <- phyloseq::plot_ordination(ps0,ps,color = factor) +
        ggplot2::geom_point(alpha = 0.8,size = 3) +
        ggplot2::stat_ellipse(type = "norm") +
        ggplot2::scale_color_brewer(palette = "Dark2") +
        ggplot2::ggtitle(paste0(distance[i]),subtitle = paste0("Anosim: p = ",as.character(res[1]),",R2 = ", " ",as.character(res[2])))

      plist[[i]] <- p

    }


  }
  if (length(plist) == 1){
    ggp <- ggarrange(plist[[1]])
  }else if(length(plist) == 2){
    ggp <- ggarrange(plist[[1]],plist[[2]],ncols = 2,common.legend = TRUE)
  }else if (length(plist) == 3){
    ggp <- ggarrange(plist[[1]],plist[[2]],plist[[3]],ncols = 2,nrows = 2,common.legend = TRUE)
  }

  return(annotate_figure(ggp,
                         top = text_grob(paste0(title))))


}



