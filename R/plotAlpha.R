#' @title 16s pipeline plot helper --- plot alpha
#'
#' @description helps generating the plots
#'
#' @param phyloseq object, x, binaryGroup,paired,parametric
#'
#' @return ggplots
#'
#' @examples plotAlpha(ps0,x = "Extraction.Method",binaryGroup = TRUE,paired = TRUE,parametric = TRUE)
#'
#' @export
#'
library(ggplot2)
library(phyloseq)
library(ggpubr)

plotAlpha <- function(ps0,
                      x = NULL,
                      facet = list(),
                      paired = TRUE,
                      parametric = TRUE){

  phyloseq::sample_data(ps0)[[x]] <- as.factor(phyloseq::sample_data(ps0)[[x]])

  ps.diversity <- phyloseq::estimate_richness(ps0,measures = c("Shannon","Observed"))
  rownames(ps.diversity) <- c(sub("\\X","",rownames(ps.diversity)))
  ps.meta <- as.data.frame(phyloseq::sample_data(ps0))

  if (!identical(rownames(ps.diversity),rownames(ps.meta))) {
    stop("Error: the rownames of diversity table and meta table should be consistent" )
  }
  ps.meta.rich <- merge(ps.diversity,ps.meta, by = "row.names")




  # unpaired
  if (!paired){
    print("non-paired")

    # non-paramatric
    if (!parametric){     # use wilcox.test(default one)
      print("non-parametric")

      length <-  length(levels(as.factor(ps.meta.rich[[x]])))
      my_comparisons <- list()

      count = 1
      for (i in 1 : length){

        j = i+1
        while (j<= length){
          temp <- c(levels(as.factor(ps.meta.rich[[x]]))[i],levels(as.factor(ps.meta.rich[[x]]))[j])
          print(temp)
          my_comparisons[[count]] <- temp
          count <- count +1
          j <- j+1
        }

      }

      print(my_comparisons)


      p.observed <- ggplot2::ggplot(ps.meta.rich, aes(x = ps.meta.rich[[x]], y = Observed,group = ps.meta.rich[[x]])) +
        geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
        #facet_wrap(~Cohort~timeLabel)+
        ggpubr::stat_compare_means(comparisons = my_comparisons,label.y = max(ps.meta.rich$Observed)+10) +
        labs(y = "Observed Richness", x = "") +
        ylim(min(ps.meta.rich$Observed),max(ps.meta.rich$Observed)+60)+
        theme(legend.position = "bottom")

      p.shannon <- ggplot2::ggplot(ps.meta.rich, aes(x = ps.meta.rich[[x]], y = Shannon,group = ps.meta.rich[[x]])) +
        #facet_wrap(~Cohort~timeLabel)+
        geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
        ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = max(ps.meta.rich$Shannon)+0.5) +
        labs(y = "Shannon Index", x = "") +
        ylim(min(ps.meta.rich$Shannon),max(ps.meta.rich$Shannon)+1)+
        theme(legend.position = "bottom")

      if(length(facet) > 1){

        formula <- as.formula(paste0("~",facet[1],"~",facet[2]))

      }else if(length(facet) == 1){
        formula <- as.formula(paste0("~",facet[1]))

      }

      p.observed <- p.observed + ggplot2::facet_wrap(formula)
      p.shannon <-  p.shannon + ggplot2::facet_wrap(formula)




    }else if(paramatric){     #use t.test

      print("parametric")

      #modified this
      my_comparisons <- list()

      length <-  length(levels(as.factor(ps.meta.rich[[x]])))

      count = 1
      for (i in 1 : length){

        j = i+1
        while (j<= length){
          temp <- c(levels(as.factor(ps.meta.rich[[x]]))[i],levels(as.factor(ps.meta.rich[[x]]))[j])
          print(temp)
          my_comparisons[[count]] <- temp
          count <- count +1
          j <- j+1
        }

      }

      print(my_comparisons)


      p.observed <- ggplot2::ggplot(ps.meta.rich, aes(x = ps.meta.rich[[x]], y = Observed, group = ps.meta.rich[[x]])) +
        geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
        #facet_wrap(~Cohort~timeLabel)+
        ggpubr::stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(ps.meta.rich$Observed)+10) +
        labs(y = "Observed Richness", x = "") +
        ylim(min(ps.meta.rich$Observed),max(ps.meta.rich$Observed)+60)+
        theme(legend.position = "bottom")

      p.shannon <- ggplot2::ggplot(ps.meta.rich, aes(x = ps.meta.rich[[x]], y = Shannon, group = ps.meta.rich[[x]])) +
        #facet_wrap(~Cohort~timeLabel)+
        geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
        ggpubr::stat_compare_means(method = "t.test",comparisons = my_comparisons, label.y = max(ps.meta.rich$Shannon)+0.5) +
        labs(y = "Shannon Index", x = "") +
        ylim(min(ps.meta.rich$Shannon),max(ps.meta.rich$Shannon)+1)+
        theme(legend.position = "bottom")

      if(length(facet) > 1){

        formula <- as.formula(paste0("~",facet[1],"~",facet[2]))

      }else if(length(facet) == 1){
        formula <- as.formula(paste0("~",facet[1]))

      }

      p.observed <- p.observed + ggplot2::facet_wrap(formula)
      p.shannon <-  p.shannon + ggplot2::facet_wrap(formula)



    }



  }else if(paired){


    # non-paramatric
    if (!parametric){  #use wilcox.test (default one)


      my_comparisons <- list()

      length <-  length(levels(as.factor(ps.meta.rich[[x]])))

      count = 1
      for (i in 1 : length){

        j = i+1
        while (j<= length){
          temp <- c(levels(as.factor(ps.meta.rich[[x]]))[i],levels(as.factor(ps.meta.rich[[x]]))[j])
          print(temp)
          my_comparisons[[count]] <- temp
          count <- count +1
          j <- j+1
        }

      }

      print(my_comparisons)


      p.observed <- ggplot2::ggplot(ps.meta.rich, aes(x = ps.meta.rich[[x]], y = Observed, group = ps.meta.rich[[x]])) +
        geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
        #facet_wrap(~Cohort~timeLabel)+
        ggpubr::stat_compare_means(comparisons = my_comparisons,label.y = max(ps.meta.rich$Observed)+10,paired = TRUE) +
        labs(y = "Observed Richness", x = "") +
        ylim(min(ps.meta.rich$Observed),max(ps.meta.rich$Observed)+60)+
        theme(legend.position = "bottom")

      p.shannon <- ggplot2::ggplot(ps.meta.rich, aes(x = ps.meta.rich[[x]], y = Shannon, group = ps.meta.rich[[x]])) +
        #facet_wrap(~Cohort~timeLabel)+
        geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
        ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = max(ps.meta.rich$Shannon)+0.5,paired = TRUE) +
        labs(y = "Shannon Index", x = "") +
        ylim(min(ps.meta.rich$Shannon),max(ps.meta.rich$Shannon)+1)+
        theme(legend.position = "bottom")

      if(length(facet) > 1){

        formula <- as.formula(paste0("~",facet[1],"~",facet[2]))

      }else if(length(facet) == 1){
        formula <- as.formula(paste0("~",facet[1]))

      }

      p.observed <- p.observed + ggplot2::facet_wrap(formula)
      p.shannon <-  p.shannon + ggplot2::facet_wrap(formula)




    }else if(paramatric){   # use t.test

      my_comparisons <- list()

      length <-  length(levels(as.factor(ps.meta.rich[[x]])))

      count = 1
      for (i in 1 : length){

        j = i+1
        while (j<= length){
          temp <- c(levels(as.factor(ps.meta.rich[[x]]))[i],levels(as.factor(ps.meta.rich[[x]]))[j])
          print(temp)
          my_comparisons[[count]] <- temp
          count <- count +1
          j <- j+1
        }

      }

      print(my_comparisons)


      p.observed <- ggplot2::ggplot(ps.meta.rich, aes(x = ps.meta.rich[[x]], y = Observed, group = ps.meta.rich[[x]])) +
        geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
        #facet_wrap(~Cohort~timeLabel)+
        ggpubr::stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(ps.meta.rich$Observed)+10,paired = TRUE) +
        labs(y = "Observed Richness", x = "") +
        ylim(min(ps.meta.rich$Observed),max(ps.meta.rich$Observed)+60)+
        theme(legend.position = "bottom")

      p.shannon <- ggplot2::ggplot(ps.meta.rich, aes(x = ps.meta.rich[[x]], y = Shannon, group = ps.meta.rich[[x]])) +
        #facet_wrap(~Cohort~timeLabel)+
        geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
        ggpubr::stat_compare_means(method = "t.test",comparisons = my_comparisons, label.y = max(ps.meta.rich$Shannon)+0.5,paired = TRUE) +
        labs(y = "Shannon Index", x = "") +
        ylim(min(ps.meta.rich$Shannon),max(ps.meta.rich$Shannon)+1)+
        theme(legend.position = "bottom")

      if(length(facet) > 1){

        formula <- as.formula(paste0("~",facet[1],"~",facet[2]))

      }else if(length(facet) == 1){
        formula <- as.formula(paste0("~",facet[1]))

      }

      p.observed <- p.observed + ggplot2::facet_wrap(formula)
      p.shannon <-  p.shannon + ggplot2::facet_wrap(formula)




    }


  }



  return(ggarrange(p.observed,p.shannon,nrow = 2,ncol = 1))

}
