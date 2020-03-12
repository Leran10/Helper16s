#' @title 16s helper functions --- assess data
#'
#' @description help generate 16s related plots
#'
#' @param phyloseq object,category,hline
#'
#' @return ggplot plot
#'
#' @examples assessData(ps0,"INFECTED",100000)
#'
#' @export

library(tidyverse)
library(phyloseq)
library(ggplot2)
library(plyr)
library(dplyr)
library(xlsx)
library(stringr)
library(ggpubr)
library(vegan)

assessData <- function(ps0,category = NULL, hline = NULL){

  df <- merge(data.frame(sample_sums(physeq)),sample_data(physeq),by = "row.names") %>%
    mutate(Sample = "Sample")


  p <- ggboxplot(df,x = category,y = "sample_sums.physeq.",color = category)+
    geom_point()+
    geom_jitter(width = 0.2) +
    labs(y="ASV",x = "")

  if (!is.null(hline)){
    p <- p + geom_hline(yintercept=hline, linetype="dashed", color = "black")
  }

  return(p)

}


