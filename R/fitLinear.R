#' @title 16s pipeline plot helper --- linear stats
#'
#' @description helps generating the linear test results
#'
#' @param phyloseq object
#'
#' @return table
#'
#' @examples fitLinear(ps0)
#'
#' @export
#'
#'
fitLinear <- function(ps0,x = NULL,z = NULL,Cohort = NULL,method = NULL){

  diversity <-estimate_richness(ps0, measures = c("Observed","Shannon"))
  row.names(diversity) <- c(sub('\\X',"",rownames(diversity)))
  df.richness <- merge(as.data.frame(sample_data(ps0)),diversity,by = "row.names")

  # observe
  p.observed <- ggplot(df.richness, aes(x = as.numeric(df.richness[[x]]), y = as.numeric(df.richness[["Observed"]]), color = Cohort)) +
    stat_smooth(method = "lm") +
    scale_y_log10() +
    scale_x_log10() +
    labs(y = "Richness", x = x, title = paste0("lm"," model for ",x = x)) +
    geom_jitter(size = 2, alpha = 0.5, width = 0.2)


  p.shannon <- ggplot(df.richness, aes(x = as.numeric(df.richness[[x]]), y = as.numeric(df.richness[["Shannon"]]), color = Cohort)) +
    stat_smooth(method = "lm") +
    scale_x_log10() +
    labs(y = "Shannon", x = x, title = paste0("lm"," model for ",x)) +
    geom_jitter(size = 2, alpha = 0.5, width = 0.2)


  if (!is.null(z)){

    print("===========================================Observed===============================================")
    cat("\n\n")

    print("===========correlation===============")
    print(cor.test(df.richness[["Observed"]], df.richness[[x]], method = "spearman"))

    print("============if Hep.B.M12 can predict observed============")
    formula <- as.formula(paste0("Observed","~",x,"*",z))
    print(summary(aov(lm.fit <- lm(formula, data = df.richness))))

    print("=============delete the intersection==============")
    formula <- as.formula(paste0("Observed","~",x,"*",z, "- 1"))
    print(summary(lm(formula, data = df.richness)))


    print("=============with intersection============")
    formula <- as.formula(paste0("Observed","~",x,"*",z))
    fit1 <- aov(formula, data = df.richness)
    print(summary(fit1))


    print("============No intersection===============")
    formula <- as.formula(paste0("Observed","~",x,"+",z))
    fit2 <- aov(formula, data = df.richness)
    print(summary(fit2))

    print("============compare model1 and model2===============")
    print(anova(fit1,fit2))


    print("===========================================Shannon===============================================")
    cat("\n\n")

    print("===========correlation===============")
    print(cor.test(df.richness[["Shannon"]], df.richness[[x]], method = "spearman"))

    print("============if Hep.B.M12 can predict Shannon============")
    formula <- as.formula(paste0("Shannon","~",x,"*",z))
    print(summary(aov(lm.fit <- lm(formula, data = df.richness))))

    print("=============delete the intersection==============")
    formula <- as.formula(paste0("Shannon","~",x,"*",z, "- 1"))
    print(summary(lm(formula, data = df.richness)))


    print("=============with intersection============")
    formula <- as.formula(paste0("Shannon","~",x,"*",z))
    fit1 <- aov(formula, data = df.richness)
    print(summary(fit1))


    print("============No intersection===============")
    formula <- as.formula(paste0("Shannon","~",x,"+",z))
    fit2 <- aov(formula, data = df.richness)
    print(summary(fit2))

    print("============compare model1 and model2===============")
    print(anova(fit1,fit2))

  }else if (is.null(z)){

    print("============================================Observed=========================================")
    cat("\n\n")
    #Observed
    print("===============correlation===============")
    print(cor.test(df.richness[["Observed"]], df.richness[[x]], method = "spearman"))

    print("============if x can predict observed============")
    formula <- as.formula(paste0("Observed","~", x))
    print(summary(lm(formula, data = df.richness)))


    print("============================================Shannon==========================================")
    cat("\n\n")
    print("================correlation ===============")
    print(cor.test(df.richness[["Shannon"]], df.richness[[x]], method = "spearman"))

    print("============if x can predict Shannon============")
    formula <- as.formula(paste0("Shannon","~", x))
    print(summary(lm(formula, data = df.richness)))
  }

  ggarrange(p.observed,p.shannon,ncol = 2,common.legend = TRUE)


}
