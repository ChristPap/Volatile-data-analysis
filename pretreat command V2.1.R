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
  # scale = "auto" / "pareto" / "none"

  pretreat <- function(x, transf, center, scale) {
    
    # (I) Name the analysis
    ## Set the Transformation type
    if(transf == "sqrt") {
      a <- "SqRt Transf"
    } else if(transf == "log2") {
      a <- "log2 Transf"
    } else if(transf == "log10") {
      a <- "log10 Transf"
    } else if(transf == "none") {
      a <- "No Transf"
    }
    ## Set the Centering type
    if(center == TRUE) {
      b <- "Centering"
    } else if(center == FALSE) {
      b <- "No centering"
    } 
    ## Set the Scaling type
    if(scale == "auto") {
      c <- "Autoscaling"
    } else if(scale == "pareto") {
      c <- "Paretoscaling"
    } else if(scale == "none") {
      c <- "No Scaling"
    }
    ## Prepare the analysis name
    mytrasf <- paste(a,b,c, sep = ",")
    
    # (II) Perform data transformation
    if(transf == "none") {
      transf.dt <- x
    } else if (transf == "sqrt") {
      transf.dt <- sqrt(x)
    } else if (transf == "log2") {
      # In order to perform log2 transformation we first need to deal with negative and 0 values. For this reason we will add to the values of our dataset the 1/10th of the minimum value of the table
      min.val <- min(abs(x[x!=0]))/10
      dt_prelog <- x + min.val
      transf.dt <- log2(dt_prelog)
    } else if (transf == "log10") {
      # In order to perform log10 transformation we first need to deal with negative and 0 values. For this reason we will add to the values of our dataset the 1/10th of the minimum value of the table
      min.val <- min(abs(x[x!=0]))/10
      dt_prelog <- x + min.val
      transf.dt <- log10(dt_prelog)
    }
    
    # (III) Perform data scaling with(out) centering
    if(center == TRUE & scale == "auto" ) {
      # Autoscaling, also known as standardization, is when mean-centered values (values from each of which the mean is subtracted) are divided with their standard deviation
      scaled.transf.dt <- scale(transf.dt, center = T, scale = T)
    } else if(center == TRUE & scale == "pareto") {
      # Paretoscaling is when mean-centered values are divided by the square root of the standard deviation
      transf.dt <- apply(transf.dt, 2, function(x) x - mean(x) )
      scaled.transf.dt <- apply(transf.dt, 2, function(x) ( x - mean(x) ) / sqrt(sd(x)) )
    } else if(center == TRUE & scale == "none") { 
      scaled.transf.dt <- scale(transf.dt, center = T, scale = F)
    } else if(center == FALSE & scale == "auto") {
      # Autoscaling, also known as standardization, is when mean-centered values (values from each of which the mean is subtracted) are divided with their standard deviation
      scaled.transf.dt <- scale(transf.dt, center = F, scale = T)
    } else if(center == FALSE & scale == "pareto") {
      # Paretoscaling is when mean-centered values are divided by the square root of the standard deviation
      scaled.transf.dt <- apply(transf.dt, 2, function(x) ( x - mean(x) ) / sqrt(sd(x)) )
    } else if (center == FALSE & scale == "none") {
      scaled.transf.dt <- transf.dt
    }
    
    output <- list(mytrasf = mytrasf, transf_mydt = scaled.transf.dt)    
    return(output)
    
  }
