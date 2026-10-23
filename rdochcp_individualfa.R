# title: "Factor analysis for individuals (training model)"
# author: "Shaun Quah"
# output: "individual_tfmri_csv/output/", sub_id, "_FAoutput.rds"

# Load libraries
library(psych) 
library(lavaan)
library(dplyr)
library(openxlsx)
library(stringr)
library(GPArotation)

# Set pwd with individual_tfmri_csv folder before starting
setwd()

# List all files in the directory matching the pattern
pattern <- "_all_zstat1_gordon_reshaped.csv$"
files <- list.files(path = "individual_tfmri_csv/", pattern = pattern, full.names = TRUE)

# Extract subject IDs from file names
subject_ids <- str_extract(files, "\\d+")

#load hcp contrast subset
rdoc.contrastlist1 <- read.xlsx("individual_tfmri_csv/Contrasts_Domains.xlsx", cols = 4:5)

# Loop over each subject ID
for (sub_id in subject_ids) {
    filename <- paste0("individual_tfmri_csv/", sub_id, "_all_zstat1_gordon_reshaped.csv")
    
    data0 <- read.csv(filename)
    data0 <- data0[, order(names(data0))] # order alphabetically
    data0[] <- lapply(data0, as.numeric) # ensure all number is numeric and convert to matrix
    data0 <- apply(data0, 2, scale) # scale/standardize the data
    data0 <- head(data0, -5) # remove last 5 rows of the data to match HarvardOxford subcortical regions in original dataset
    
    # Subsets
    # Subset with variations, no average and deltas (24)
    columns_to_keep <- c("L1", "L2", "R1", "R2", "S1", "S2", "E1", "E2", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "M1", "M2", "M3", "M4", "M5", "M6", "G1", "G2")
    data1 <- data0[, columns_to_keep]

    
    perform_efa <- function(data) {
    tryCatch({
      pa <- suppressMessages(fa.parallel(data, fm="pa", fa="fa"))
      efa.fit <- fa(data, nfactors=pa$nfact, rotate="oblimin", fm="pa")
      loadings <- efa.fit$loadings[order(row.names(efa.fit$loadings)), ]
      loadings[loadings < .4 & loadings > -.4] <- 0
      return(loadings)
    }, error = function(e) {
      message("EFA Error: ", e)
        return(matrix(NA, nrow = ncol(data), ncol = 1)) # Return empty matrix if error occurs
    })
    }
    loadings1 <- perform_efa(data1)

    # CFA model based on EFA
    create_eDDcfamodel <- function(loadings) {
      # Extract variables with nonzero loadings from the loadings table
      nonzero_loadings <- loadings[rowSums(abs(loadings)) > 0, ]
      
      # Create a list to store factor definitions
      factor_definitions <- list()
      
      # Iterate through the columns (factors) with nonzero loadings
      for (i in 1:ncol(nonzero_loadings)) {
        # Get the factor name
        factor_name <- colnames(nonzero_loadings)[i]
        # Get the variable names with nonzero loadings for this factor
        variables <- rownames(nonzero_loadings)[nonzero_loadings[, i] != 0]
        
        # Only add the factor if it has variables assigned
        if (length(variables) > 0) {
          # Create the factor definition for this factor
          factor_definition <- paste(factor_name, "=~", paste(variables, collapse = " + "))
          # Add the factor definition to the list
          factor_definitions[[i]] <- factor_definition
        }
      }
      
      # Create the 'g' factor definition
      all_variables <- rownames(loadings)
      g_definition <- paste("g =~", paste(all_variables, collapse = " + "))
      
      # Remove any empty factor definitions
      factor_definitions <- factor_definitions[sapply(factor_definitions, nchar) > 0]
      
      # Combine the factor definitions into a single cfa model specification
      ecfa_fit_model <- c(g_definition, factor_definitions)
      ecfa_fit_model <- paste(ecfa_fit_model, collapse = "\n")
      
      return(ecfa_fit_model)
    }
    
    ecfa.fit.model1 <- create_eDDcfamodel(loadings1)

    # RDoC based on task description
    create_eRcfamodel <- function(rdoc.contrastlist) {
      rdoc.fit.model <- rdoc.contrastlist %>%
        group_by(Domains) %>%
        summarise(Contrasts = paste(Contrasts, collapse = " + ")) %>%
        mutate(FactorString = paste(Domains, "=~", Contrasts)) %>%
        pull(FactorString) %>%
        paste(collapse = "\n")
      
        return(rdoc.fit.model)
    }
    
    erdoc.fit.model1 <- create_eRcfamodel(rdoc.contrastlist1)
    
    # run DD CFA
    run_DDcfamodel <- function(cfa.fit.model, data) {
      tryCatch({
        cfa.fit <- cfa(cfa.fit.model, data, estimator = "MLR", std.lv = TRUE, orthogonal = TRUE, check.gradient = FALSE)
        
        cfa.fit.loadings <- inspect(cfa.fit,what="std")$lambda
        cfa.fit.loadings <- cfa.fit.loadings[order(row.names(cfa.fit.loadings)), ]
        cfa.fit.loadings <- cfa.fit.loadings[, order(colnames(cfa.fit.loadings))]
        
        cfa.fit.fscores <- lavPredict(cfa.fit, method = "Bartlett")
        cfa.fit.measures <- tryCatch({
          fitMeasures(cfa.fit, c("rmsea.robust", "cfi.robust", "tli.robust", "aic", "bic"))
        }, error = function(e) {
          NA
        })
        
        return(list(loadings = cfa.fit.loadings, fscores = cfa.fit.fscores, fit.measures = cfa.fit.measures))
      }, error = function(e) {
        # If an error occurs, return a list of NAs
        return(list(loadings = NA, fscores = NA, fit.measures = NA))
      })
    }
    
    edd.fit.output1 <- run_DDcfamodel(ecfa.fit.model1, data1)
    
    # run RDoC CFA
    run_Rcfamodel <- function(rdoc.fit.model, data) {
      tryCatch({
        rdoc.fit <- cfa(rdoc.fit.model, data, estimator = "MLR", std.lv = TRUE, check.gradient = FALSE)
        
        rdoc.fit.loadings <- inspect(rdoc.fit,what="std")$lambda
        rdoc.fit.loadings <- rdoc.fit.loadings[order(row.names(rdoc.fit.loadings)), ]
        rdoc.fit.loadings <- rdoc.fit.loadings[, order(colnames(rdoc.fit.loadings))]
        
        rdoc.fit.fscores <- lavPredict(rdoc.fit, method = "Bartlett")
        rdoc.fit.measures <- tryCatch({
          fitMeasures(rdoc.fit, c("rmsea.robust", "cfi.robust", "tli.robust", "aic", "bic"))
        }, error = function(e) {
          NA
        })
      
        return(list(loadings = rdoc.fit.loadings, fscores = rdoc.fit.fscores, fit.measures = rdoc.fit.measures))
      }, error = function(e) {
        # If an error occurs, return a list of NAs
        return(list(loadings = NA, fscores = NA, fit.measures = NA))
      })
    }
    
    erdoc.fit.output1 <- run_Rcfamodel(erdoc.fit.model1, data1)
    
    # Add folder to save in
    outputlist <- list(edd.fit.output1 = edd.fit.output1, erdoc.fit.output1 = erdoc.fit.output1)
    outputname <- paste0("individual_tfmri_csv/output/", sub_id, "_FAoutput.rds")
    saveRDS(outputlist, outputname)
}