
################################ Load the packages #############################

library("MEGENA") # For the Multiscale Embedded Gene Co-expression Network Analysis
library("readr") # This is a package to read large datasets quickly
library("tibble") # For data orgnaization
library("doParallel") # For parallel processing
library("foreach") # For parallel processing
library("org.Hs.eg.db") # For gene annotation


############################## Create the functions ############################


#~~~~~~~~~~~~~~~~~~~~~~ Function 1 -- Cleans residuals ~~~~~~~~~~~~~~~~~~~~~~~~#

clean_residuals <- function(df) {
  library("org.Hs.eg.db")
  library("dplyr")
  df1 = df %>% dplyr::rename("Ensemble_ID" = ...1)
    df1$Ensemble_ID <- sub("\\.\\d+", "", df1$Ensemble_ID)
    df1 = df1 %>% dplyr::rename(Gene_ID = Ensemble_ID)
      rownames(df1) = df1$Ensemble_ID
      genes = c(df1$Gene_ID)
        df2 <- select(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", 
                      columns = c("ENSEMBL", "SYMBOL"))
        df2 = df2[!duplicated(df2$ENSEMBL), ]
        df3 <- df2 %>% dplyr::rename(Gene_ID = ENSEMBL)
      df4 <- dplyr::left_join(df3, df1, by = "Gene_ID") %>%
      dplyr::select(SYMBOL, everything()) %>%
    dplyr::select(-c(Gene_ID)) %>% tidyr::drop_na() %>% as.data.frame()
    df5 <- df4 %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
  tibble::column_to_rownames("SYMBOL") 
df5 = df5 %>% head(80)  # I do this to shorten my own code
df6 = df5 %>% as.matrix()
return(df6)
}

nimh_matrix_race <-clean_residuals(residuals_nimh_race)

#~~~~~~~~~~~~~~~~~~~~~~~ Function 2 -- Runs the MEGENA ~~~~~~~~~~~~~~~~~~~~~~~~#

start_time <- Sys.time()
run_MEGENA <- function(nimh_matrix_race, n.cores = 12, doPar = TRUE, 
                       method = "pearson", FDR.cutoff = 0.05, 
                       module.pval = 0.05, hub.pval = 0.05, cor.perm = 10,
                       hub.perm = 100, annot.table = NULL, id.col = 1, 
                       symbol.col = 2) {
  # load libraries
  library("MEGENA")
  library("dplyr")
  library("doParallel")
  library("foreach")
  
  # Calculating correlations 
  ijw <- MEGENA::calculate.correlation(nimh_matrix_race, doPerm = cor.perm,
                                       output.corTable = FALSE, 
                                       output.permFDR = FALSE, method = method)
  
  # Calculating the Planar Filter Network (PFN)
  # register multiple cores if needed: note that set.parallel.backend() is deprecated. 
  run.par <- doPar & (getDoParWorkers() == 1) 
  if (run.par) {
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    # check how many workers are there
    cat(paste("number of cores to use:", getDoParWorkers(), "\n", sep = ""))
  }
  
  el_1 <- MEGENA::calculate.PFN(ijw[,1:3], doPar = doPar, 
                                num.cores = n.cores, 
                                keep.track = FALSE)
  el_1 <- na.omit(el_1)
  g <- igraph::graph.data.frame(el_1, directed = FALSE)
  g <- na.omit(g)
  
  # perform MCA clustering.
  MEGENA.output <- MEGENA::do.MEGENA(g, mod.pval = module.pval, hub.pval = hub.pval, 
                                     remove.unsig = TRUE,
                                     min.size = 10, max.size = vcount(g)/2, 
                                     doPar = doPar, num.cores = n.cores,
                                     n.perm = hub.perm, save.output = FALSE)
  
  # unregister cores as these are not needed anymore.
  if (getDoParWorkers() > 1) {
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  
  summary.output <- MEGENA.ModuleSummary(MEGENA.output, mod.pvalue = module.pval, 
                                         hub.pvalue = hub.pval,
                                         min.size = 10, max.size = vcount(g)/2, 
                                         annot.table = annot.table, 
                                         id.col = id.col, 
                                         symbol.col = symbol.col, 
                                         output.sig = TRUE)
  
  if (!is.null(annot.table)) {
    # update annotation to map to gene symbols
    V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name, 
                                                       annot.table[[id.col]])], 
                       V(g)$name, sep = "|")
    summary.output <- output[c("mapped.modules","module.table")]
    names(summary.output)[1] <- "modules"
  }
  result_list <- list(summary_output = summary.output, g_network = g)
  names(result_list) <- c("module_summary", "network_graph")
  
  end_time <- Sys.time()
  cat("Time taken to complete: ", end_time - start_time)
  return(result_list)
}

results <- run_MEGENA(nimh_matrix_race)


#~~~~~~~~~~~~~ Function 3 -- Runs functions 1 and 2 in parallel ~~~~~~~~~~~~~~~#

# set up parallel backend
cl <- makeCluster(6)
registerDoParallel(cl)

# define list of data frames
df_list <- list(residuals_nimh_race, residuals_ppm_race,  # Race
                residuals_nimh_dx,residuals_ppm_dx, # Dx
                residuals_nimh_interaction, residuals_ppm_interaction) # Interaction

start_time <- Sys.time()
    # apply clean_residuals and run_MEGENA to each data frame in parallel
    result_list <- foreach::foreach(df = df_list, 
                                    .packages = c("MEGENA", "dplyr", 
                                                  "doParallel", "foreach", 
                                                  "org.Hs.eg.db")) %dopar% {
      # apply clean_residuals
      cleaned_df <- clean_residuals(df)
      # apply run_MEGENA
      meg_output <- run_MEGENA(cleaned_df)
      # return result 
      list(cleaned_df = cleaned_df, meg_output = meg_output)
    }
  # Get the total run time
  end_time <- Sys.time()
  cat("Time taken to complete: ", end_time - start_time)
  # stop parallel backend
  stopCluster(cl)
  
# name each output in the result_list with its corresponding data frame name
names(result_list) <- c("residuals_nimh_race_megena", "residuals_ppm_race_megena", 
                          "residuals_nimh_dx_megena", "residuals_ppm_dx_megena",
                          "residuals_nimh_interaction_megena", "residuals_ppm_interaction_megena")


################################ Save the file #################################
saveRDS(result_list, file = "my_summary_output.RDS")





































