


#' @title Download RNAseq data from TCGA
#' @description
#' This helper function downloads RNAseq data for the indicated dataset.
#' The raw data will be stored in a directory "raw.data.dir" and an rda version will be saved in 'rda.dir'
#'
#' @param dataset.name the name of the dataset to download e.g TCGA-GBM, TCGA-HNSC
#' @param raw.data.dir the directory to save raw data
#' @param rda.dir the directory to save the processed R compatiable data to
#' @return counts table of raw RNAseq values
#' @export
#' @examples
#' download.TCGA.RNAseq("TCGA-BRCA")

download.TCGA.RNAseq <- function(dataset.name,
                                 raw.data.dir = "GDCdata/",
                                 rda.dir = "data/"){

  # create directories
  if(!dir.exists(raw.data.dir)){
    dir.create(raw.data.dir, recursive = TRUE)
  }

  if(!dir.exists(rda.dir)){
    dir.create(rda.dir, recursive = TRUE)
  }


  print(paste0("downloading dataset ", dataset.name))


  query.raw.expression.hg38 <- TCGAbiolinks::GDCquery(project = dataset.name,
                                                      data.category = "Transcriptome Profiling",
                                                      workflow.type = "STAR - Counts",
                                                      data.type = "Gene Expression Quantification",
                                                      experimental.strategy = "RNA-Seq")

  TCGAbiolinks::GDCdownload(query.raw.expression.hg38,
                            files.per.chunk = 6,
                            directory = raw.data.dir)

  TCGA.raw.counts <- TCGAbiolinks::GDCprepare(query = query.raw.expression.hg38,
                                              summarizedExperiment = FALSE,
                                              save = TRUE,
                                              save.filename = paste0(rda.dir, dataset.name, "_Transcriptome_Profiling_raw_counts.rda"))

  return(TCGA.raw.counts)

}


#' @title Load data from an rda file
#' @description
#' load rda file that is saved from a download.TCGA.RNAseq function call
#' @param path.to.rda path to rda file
#' @return TCGA.raw.counts
#' @examples
#' load.rda("/home/user/Desktop/TCGA.raw.counts.rda")

load.rda <- function(path.to.rda){

  load(file = path.to.rda)
  TCGA.raw.counts <- data

  return(TCGA.raw.counts)

}








#' @title Preprocess TCGA RNAseq data
#' @description
#' Process the output of download.TCGA.RNAseq function by filtering for protein coding genes only, normalising by GC content, and quantile filtering for low expressed genes
#' @param input.data output from download.TCGA.RNAseq function
#' @param dataset.name character string containing the name of the dataset
#' @param quantile.filt logical indicating whether to filter the data by quantile
#' @param save.normed.data logical indicating whether to save the normalised data to file
#' @param output.dir character string containing the directory to save the normalised data to
#' @return  data.frame containing the normalised data
#' @examples
#' raw.data <- download.TCGA.RNAseq("TCGA-BRCA")
#' or
#' raw.data <- load.rda("path.to.rda.file")
#'
#' normed.data <- preprocess.TCGA(input.data = raw.data, dataset.name = "TCGA_BRCA")
#'
preprocess.TCGA <- function(input.data,
                            dataset.name = NULL,
                            quantile.filt = TRUE,
                            save.normed.data = TRUE,
                            output.dir = "data/"){


  # check that data.set name is defined if save.normed.data is set to TRUE
  if(save.normed.data & is.null(dataset.name)){
    stop("save.normed.data is set to TRUE but dataset.name is NULL, please define a dataset.name as this will be used in saving the file")
  }



  # filtering dataset
  print("Filtering dataset for protein coding etc.")

  filt.data <- input.data %>%
    dplyr::filter(.data$gene_type == "protein_coding") %>%
    dplyr::mutate(gene_sum = rowSums(dplyr::select(., -c(gene_id, gene_name, gene_type)))) %>%
    dplyr::arrange(dplyr::desc(.data$gene_sum)) %>%
    dplyr::distinct(.data$gene_name, .keep_all = TRUE) %>%
    tibble::column_to_rownames(var = "gene_name") %>%
    dplyr::select(-c(gene_id, gene_type, gene_sum)) %>%
    dplyr::select(tidyselect::matches("tpm_unstranded_.*"))

  # note that dataset contains unstranded, stranded_first, stranded_second, tpm_unstranded, fpkm_unstranded, fpkm_uq_unstranded counts for each sample
  # use tpm as tpm is best for comparing expression between samples

  print("Normalising by GC content")

  normed.data <- TCGAbiolinks::TCGAanalyze_Normalization(tabDF = filt.data,
                                                         geneInfo = TCGAbiolinks::geneInfo, # geneInfo or geneInfoHT
                                                         method = "gcContent")


  if(quantile.filt){
    # filter genes to reduce low expressed genes and false positive results
    # quantile filter of genes
    print("Quantile Filtering dataset")

    normed.data <- TCGAbiolinks::TCGAanalyze_Filtering(tabDF = normed.data,
                                                       method = "quantile",
                                                       qnt.cut =  0.25)

  }




  if(save.normed.data){
    print("Writing normalised & Filtered data to file")

    # Writing Normalised data
    utils::write.table(normed.data, paste0(output.dir, "Normalised_filtered_", dataset.name, "_RNAseq_data.txt"),
                sep = "\t", row.names = T, quote = F)


  }

  return(normed.data)

}





#' @title threshold.goi
#' @description This function takes a gene of interest (goi) and thresholds the expression data to identify patients with high and low expression of the goi.
#' @param input.data A matrix of gene expression data.
#' @param goi A gene of interest.
#' @return A list of patients with high and low expression of the goi.
#' @examples
#' threshold.goi(input.data = expression.data, goi = "NKG7")

threshold.goi <- function(input.data, goi){


if(length(goi) == 1){

# extract goi expression data
  goi.expression <- input.data[goi, ]

  # get some expression summary statistics
  summary.value <- summary(goi.expression)

  # define upper and lower quartiles of expression
  upper.quartile <- floor(summary.value[5])
  lower.quartile <- ceiling(summary.value[2])

# use threshold to define high expressing patients
  High.exp.logic <- goi.expression >= upper.quartile
  print(paste0("n with high expression of ", goi, " = ", sum(High.exp.logic)))

  # use threshold to define low expressing patients
  Low.exp.logic <- goi.expression <= lower.quartile
  print(paste0("n with low expression of ", goi, " = ", sum(Low.exp.logic)))

  # extract patient IDs
  patients.high <- names(goi.expression)[High.exp.logic]
  patients.low <- names(goi.expression)[Low.exp.logic]

  # create list for exporting
  output.list <- list(patients.high = patients.high,
                      patients.low = patients.low)

  return(output.list)

}else if(length(goi) > 1){


  print("need to add code for multiple goi thresholding")
}

}




#' @title Annotate clinical data with gene expression signature
#' @description Annotate clinical data with gene expression signature
#' @param patients.list list of patients with high and low gene expression signature
#' @param dataset.name name of dataset to retrieve clinical data from
#' @return dataframe of clinical data annotated with gene expression signature
#' @examples
#' annotate.clinical.data(patients.list, dataset.name)

annotate.clinical.data <- function(patients.list,
                                   dataset.name){

  print("Downloading clinical data")

  clin.data <- TCGAbiolinks::GDCquery_clinic(dataset.name, "clinical")

  print("Identifying samples by Gene expression signature")

  # identifying samples by signature
  high.patients <- stringr::str_extract(patients.list[["patients.high"]],
                                        pattern = "TCGA-..-....")

  # ensure only unique entries of patients
  high.patients <- unique(high.patients)

  patients.high.logic <- clin.data$bcr_patient_barcode %in% high.patients


  low.patients <- stringr::str_extract(patients.list[["patients.low"]],
                                       pattern = "TCGA-..-....")

  # ensure only unique entries of patients
  low.patients <- unique(low.patients)

  patients.low.logic <- clin.data$bcr_patient_barcode %in% low.patients

  # sanity checks
  # ensure number of patients clinical data recovered is similar to number of patients data requested

  print(paste0("length of patients annotated as high signature score = ", length(high.patients)))

  print(paste0("length of patients annotated as low signature score = ", length(low.patients)))

  print(paste0("retrieved clinical data for ", sum(patients.high.logic), " patients annotated as high goi"))

  print(paste0("retrieved clinical data for ", sum(patients.low.logic), " patients annotated as high goi"))

  # Add annotation to clinical data table
  print("Add annotation to clinical data table")

  clin.data$Signature[patients.high.logic] <- paste0("High")

  clin.data$Signature[patients.low.logic] <- paste0("Low")


  return(clin.data)

}


#' @title plot the overall survival

#' @description
#' Generate Kaplan-Meier Overall Survival plot
#'
#' @param annotated.clin.data output from annotate.clinical.data function
#' @param dataset.name name of TCGA dataset
#' @param output.dir file.path where plots should be saved
#' @param with.risk logical, should risk table be included
#' @param with.conf logical should confidence interval be included
#'
#' @return
#' saves an OS plot to the specified location

plot.OS <- function(annotated.clin.data,
                    dataset.name,
                    goi,
                    output.dir = "figures/",
                    with.risk = TRUE,
                    with.conf = TRUE){

  # create output directory
  if(!dir.exists(output.dir)){
    dir.create(output.dir, recursive = TRUE)
  }


  print("Generating OS plot")

  if(with.risk & with.conf){
    title.var <- ""
  }else if(!with.risk & with.conf){
    title.var <- "_noRisk"
  }else if(with.risk & !with.conf){
    title.var <- "_noCI"
  }else if(!with.risk & !with.conf){
    title.var <- "_noRisk_noCI"
  }


  # Full plot with CI and risk table
  TCGAbiolinks::TCGAanalyze_survival(annotated.clin.data,
                                     clusterCol = "Signature",
                                     filename = paste0(output.dir, dataset.name, "_overall_survival", title.var, ".pdf"),
                                     risk.table = with.risk,
                                     color = c("Red", "Blue"),
                                     main = paste0("Kaplan-Meier Overall Survival ", dataset.name, " ", goi, " High vs. Low"),
                                     conf.int = with.conf,
                                     height = 10,
                                     width = 10)


  print(paste0("OS Plot saved to file.path = ", output.dir))

}






#HNSC.raw.data <- load.rda("TCGA-HNSC_Transcriptome_Profiling_raw_counts.rda")
#normed.data <- preprocess.TCGA(input.data = HNSC.raw.data, save.normed.data = F)


#' @title Plot overall survival curve for a selected TCGA dataset based on gene signature
#' @description
#' This function will download TCGA RNAseq data, process it, threshold patients based on provided gene signature
#' then plot an overall surival curve
#' @param TCGA.dataset A character string indicating the TCGA dataset to be used e.g "TCGA-HNSC"
#' @param gene.signature A character string indicating the gene of interest or vector of strings indicating gene signature
#' @param gene.name A character string indicating the name of the gene of interest
#' @param rda.file.path A character string indicating the path to the rda file containing the raw data
#' @param plot.dir File path for where to save OS plot
#' @param data.dir File path for where to save data to
#'
#' @details
#' gene.name param allows the user to define a different name if the gene id is an uncharacteristic name.
#' For example, gene.signature = "CD226" but perhaps you want plots to be named "DNAM1" therefore set gene.name = "DNAM1"
#' If no gene.name is set, it automatically inherits gene.signature
#'
#' Provide a rda.file.path if data has already been downloaded. this will save signficant time
#'
#' @return
#' A plot of overall survival is saved to the desired output location
#' @export
#' @examples
#' TCGA.OS(TCGA.dataset = "TCGA-HNSC",
#'         gene.signature = "NKG7",
#'         rda.file.path = "TCGA-HNSC_Transcriptome_Profiling_raw_counts.rda")

TCGA.OS <- function(TCGA.dataset,
                    gene.signature,
                    gene.name = NULL,
                    rda.file.path = NULL,
                    plot.dir = "figures/",
                    data.dir = "data/"){

  # if gene.name is not defined then let gene.name = goi
  if(is.null(gene.name)){
    gene.name <- gene.signature
  }

  # Check if rda file path is provided, if yes then load data, if not then download data
  if(is.null(rda.file.path)){
    print("no file path provided, data will be downloaded")

    raw.data <- download.TCGA.RNAseq(dataset.name = TCGA.dataset,
                                     raw.data.dir = paste0(data.dir, "GDCdata/"),
                                     rda.dir = data.dir)

  }else{
    print("file path provided, data will loaded from file")

    raw.data <- load.rda(path.to.rda = rda.file.path)

  }

  # process dataset
  print("Processing and normalising dataset")

  normed.data <- preprocess.TCGA(input.data = raw.data,
                                 dataset.name = TCGA.dataset,
                                 quantile.filt = TRUE,
                                 save.normed.data = TRUE,
                                 output.dir = data.dir)


  # Extract expression, define threshold, and annotate patients as high or low exp
  print("Threshold and define high and low expression patients")

  threshold.list <- threshold.goi(input.data = normed.data,
                                  goi = gene.signature)


  # get clinical data and annotate patients as high or low exp
  print("Getting clinical data and defining patients as high or low")

  clinical.data <- annotate.clinical.data(patients.list = threshold.list,
                                          dataset.name = TCGA.dataset)

  # plotting overall survival
  print("Plotting OS")

  plot.OS(annotated.clin.data = clinical.data,
          dataset.name = TCGA.dataset,
          goi = gene.signature,
          output.dir = plot.dir,
          with.risk = TRUE,
          with.conf = TRUE)

}





























