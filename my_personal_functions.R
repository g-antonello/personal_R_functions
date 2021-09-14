#################################
# decluttering function based on having "tmp" in the name, not too beautiful, but works
declutter.environment.tmp_test <- function(){
  # remove those containing 
  # "tmp"     | OK
  # "test"    | OK
  #  nchar==1 | OK
  objs.t0 <- ls(envir = parent.frame(1))
  obj_to_del <- grep("tmp|test", ls(envir = parent.frame(1)), value = TRUE)
  obj_to_del <- obj_to_del[2:length(obj_to_del)]
  obj_to_del_final <- c(obj_to_del, ls(envir = parent.frame(1))[nchar(ls(envir = parent.frame(1))) ==1])
  rm(list = obj_to_del_final, envir = parent.frame(1))
  objs.t1 <- ls(envir = parent.frame(1))
  cat("elements removed:\n")
  return(objs.t0[!(objs.t0 %in% objs.t1)])
}

##################################
# palette function

getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))  # change the palette as well as the number of colors will change according to palette.

# function to substitute metadata in phyloseq objects. the two data frames must have row names matching sample names, and the sample order doesn't matter, because the functions will reorder the new metatadata rows, to match the sample order in the phyloseq object.
substitute_metadata <- function(physeq, new_metadata){
  if(any(!rownames(new_metadata)%in%sample_names(physeq))){ # this is like asking: is any item of the vector NOT true? if responds true, it means rownames are not all there.
    stop ("not all row names of the metadata are in the sample names of the phyloseq object")
  }
  # reorganize metadata
  new_metadata <- new_metadata[match(sample_names(physeq), rownames(new_metadata)),]
  # remove old metadata
  physeq@sam_data <- NULL
  sample_data(physeq) <- sample_data(new_metadata)
  return(physeq)
}

############################################
# ref Sudarshan Shetty, Gerben Hermes, Wageningen University and Research
getPalette = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired")) 

show_palette_colors <- function(html_codes){
  image(1:length(html_codes), 1, as.matrix(1:length(html_codes)), 
        col = html_codes, 
        xlab = "", ylab = "")}

# option 2
StatN <- ggproto("StatN", Stat,
                 required_aes = c("x", "y"), 
                 compute_group = function(data, scales) {
                   y <- data$y
                   y <- y[!is.na(y)]
                   n <- length(y)
                   data.frame(x = data$x[1], y = min(y), label = paste0("n=", n))
                 }
          )
stat_n <- function(mapping = NULL, data = NULL, geom = "text", 
                   position = "identity", inherit.aes = TRUE, show.legend = NA, 
                   na.rm = FALSE, ...) {
  ggplot2::layer(stat = StatN, mapping = mapping, data = data, geom = geom, 
                 position = position, inherit.aes = inherit.aes, show.legend = show.legend, 
                 params = list(na.rm = na.rm, ...))
}
######################################################################
# merge a dataframe containing an ASV column and adds some levels of phylogeny based on a phyloseq object

add_ASV_name_to_ASV_results <- function(df,col_ASV,  physeq, tax_lvls = "all") {
  tax <- data.frame(physeq@tax_table@.Data)
  if(tax_lvls != "all"){
    tax <- tax[, c(col_ASV,tax_lvls)]
  }
  m <- merge(df, tax, by = col_ASV)  
  return(m)
}

#####################################################################
# function to detect outliers in data prior to ggplot
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.2 * IQR(x) | x > quantile(x, 0.75) + 1.2 * IQR(x))
}

#####################################################################
# function taking a table as input and adding percentage values 
add_percent_columns <- function(tb){
  new_tb <- matrix(0,nrow = nrow(tb), ncol = ncol(tb)*2)
  rownames(new_tb) <- rownames(tb)
  colnames(new_tb) <- rep(c("abs", "%"), ncol(tb))
  for(i in 1:ncol(tb)){

    new_tb[,c(2*i-1, 2*i)] <- cbind(tb[,i], round(tb[,i]/sum(tb[,i])*100, digits = 1))

  }
  
  return(new_tb)
}
#####################################################################
# function to generate demographic table ------------------ NOT WORKING YET!!
make_demogr_table <- function(df, row_variables_as_are, row_var_names_nicer, col_variable, col_var_names_nicer, add_perc_column, caption){
  print("Doesn't yet support NAs in the variables of interest")
  # make list of tables, transformed as matrices
  l <- lapply(row_variables_as_are, function(x) as.matrix(table(df[[x]], df[[col_variable]]))) # for each row variable, build a table, in a matrix form
  
  # in case you want a percent column next to the original column
  if (add_perc_column){
    l <- lapply(l, add_percent_columns)
    headr <- c("1",rep("2", length(col_var_names_nicer))) # initiate a header character vector, with length = 2*columns of the input table
    final_table_matr <- do.call(rbind, l) # bind all tables in the "l" list
  }else{
    # even if you don't want the relative population, still you need a header, but with the same length of the initial matrix
    headr <- c("1",rep("1", length(col_var_names_nicer)))
    final_table_matr <- do.call(rbind, l) # bind all tables in the "l" list
    colnames(final_table_matr) <- NULL
  }
  
  
  # initiate the header, the first column will be empty
  names(headr) <- c(" ", col_var_names_nicer)
  
  library(kableExtra)
  kb <- kableExtra::kbl(final_table_matr[rowSums(final_table_matr) != 0,],format = "html", caption = caption) %>% 
    kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
    add_header_above(headr)

  row_point <- 1
  row_names_both <- cbind(row_variables_as_are, row_var_names_nicer)
  for(r in 1:nrow(row_names_both)){
    kb <- pack_rows(kable_input = kb,
                    group_label = row_names_both[r,2],
                    start_row = row_point,
                    end_row = length(na.omit(unique(as.character(df[[row_names_both[r,1]]])))),
                    label_row_css = "background-color: #666; color: #fff;",indent = F
                    )

    row_point <-  row_point + length(na.omit(unique(as.character(df[[row_names_both[r,1]]]))))
    #print(length(na.omit(unique(as.character(df[[r]])))))

  }
  return(kb)

}


######################################################################
# function to bin a numeric variable into its upper value
bin_integer_into_integer <- function(x, multiples, include_highest = TRUE){
  brks <- seq(0, max(x)+multiples, multiples) # divided by 7 and a very small value, to include 
  x_fact <- as.character(cut(x, breaks = brks)) %>%
    strsplit(split = "\\(|\\]") %>% sapply("[", 2) %>% 
    strsplit(split = ",", fixed = TRUE) %>% sapply("[", 2) %>% 
    as.numeric()
  x_fact_divided <- x_fact/multiples
  return(x_fact_divided)
}
  
######################################################################
# convert data.frame columns all at once
convert_col_classes <- function(df,types){
  out <- lapply(1:length(df),FUN = function(i){FUN1 <- switch(types[i],character = as.character,numeric = as.numeric,factor = as.factor); FUN1(df[,i])})
  names(out) <- colnames(df)
  as.data.frame(out,stringsAsFactors = FALSE)
}

######################################################################
# boxplot of beta diversities within- and between- two levels of a factor

beta_div_betw_withn <- function(physeq, dist, variab, keep_sample_lists = TRUE){

  # preparatory part
  if(class(dist) %in% c("matrix", "dist")){
    print("using dist as distance matrix/distance")
    dist_mat <- as.matrix(dist)
  }
  
  if(class(dist) == "character"){
    print("need to compute the beta diversity")
    dist <- distance(physeq, method = dist)
    dist <- as.matrix(dist)
  }
  
  #get only upper triangle values, remove also the diagonal of 0
  dist_upper <- dist_mat
  dist_upper[lower.tri(dist_upper, diag=TRUE)] <- NA
  
  # melt, remove NAs, and rename the columns
  dist_upper_molten <- reshape2::melt(dist_upper) %>% na.omit()
  colnames(dist_upper_molten) <- c("samp1", "samp2", "value")
  
  #list of sample names for each category in "variab" 
  split_samp_names <- meta(physeq) %>% 
    rownames_to_column("samp_names") %>% 
    select("samp_names", variab) %>% 
    dplyr::group_split(!!rlang::sym(variab), .keep = FALSE) %>% 
    lapply(function(d) d$samp_names)
  names(split_samp_names) <- unique(as.character(meta(physeq)[[variab]]))

  # actual value extraction part
  # get within-group diversity indices
  within_divs <- lapply(split_samp_names, function(s)
    dist_upper_molten %>% 
      filter(samp1 %in% s & samp2 %in% s) %>% 
      .$value
  )
  
  # get between-group diversity indices
  between_list_samples_and_values <- combn(names(split_samp_names), 2, simplify = FALSE) %>% 
    lapply(function(n) split_samp_names[n]) 
  
  names(between_list_samples_and_values) <- sapply(between_list_samples_and_values, function(n) paste(names(n), collapse = " VS "))
  
  between_list_samples_and_values <- between_list_samples_and_values %>% 
    lapply(function(n) dist_upper_molten %>% filter((samp1 %in% n[[1]] & samp2 %in% n[[2]])|(samp2 %in% n[[1]] & samp1 %in% n[[2]]))) 
  
  # save the results as a list
  if(keep_sample_lists){
    betadiv_results <- list(
      samples_lists = split_samp_names,
      within = within_divs,
      between = between_list_samples_and_values
    )}else{
      betadiv_results <- list(
        within = within_divs,
        between = lapply(between_list_samples_and_values, function(df) df$value)
      )
    }
  
  
  return(betadiv_results)
}

######################################################################
# check if a vector is made of nucleotide sequences (also works with RNA and nucleotides = "n" or "N")
is_nucleotide <- function(x){
  l <- sapply(strsplit(x, "*"), function(i) as.character(i) %in% 
                c("a", "A", "c", "C", "g", "G", "t", "T", "u", "U", "n", "N"))
  s <- sapply(l, function(i) all(all(i == TRUE) ==TRUE))
  return(all(s==TRUE))
  
}
###########################################################################################
# generate a fasta file out of reference sequences in your dataset
phyloseq_to_fasta <- function(physeq, seq_location = "tax_table", destination_file = NULL){
  if(seq_location == "tax_table"){
    tx <- as.data.frame(physeq@tax_table@.Data)
    pos <- sapply(tx, is_nucleotide)
    seqs <- as.character(tx[, pos])
    names(seqs) <- rownames(tx)
    
    if (is.null(destination_file)){
      return(seqs)
      }else{
    
    cat(file = destination_file, paste(">", names(seqs), "\n", seqs, "\n", sep = ""), sep = "")
    cat(paste("\nfasta file written: '", tools::file_path_as_absolute(destination_file), "'", sep = ""))
      } 
  }
  if(seq_location == "refseq"){
    seqs <- refseq(physeq)
    return(seqs)
    cat("\nrefseq extraction is not implemented yet")
  }
}

######################################################################
# sample conditional means and standard deviation
samp_cond_means_sd <- function(x, g, na.rm = TRUE){
  x <- x[!is.na(g)]
  g <- g[!is.na(g)]
  mean <- sapply(unique(g), function(i) mean(x[grepl(i,g)], na.rm = TRUE))
  std <- sapply(unique(g), function(i) sd(x[grepl(i,g)], na.rm = TRUE))
  names(mean) <- unique(g)
  return(rbind(mean, std))
}

################################################################
# sample conditional medians and interquartile range
samp_cond_medians_IQR <- function(x, g, na.rm = TRUE){
  x <- x[!is.na(g)]
  g <- g[!is.na(g)]
  median <- sapply(unique(g), function(i) median(x[grepl(i,g)], na.rm = TRUE))
  iqr <- sapply(unique(g), function(i) IQR(x[grepl(i,g)], na.rm = TRUE))
  names(median) <- unique(g)
  return(rbind(median, iqr))
}

##################################################################
# geometric mean of a data matrix
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

###############################################################################
# clone a directory and subdirectory's structure into another location
# can specify both an existing and a non-existing target directory location (end_dir)

clone_dir_structure <- function(start_dir, end_dir){
  
  if (!endsWith(start_dir, "/")){ # need to remove "/" if it is at the end of the input text
    start_dir <- paste(start_dir, "/", sep = "")
  } # strip that character from the end of the start_dir
  
  d <- list.dirs(start_dir)
  d_clean <- stringr::str_remove(d, start_dir)
  d_clean <- d_clean[nchar(d_clean) != 0]
  
  while(any(startsWith(d_clean, "/"))){
    d_clean[startsWith(d_clean, "/")] <- substr(x = d_clean[startsWith(d_clean, "/")], 
                                                start = 2, 
                                                stop = nchar(d_clean[startsWith(d_clean, "/")]))
  }
  
  end_directories <- paste(paste(end_dir, d_clean, sep = ""), "/", sep = "")  
  sapply(end_directories, function(i) dir.create(i, recursive = TRUE))
}


transform_microbiota_ga <- function(physeq, transform){
  # transform otu counts
  if(transform %in% c("asinh", "arcsinh")){
    if(taxa_are_rows(physeq)){
      transformed_otu <- apply(as.matrix(otu_table(physeq)), 1, asinh)
    }else{
      transformed_otu <- apply(as.matrix(otu_table(physeq)), 2, asinh)
    }
    physeq_transf <- physeq
    otu_table(physeq_transf) <- otu_table(transformed_otu,taxa_are_rows = taxa_are_rows(physeq))
    # substitute otu table with newly transformed one
  }else{
    physeq_transf <- microbiome::transform(physeq, transform = transform)
  }
  
  return(physeq_transf)
}
####################################################################
boxplot_taxon <- function(physeq, taxon, trait, fill = NULL, colour = NULL, transform= "clr", jitter_alpha = NULL){
  
  physeq_transf <- transform_microbiota_ga(physeq = physeq, transform = transform)
  
  # prepare data
  if (taxa_are_rows(physeq_transf)){
    otutab <- t(physeq_trasnf@otu_table@.Data)
    
  }else{ #taxa are not rows
    otutab <- physeq_transf@otu_table@.Data
  }
  
  colnames(otutab) <- make.names(colnames(otutab))
  
  identical(rownames(as.data.frame(physeq_transf@sam_data)), rownames(otutab))
  data_merged <- as.data.frame(cbind(otutab, as.data.frame(physeq_transf@sam_data)))
  data_merged[[trait]] <- as.factor(data_merged[[trait]])
  if (!is.null(jitter_alpha)){
    boxp <- ggplot(data = data_merged, aes_string(x = trait, y = make.names(taxon))) + geom_boxplot(aes_string(colour=colour, fill = fill)) + geom_jitter(alpha = jitter_alpha)  + ggtitle(paste("trasnformation method: ", transform, sep = ""))
  }else{
    boxp <- ggplot(data = data_merged, aes_string(x = trait, y = make.names(taxon))) + geom_boxplot(aes_string(colour=colour, fill = fill)) + ggtitle(paste("transformation method: ", transform, sep = ""))
  }
  
  return(boxp)
}


##########################################################################################
# same as above, but scatterplot

scatterplot_taxon <- function(physeq, taxon, trait, fill = NULL, colour = NULL, transform= "clr", jitter_alpha = NULL){
  # transform otu counts
  physeq_transf <- microbiome::transform(physeq, transform = transform) 
  # prepare data
  if (taxa_are_rows(physeq_transf)){
    otutab <- t(physeq_trasnf@otu_table@.Data)
    
  }else{ #taxa are not rows
    otutab <- physeq_transf@otu_table@.Data
  }
  identical(rownames(as.data.frame(physeq_transf@sam_data)), rownames(otutab))
  data_merged <- as.data.frame(cbind(otutab,as.data.frame(physeq_transf@sam_data)))
  if (!is.null(jitter_alpha)){
    scatpl <- ggplot(data = data_merged, aes_string(x = trait, y = taxon)) + geom_point(aes_string(colour=colour, fill = fill)) + geom_jitter(alpha = jitter_alpha)  + ggtitle(paste("trasnformation method: ", transform, sep = ""))
  }else{
    scatpl <- ggplot(data = data_merged, aes_string(x = trait, y = taxon)) + geom_point(aes_string(colour=colour, fill = fill)) + ggtitle(paste("transformation method: ", transform, sep = ""))
  }
  
  return(scatpl)
}

##########################################################################################
############ Function to remove na's from wanted variables in a phyloseq object

drop_NA_from_phyloseq_vars <- function(physeq, vars){
  physeq_cleaned <- physeq
  for(v in vars){
    print(v)
    samples_to_save <- meta(physeq_cleaned) %>% 
      filter(!is.na(!!rlang::sym(v))) %>% 
      rownames()
    tmp <- sample_names(physeq_cleaned) %in% samples_to_save
    physeq_cleaned <- prune_samples(samples = samples_to_save, x =  physeq_cleaned)
  }

  
  return(physeq_cleaned)
}


##########################################################################################
# plot_bar function, with customizable border
  
plot_bar2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA){
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

##################################################################################################
# to write all combinations of results for a DESeq2 object and a trait of interest
write_deseq_results <- function(deseq_results, trait, target_dir, physeq_obj = NULL, physeq_aggr_lvl, taxon_lvls_added = "all", sort_p.val = TRUE, shorten_name = FALSE){
  combin <- (t(combn(levels(deseq_results[[trait]]),2)))
  library(DESeq2)
  for (l in 1: nrow(combin)){
    if(sort_p.val){
      if(!is.null(physeq_obj)){ # if you want to arrange p value and add taxonomic info
        write.csv(results(deseq_results, contrast =c(trait, combin[l,1], combin[l,2])) %>% as.data.frame() %>% 
                    tibble::rownames_to_column(physeq_aggr_lvl) %>%
                    add_ASV_name_to_ASV_results(., col_ASV = physeq_aggr_lvl, physeq = physeq_obj, tax_lvls = taxon_lvls_added)%>% arrange(padj),
                  file = paste(target_dir, ifelse(shorten_name,"", "data=chrismb; "), "trait=", trait,"; ","contrast=", combin[l,1], "_vs_", combin[l,2], ".csv", sep = ""))
      } else{ # arrange p value but no taxonomic info
        write.csv(results(deseq_results, contrast =c(trait, combin[l,1], combin[l,2])) %>% as.data.frame() %>% 
                    tibble::rownames_to_column(physeq_aggr_lvl) %>%
                    arrange(padj),
                  file = paste(target_dir, ifelse(shorten_name,"", "data=chrismb; "), "trait=", trait,"; ","contrast=", combin[l,1], "_vs_", combin[l,2], ".csv", sep = ""))
      }
    } else{ # do not arrange p value
      if(!is.null(physeq_obj)){ # do not arrange p value but add ASV info
        write.csv(results(deseq_results, contrast =c(trait, combin[l,1], combin[l,2])) %>% as.data.frame() %>%
                    tibble::rownames_to_column(physeq_aggr_lvl)%>%
                    add_ASV_name_to_ASV_results(., col_ASV = physeq_aggr_lvl, physeq = physeq_obj, tax_lvls = taxon_lvls_added),
                  file = paste(target_dir, ifelse(shorten_name,"", "data=chrismb; "), "trait=", trait,"; ","contrast=", combin[l,1], "_vs_", combin[l,2], ".csv", sep = ""))
      } else{ # do not arrange p value nor add ASV info
        write.csv(results(deseq_results, contrast =c(trait, combin[l,1], combin[l,2])) %>% as.data.frame()%>%
                    tibble::rownames_to_column(physeq_aggr_lvl),file = paste(target_dir, ifelse(shorten_name,"", "data=chrismb; "), "trait=", trait,"; ","contrast=", combin[l,1], "_vs_", combin[l,2], ".csv", sep = ""))
      }
      
    }
    
  }
  
}

############
# this function does the same as the write_deseq... but returns you directly the data frame list 
deseq_results_into_list <- function(deseq_results, trait, physeq_obj = NULL, physeq_aggr_lvl, taxon_lvls_added = "all", sort_p.val = TRUE, shorten_name = FALSE){
  combin <- (t(combn(levels(deseq_results[[trait]]),2)))
  library(DESeq2)
  
  results_list <- list()
  for (l in 1: nrow(combin)){
    if(sort_p.val){
      if(!is.null(physeq_obj)){ # if you want to arrange p value and add taxonomic info
        results_list[[paste(combin[l,1], "_vs_", combin[l,2], sep = "")]] <- results(deseq_results, contrast =c(trait, combin[l,1], combin[l,2])) %>% as.data.frame() %>% 
                    tibble::rownames_to_column(physeq_aggr_lvl) %>%
                    add_ASV_name_to_ASV_results(., col_ASV = physeq_aggr_lvl, physeq = physeq_obj, tax_lvls = taxon_lvls_added)%>% arrange(padj)
      } else{ # arrange p value but no taxonomic info
        results_list[[paste(combin[l,1], "_vs_", combin[l,2], sep = "")]] <- results(deseq_results, contrast =c(trait, combin[l,1], combin[l,2])) %>% as.data.frame() %>% 
                    tibble::rownames_to_column(physeq_aggr_lvl)
      }
    } else{ # do not arrange p value
      if(!is.null(physeq_obj)){ # do not arrange p value but add ASV info
        results_list[[paste(combin[l,1], "_vs_", combin[l,2], sep = "")]] <- results(deseq_results, contrast =c(trait, combin[l,1], combin[l,2])) %>% as.data.frame() %>%
                    tibble::rownames_to_column(physeq_aggr_lvl)%>%
                    add_ASV_name_to_ASV_results(., col_ASV = physeq_aggr_lvl, physeq = physeq_obj, tax_lvls = taxon_lvls_added)
      } else{ # do not arrange p value nor add ASV info
        results_list[[paste(combin[l,1], "_vs_", combin[l,2], sep = "")]] <-  results(deseq_results, contrast =c(trait, combin[l,1], combin[l,2])) %>% as.data.frame()%>%
                    tibble::rownames_to_column(physeq_aggr_lvl)
      }
      
    }
    
  }
  return(results_list) 
}

###################################################################################
phyloseq_to_edgeR <- function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}


#############################################################################

save_tags <- function (tags, file, selfcontained = F, libdir = "./lib") 
{
  if (is.null(libdir)) {
    libdir <- paste(tools::file_path_sans_ext(basename(file)), 
                    "_files", sep = "")
  }
  htmltools::save_html(tags, file = file, libdir = libdir)
  if (selfcontained) {
    if (!htmlwidgets:::pandoc_available()) {
      stop("Saving a widget with selfcontained = TRUE requires pandoc. For details see:\n", 
           "https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md")
    }
    htmlwidgets:::pandoc_self_contained_html(file, file)
    unlink(libdir, recursive = TRUE)
  }
  return(file)
}
###########################################################################
# from GitHub/jasdumas, https://github.com/jasdumas/data-exploreR/blob/master/www/flattenCorrMatrix.R

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
############################################################################

find_char_in_string <- function(str, chars_to_check, is.there = TRUE){
  if(is.there){
    v <- chars_to_check[sapply(chars_to_check, function(i) grepl(i, str))]
    return(v)  
  }else{
    v <- chars_to_check[!sapply(chars_to_check, function(i) grepl(i, str))]
    
  }
  
}

##########################################################################
extract_corncob_results_catecorical <- function(res, variable_of_interest, contrast_lvl ="choose another level of the factor, remember that the first is the baseline", tax_lvl, differential_test = "choose between 'Abundance', 'Variability'"){
  p <- plot(res, tax_lvl)$data
  pvals_df <- data.frame(taxa = names(res$p), 
                    pval = res$p, 
                    padj = res$p_fdr)
  
  var_contrasts <- strsplit(p$variable, "\n") %>% sapply("[", 1) %>% strsplit(variable_of_interest, fixed = T) %>% sapply("[", 2)
  differential_operation <- strsplit(p$variable, "\n") %>% sapply("[", 2) %>% strsplit(" ") %>% sapply("[", 2) 
  p$variable <- var_contrasts
  p$type_of_differential_operation <- differential_operation
  
  p <- merge(p, pvals_df, all.x = TRUE, by = "taxa")
  colnames(p) <- c(tax_lvl,"log2FoldChange", "stderror-", "stderror+", "contrast", "type_of_differential_operation", "pval", "padj")
  p <- p[,c(1:4,7,8,5,6)]
  
  # filter the type of differential analysis wanted
  p_subs <- p %>% filter(contrast == contrast)
  if(differential_test == "both"){
    p %>% filter(contrast == contrast) %>%
      return()
  } else{
    p %>% filter(type_of_differential_operation == differential_test) %>% filter(contrast == contrast_lvl) %>%
      return()
  }
  
}


#################################################
### custom correlation matrix and p value
cor.pmat <- function(x,method = corr.method, ...) {
  mat <- as.matrix(x)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  colnames(p.mat) <- colnames(mat)
  rownames(p.mat) <- colnames(mat)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      #print(tmp$p.value)
      p.mat[i, j] <- tmp$p.value
      p.mat[j, i] <- tmp$p.value
    }
  }
  
  return(as.matrix(p.mat))
}

corr_matrix_and_signif_custom <- function(mat_or_df, corr.method){

  l <- list(matrix_r = cor(mat_or_df, method = corr.method),
            matrix_p = cor.pmat(mat_or_df, method = corr.method)
            )
  return(l)
}

#################################################################
## extract corncob results for the significant taxa
extract_corncob_results_continuous <- function(corncob_res, which_differential = "DA", variab_estim_to_extract, taxa_to_extract = "significant"){
  if (which_differential == "both"){
    extrc <- paste(c("mu", "phi"), variab_estim_to_extract, sep = ".") 
  }
  if(which_differential == "DA"){
    extrc <- paste("mu", variab_estim_to_extract, sep = ".")  
  }
  if(which_differential == "DV"){
    extrc <- paste("phi", variab_estim_to_extract, sep = ".")  
  }
  
  if(taxa_to_extract == "significant"){
  tmp <- corncob_res$significant_models
  names(tmp) <- corncob_res$significant_taxa
  
  padj_df <- corncob_res$p_fdr[names(corncob_res$p_fdr) %in% corncob_res$significant_taxa] %>% 
    as.data.frame() %>% 
    rownames_to_column("Genus") %>% 
    dplyr::rename(., padj = .)
 
  tmp_signif_coeffs_corncob <- lapply(tmp , coefficients) %>% 
    lapply(function(i) as.data.frame(i) %>% rownames_to_column("variab")) %>% 
    lapply(function(i) i[i[,1] %in% extrc,]) %>%
    bind_rows(.id = "Genus")
  
  
  colnames(tmp_signif_coeffs_corncob) <- c("Genus", "estim_name","estimate", "estimate.SE", "t_value", "pvalue")
  
  extracted_output <- merge(tmp_signif_coeffs_corncob, padj_df, by = "Genus", all.x = TRUE, sort = FALSE) %>% 
    mutate(method = "corncob")
  } 
  if (taxa_to_extract == "all"){
    tmp <- corncob_res$all_models
    names(tmp) <- names(corncob_res$p)
    
    tmp_signif_coeffs_corncob <- lapply(tmp , coefficients) %>% 
      lapply(function(i) as.data.frame(i) %>% rownames_to_column("variab")) %>% 
      lapply(function(i) i[i[,1] %in% extrc,]) %>%
      bind_rows(.id = "Genus")
    
    
    colnames(tmp_signif_coeffs_corncob) <- c("Genus", "estim_name","estimate", "estimate.SE", "t_value", "pvalue")
    
    extracted_output <- tmp_signif_coeffs_corncob %>% 
      mutate(method = "corncob")
  }
  
  return(extracted_output)
}

###################################################

open.plot.directory <- function(plot_path){
  
  if (.Platform$OS.type == "windows"){# in case system is windows we need to change "/" into "\"
    if(!(str_sub(plot_path, start = nchar(plot_path), end = nchar(plot_path)) == "/")){ # if the directory name given ends with "/", then we can play with it without manipulation, otherwise, we need to do some managing
      s <- strsplit(plot_path, "/")[[1]]
      s_sub <- s[1:length(s)-1]
      plot_path_windows <- paste(s_sub, collapse = "\\") %>% paste("\\", sep = "")
    } else{ # if it ends with "/", just translate the "/" character
        plot_path_windows <- gsub(pattern = "/", replacement = "\\\\", x = plot_path)
    }
    
    shell.exec(plot_path_windows)
  } else{ # if it is not windows
      if(!endsWith(plot_path, "/")){# if the path does not end with "/", then cut to the first piece ending with it, 
                                    # because we need to open directories
        s <- strsplit(plot_path, "/")[[1]] 
        s_sub <- s[1:length(s)-1] 
        plot_path <- paste(s_sub, collapse = "/") %>% paste0("/")
      }
    # correct the issue about the spaces, because bash needs spaces to be escaped.
    # don't get scared by the double \\, because it is R's way of displaying a single one
    plot_path_unix <- gsub(pattern = " ", replacement = "\\ ", x = plot_path,fixed = TRUE)
    system(paste(Sys.getenv("R_BROWSER"), plot_path_unix))
  }
}

###############################################################################NOT READY
# function to merge metadata of phyloseq to ordintated distance
metadata_and_ordination <- function(ord, phy, n_axes){
  variance_explained <- round(ord$values[,2]*100, 1)
  vectors <- ord$vectors[,1:n_axes]
  colnames(vectors) <- paste(colnames(vectors), " (",variance_explained[n_axes], "%)", sep = "") 
  return(vectors)
}

###################################################################
save_tmp_ggfigure <- function(Robj, format, ...){
  tmp_hour <- Sys.time() %>% as.character() %>% strsplit(" ") %>% sapply("[", 2) 
  tmp_hour_nicer <- paste("time=", gsub(x = tmp_hour, pattern = ":", replacement = ".", fixed = T), sep = "")
  final_name <- paste("date=", Sys.Date()," ; ", tmp_hour_nicer, sep = "")
  # tmp_name <- 
  suppressWarnings(dir.create("~/../Desktop/tmp/R_tmp_figures/", recursive = TRUE))
  ggsave(plot = Robj, filename = paste("~/../Desktop/tmp/R_tmp_figures/",final_name, ".", format, sep = ""), dpi = ..., height = ..., width = ...)
}
###################################################################

phy_to_ped <- function(physeq, individual_id, taxa, trasform, covariates, dest_directory){
  
  # transform sample counts, 
  # allowed ones are inside microbiome::transform and the arcsinh
  physeq <- transform_microbiota_ga(physeq = physeq, transform = transform)  
  
  # get otu table
  if(taxa_are_rows(physeq)){
    otutable <- t(as.matrix(otu_table(physeq)))
  }else{
    otutable <- as.matrix(otu_table(physeq))
  }
  otutable <- data.frame(otutable)
  
  # get metadata table
  metatable <- microbiome::meta(physeq)
  
  # get only variables of interest
  metatable_onlyvarsofInterest <- metatable[,c(rep(individual_id, 4), 
                                          covariates)] %>% 
    cbind(otutable[,taxa]) %>% 
    data.frame() %>% 
    mutate(sex = ifelse(sex == "male", 1, 2))
  
  colnames(metatable_onlyvarsofInterest) <-  c("#FAM_ID", "IND_ID", "FAT_ID", "MOT_ID", covariates, taxa)
  rownames(metatable_onlyvarsofInterest) <- NULL
  
  if(!is.null(dest_directory)){
    write.table(x = metatable_onlyvarsofInterest, file = paste(dest_directory, "/", "mbGWAS_PED; taxa=", paste(taxa, collapse ="&"), ".ped", sep = ""), quote = FALSE)
  }
  return(metatable_onlyvarsofInterest)
}

###################################################################
