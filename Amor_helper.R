convert_tibble_to_dataframe <- function(x){
  res <- as.data.frame(x[, -1])
  rownames(res) <- x[, 1] %>% unlist()
  return(res)
}

convert_data_to_tibble <- function(x, name="Gene"){
  res <- as.data.frame(x) %>% 
     rownames_to_column(name) %>% 
     as_tibble()
  return(res)
}




### integrate enrichR for multiple groups
# input needs to contains 2 columns
# Gene 
# Group
run_enrichR_groups <- function(gene.list, dbs="GO", prefix="enrichR", 
                               by.enrich="Combined.Score", top=8, height=11, width=9, 
                               summary=FALSE){
   require(enrichR)
   require(pheatmap)
   require(RColorBrewer)

  if(dbs == "GO"){
    print(paste("Use preset pathways (GO, KEGG, REACTOME, BIOCARTA)"))
    dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021" ,
             "WikiPathway_2021_Human", "Panther_2016", "KEGG_2021_Human", "KEGG_2019_Mouse", "BioCarta_2016",
             "Reactome_2022", "CellMarker_Augmented_2021", "MSigDB_Hallmark_2020",
             "Azimuth_Cell_Types_2021", "HuBMAP_ASCTplusB_augmented_2022") 
  }else if(dbs == "TF"){
    print(paste("Use preset TF pathways)"))
    summary <- TRUE
    dbs <- c("ChEA_2022", "TRANSFAC_and_JASPAR_PWMs", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
             "ENCODE_TF_ChIP-seq_2015" ,"Transcription_Factor_PPIs", "Epigenomics_Roadmap_HM_ChIP-seq",
             "TRRUST_Transcription_Factors_2019", "TF_Perturbations_Followed_by_Expression", "ARCHS4_TFs_Coexp")
  }else{
    print(paste("User defined dbs"))
  }
  
  type.col <- gene.list %>% 
    dplyr::select(-Gene) %>% 
    names
  enrich.res <- gene.list %>% 
    dplyr::rename(Type = !!sym(type.col)) %>% 
    group_by(Type) %>% 
    nest() %>% 
    mutate(res = purrr::map(data, function(x){ Sys.sleep(1); unlist(x) %>% enrichr(., dbs)}))
  
   # split by dbs
  by_dbs <- enrich.res %>% 
    dplyr::select(-data) %>% 
    unnest(res) %>% 
    mutate(cnt = map_dbl(res, nrow)) %>% 
    dplyr::filter(cnt > 0) %>%
    dplyr::select(-cnt) %>%
    mutate(dbs = names(res)) %>%
    split(., .$dbs) %>% 
    map(~ unnest(.x, res) %>% extract_top_enrichR_res(top=top, by=by.enrich))
  names(by_dbs) <- substr(names(by_dbs), 1, 30)
  writexl::write_xlsx(by_dbs, paste0(prefix, ".bydbs.xlsx"))
  
  pdf(paste0(prefix, ".top", top, ".pdf"), height=height, width=width)
  for(name in names(by_dbs)){
    cur.dat <- by_dbs[[name]] %>% 
       dplyr::mutate(Term = gsub("[_ ]Homo.*", "", Term)) %>%
       dplyr::mutate(Term = gsub(" WP.*", "", Term)) %>%
       dplyr::mutate(Term = gsub(" \\(GO.*", "", Term)) %>%
       dplyr::filter(!duplicated(Term)) %>%
       mutate(across(-contains("Term"), function(x) ifelse(is.infinite(x), 500, x)))

    if(nrow(cur.dat) >= 2 & ncol(cur.dat) >= 3){
       x2 <- convert_tibble_to_dataframe(cur.dat)
       fontsize_row = pheatmap_fontsize_row(x2)
       color = colorRampPalette(c("grey99", "#FFFFB3", "#FF7F00", "#E41A1C"))(60)

       pheatmap(x2, main=name, fontsize_row = fontsize_row, color=color, 
                border_color = "grey88", cellwidth=30,
                clustering_method = "ward.D2")
    }
  }
  dev.off()

   # split by Type
  by_type <- enrich.res %>% 
    dplyr::select(-data) %>% 
    unnest(res) %>% 
    mutate(cnt = map_dbl(res, nrow)) %>% 
    dplyr::filter(cnt > 0) %>%
    dplyr::select(-cnt) %>%
    mutate(dbs = names(res)) %>%
    unnest(res) %>% 
    split(., .$Type) 
  names(by_type) <- substr(names(by_type), 1, 35)
  writexl::write_xlsx(by_type, paste0(prefix, ".byType.xlsx"))

  ### three columns for go_barplot, Term, Score, Pathway
  pdf(paste0(prefix, ".bar.pdf"), height=13, width=12)
  aa <- by_type %>% 
    map(~ .x %>% 
          dplyr::mutate(Term = gsub(" \\(GO.*", "", Term)) %>%
          dplyr::mutate(Term = gsub("[_ ]Homo.*", "", Term)) %>%
          dplyr::mutate(Term = gsub(" WP.*", "", Term)) %>%
          dplyr::rename(Score = Combined.Score) %>% 
          mutate(Pathway = dbs) %>% 
          dplyr::select(Term, Score, Pathway) %>%
          group_by(Term) %>% 
          arrange(desc(Score)) %>% 
          dplyr::slice(1) %>% 
          ungroup() %>% 
          group_by(Pathway) %>%
          arrange(Pathway, desc(Score)) %>%
          dplyr::slice(1:top) %>%
          arrange(Pathway, Score) %>%
          ungroup() %>%
          dplyr::mutate(Term = factor(Term, levels=Term %>% unique))
    ) %>% 
    map2(paste(names(.), by.enrich, sep="."), ., go_barplot)
  dev.off()
  
  # return the data for further processing
  if(summary){
    summarise_TF <- enrich.res %>% 
        dplyr::select(-data) %>% 
        unnest(res) %>% 
        bind_cols(tibble(dbs = rep(dbs, length.out=nrow(.)))) %>%
        unnest(res) %>% 
        arrange(desc(Combined.Score)) %>%
        mutate(Term = gsub("_.*", "", Term) %>% gsub(" .*", "", .)) %>% # Grep TF name
        distinct(Type, dbs, Term) %>%  # remove Term (TF) that have multiple entries
        group_by(Type, dbs) %>% 
        dplyr::slice(1:top) %>%  # select top Term (TF) from each Type
        ungroup() %>%
        group_by(Type, Term) %>%  # count how many times each Term (TF) appears 
        summarise(Freq = n()) %>% 
        ungroup() %>% 
        arrange(desc(Freq)) %>% 
        mutate(Term = factor(Term, levels=unique(Term))) %>%  # keep the order of the Term (TF)
        spread(Type, Freq) %>% 
        replace(is.na(.), 0) %>% 
        arrange_at(-1, dplyr::desc)
    write_csv(summarise_TF, paste0(prefix, ".summary.csv"))
  }

  return(enrich.res)
}




# v2. just to be save, loop one by one and to match, replace
replace_ModuleScore_feature_name_v2 <- function(object, feature_list, name = "Cluster"){
  for(i in length(feature_list):1){
    colnames(object@meta.data) <- gsub(paste0(name, i), names(feature_list)[i], colnames(object@meta.data))
  }
  return(object)
}
