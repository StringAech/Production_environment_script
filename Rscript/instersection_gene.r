List2ggupset <- function(list) {
  gene_pathway_membership <- fromList(list) %>%
    apply(2, as.logical)
  rownames(gene_pathway_membership) <- unique(unlist(list))
  
  tidy_pathway_member <- gene_pathway_membership %>%
    t() %>%
    as_tibble(rownames = "Pathway") %>%
    gather(Gene, Member, -Pathway) %>%
    filter(Member) %>%
    select(-Member) %>%
    group_by(Gene) %>%
    summarise(Pathways = list(Pathway))
  return(tidy_pathway_member)
}
output_intersection <- function(intersecion_set, file){
  intersecion_set$Pathways <- as.character(intersecion_set$Pathways)
  intersecion_set <- split(intersecion_set, intersecion_set$Pathways) %>% 
    lapply(\(.x) dplyr::select(.x, Gene))
  names(intersecion_set) <-
    names(intersecion_set) %>% 
    gsub("c\\(", "", .) %>% 
    gsub("\\)", "", .)
  intersecion_set <- lapply(names(intersecion_set), \(.x){
    intersecion_set[[.x]] <- intersecion_set[[.x]] %>% 
      dplyr::mutate(Group = .x)
  }) %>% 
    do.call(rbind, .)
  return(intersecion_set)
}


intersecion_set <- List2ggupset(DEG_list) %>% 
  output_intersection()

write.xlsx(x = intersecion_set, file = "~/work/asd.xlsx")