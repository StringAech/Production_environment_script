DSP_subset <- function(object, patter = NULL){
  object@anno <- object@anno %>% 
    dplyr::filter(!!rlang::parse_expr(patter))
  object@raw_exp <- object@raw_exp %>% 
    as.data.frame() %>% 
    dplyr::select(all_of(rownames(object@anno))) %>% as.matrix()
  object@exp <- object@exp %>% 
    as.data.frame() %>% 
    dplyr::select(all_of(rownames(object@anno))) %>% as.matrix()
  return(object)
}
