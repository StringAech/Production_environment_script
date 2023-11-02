tidyfunction <- function(list,ncol=ncol){
  for (i in 1:(length(list)-ncol)) {
    if(i %% ncol == 1){
      list[[i]] <- list[[i]]+
        theme(#axis.text.x = element_blank(), 
          #axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
        )
    }else{
      list[[i]] <- list[[i]]+
        theme(#axis.text.y = element_blank(), 
          #axis.ticks.y = element_blank(), 
          axis.title.y = element_blank()
        )+
        theme(#axis.text.x = element_blank(), 
          # axis.ticks.x = element_blank(), 
          axis.title.x = element_blank()
        )
    }
  }
  for (i in (length(list)-(ncol+2)):length(list)) {
    list[[i]] <- list[[i]]+
      theme(#axis.text.y = element_blank(), 
        #axis.ticks.y = element_blank(), 
        axis.title.y = element_blank()
      )
  }
  return(list)
}
ggarrange(plotlist = tidyfunction(test,ncol = 5),ncol = 5,nrow = 5)

### test code ###
plot_data <- lapply(exp_group_object@deg, function(.x,.y){
  .x@deg_enrich
})
GO_up_list <- lapply(plot_data, function(.x){return(.x$GO_UP)})
GO_down_list <- lapply(plot_data, function(.x){return(.x$GO_DOWN)})
KEGG_up_list <- lapply(plot_data, function(.x){return(.x$KEGG_UP)})
KEGG_down_list <- lapply(plot_data, function(.x){return(.x$KEGG_DOWN)})
test <- mapply(function(.x,.y){
  data <- .x@result
  data$pvalue2 <- -log10(data$pvalue)
  p <- ggplot(data,aes(x = pvalue2,y = NES))+
    geom_point(size = 3.5)+xlab("-Log10(Adjusted_pvalue)")+
    ggtitle(".y")
  return(p)
},GO_up_list,names(GO_up_list))

ggarrange(plotlist = tidyfunction(test,ncol = 5),ncol = 5,nrow = 5)
