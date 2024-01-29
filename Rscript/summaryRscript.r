#!/path/to/Rscript
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
spec <-  matrix(c('RDs_path', 'p', 1, "character", "need to input the DEG-calculation path",
                  'out_dir', 'o', 1, "character", "output file path"), byrow=TRUE, ncol=5)
opt <-  getopt(spec)
if ( is.null(opt$RDs_path) ) { opt$RDs_path  = list.files(path = "./", recursive = T, 
                                                          pattern = ".*DEG-cal.*.RData", 
                                                          full.names = F)  }
if ( is.null(opt$out_dir)  ) { opt$out_dir   = "./"    }

# if( length(opt$RDs_path) == 0 ){
#   message("\n")
#   cat("\033[31mThe related RDs file is not found in the current path! \n\033[0m\n")
#   # message("The related RDs file is not found in the current path! \n")
#   cat(paste(getopt(spec=spec, usage = T), "\n"))
#   quit()
# }
message("\n")
stop_function <- function(){stop()}
tryCatch({stop_function()}, error = function(e){
  # message("Please confirm whether the path of RData is correct, It will hold for 10 seconds: \n", 
  #         opt$RDs_path)
  # cat("\033[32mThis text is green.\033[0m\n")
  cat(paste0("\033[32mPlease confirm whether the path of RData is correct, It will hold for 10 seconds: \n",
             getwd(), "/", "\033[31m",opt$RDs_path, "\n\033[0m\n","\033[0m\n"))
  # response <- readline()
  # if (response == "n") {stop("Code execution stopped by user.")}
  Sys.sleep(10)
})
# message("\n")
load(opt$RDs_path)
DEG_group <- lapply(exp_group_object@deg, \(.x) {
  log <- .x@result_df %>% dplyr::filter(change != "NOT") %>% nrow %>% as.logical()
  return(if (log) { .x@label_string %>% gsub("的差异", "", .)})
}) %>% do.call(rbind, .) %>% as.character()
no_DEG_group <- lapply(exp_group_object@deg, \(.x) {
  log <- .x@result_df %>% dplyr::filter(change != "NOT") %>% nrow %>% as.logical()
  return(if (!log) { .x@label_string %>% gsub("的差异", "", .)})
}) %>% do.call(rbind, .) %>% as.character()
list <- list(DEG_group, no_DEG_group)
lapply(lapply(list, unlist), `length<-`, max(lengths(list))) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  setNames(c("DEG_group", "no_DEG_group")) %>%
  write.table(file = paste0(opt$out_dir, "/deg_group.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
# message("\n")
message("Start copying files to the corresponding folder....")
