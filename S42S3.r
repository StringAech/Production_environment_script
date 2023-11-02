nofilter_exp <- readRDS("~/work/DSP/CTA/YTKY0022_北大人民/DSP_CTA/results/nofilter_exp.Rds")
slot_names <- slotNames(nofilter_exp)
S3_nofilter_exp <- list()
for (slot_name in slot_names) {
  S3_nofilter_exp[[slot_name]] <- slot(nofilter_exp , slot_name)
}

filter_exp <- readRDS("~/work/DSP/CTA/YTKY0022_北大人民/DSP_CTA/results/filter_exp.Rds")
slot_names <- slotNames(filter_exp)
S3_filter_exp <- list()
for (slot_name in slot_names) {
  S3_filter_exp[[slot_name]] <- slot(filter_exp , slot_name)
}

saveRDS(S3_filter_exp,"./S3_filter_exp.Rds")


slot_names <- slotNames(exp_group_object)
S3_res <- list()
for (slot_name in slot_names) {
  S3_res[[slot_name]] <- slot(exp_group_object , slot_name)
}

results <- S3_res$deg
saveRDS(results,"./S3_deg.rds")

test <- lapply(results, function(.x){
  slot_names <- slotNames(.x)
  data <- list()
  for (slot_name in slot_names) {
    data[[slot_name]] <- slot(.x , slot_name)
  }
  return(data)
})

saveRDS(test,"./S3_deg2.rds")
