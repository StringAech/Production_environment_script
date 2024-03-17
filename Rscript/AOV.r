setwd("~/work")
library(dplyr)
library(stringr)
fold <- list.dirs(full.names = F) %>%
  str_extract(".*(?=YKKY)")

system(paste0("find . -type d -name '", "YKKY", "*'", sep = "", "> ./projectID.txt"))

fold %>% str_extract(".*YKKY[0-9]{4}") %>% unique()
system2(command = sprintf("find . -type d -name '%s*'", "YKKY"), stdout = TRUE, stderr = FALSE)
syu

fold <- system2(command = "find", args = ". -type d -name 'YKKY*'", stdout = T) %>% 
  str_extract(".*YKKY[0-9]{4}") %>% unique()
fold %>% str_extract("(?<=\\/)[^\\/]*$")
DIR <- data.frame(penv_dir = fold, 
                  DP_dir = paste0("/mnt/afs/rescro/")
                    
                    
                    
                    str_extract(fold, "(?<=\\/)[^\\/]*$"))


head(iris)
exp_group_object@anno


exp <- exp_group_object@exp %>% t
exp[1:3,1:3]
anno <- exp_group_object@anno %>% 
  dplyr::select(Group)
data <- merge(exp, anno, by = 0)
head(data)
data$Group <- as.factor(data$Group)
data <- data[, -1]
group_factors <- as.factor(data[, "Group"]) # 第一列作为因子
numeric_data <- data %>% dplyr::select(-Group) # 其余列作为数值型数据
for(i in seq_along(numeric_data)) {
  current_column <- numeric_data[, i]
  
  # 创建模型
  model <- aov(current_column ~ group_factors)
  anova_results <- anova(model)
  
}


res <- apply(dplyr::select(data, -Group) %>% select(1:3), 2, \(.x) aov(formula = .x ~ data$Group, data))
