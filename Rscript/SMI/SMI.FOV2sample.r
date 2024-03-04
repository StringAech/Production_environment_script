data <- data.frame(
  TMA = c("A1", "A2", "A3"),
  FOV_label = c("111-113", "114-116,336", "117-118,337")
)
SMI.FOV2sample <- function(data, FOV.col = "FOV_label", TMA.col = "TMA"){
  data <- data %>% dplyr::select(TMA = {{TMA.col}}, FOV_label = {{FOV.col}}, everything())
  data1 <- data %>% 
    separate_longer_delim(FOV_label, delim = ",") %>% 
    dplyr::filter(grepl("\\-", FOV_label)) %>% 
    separate(FOV_label, into = c("FOV_start", "FOV_end"), sep = "-") %>%
    mutate(FOV_start = as.integer(FOV_start),
           FOV_end = as.integer(FOV_end)) %>%
    mutate(FOV = map2(FOV_start, FOV_end, seq)) %>%
    unnest(cols = FOV) %>% 
    dplyr::select(-FOV_start, -FOV_end)
  data2 <- data %>% 
    separate_longer_delim(FOV_label, delim = ",") %>% 
    dplyr::filter(!grepl("\\-", FOV_label)) %>% 
    dplyr::select(TMA = TMA, FOV = FOV_label, everything())
  data <- rbind(data1, data2) %>% 
    arrange(TMA) %>% 
    dplyr::select(fov = FOV, sample = TMA, everything())
  return(data)
}


# SMI.FOV2sample(data, FOV.col = "FOV_label", TMA.col = "TMA")

source("../../../loading.package.R")
suppressMessages(library(purrr))
data <- read.xlsx("../YKKY0279-HsRNA-1000-20240129-CosMx SMI平台上机操作实验记录表.xlsx", sheet = 2) %>% 
  setNames(c("yuceID", "TMA", "FOV"))
SMI.FOV2sample(data, FOV.col = "FOV", TMA.col = "TMA") %>% 
  dplyr::select(1, 2) %>% 
  write.csv("../sample2FOV.csv", row.names = F, quote = F)