#### upset plot 
install.packages("UpSetR")
devtools::install_github("hms-dbmi/UpSetR")
devtools::install_github('krassowski/complex-upset')
library(UpSetR)
library(ComplexUpset)
# listInput <- list(
#   one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), 
#   two = c(1, 2, 4, 5, 10), 
#   three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))
# p <- upset(fromList(listInput), order.by = "freq")

suppressPackageStartupMessages(library(EPARS))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggsci))

load("~/work/Olink/YKKY0100/Olink_plasma/results_plasma/rendering_cache_Olink/DEG-calculation_55d5c4ba314cbbb7e3ee0dca4a9def16.RData")
names(exp_group_object@deg)

exp_list <- exp_group_object@deg %>% 
  lapply(\(.x){.x@result_df}) %>% 
  lapply(\(.x){
    .x %>% dplyr::filter(Change != "NOT") %>% rownames()
  })

a <- fromList(exp_list)
rownames(a) <- exp_list %>% unlist %>% unique()
p2 <- upset(a,colnames(a))


gene_pathway_membership <- a %>% apply(2, as.logical) %>% t


tidy_pathway_member <- gene_pathway_membership %>%
  as_tibble(rownames = "Pathway") %>%
  gather(Gene, Member, -Pathway) %>%
  filter(Member) %>%
  select(- Member) %>%
  group_by(Gene) %>%
  summarize(Pathways = list(Pathway))

p <- tidy_pathway_member %>%
  ggplot(aes(x = Pathways)) +
  geom_bar() +
  scale_x_upset(
    order_by = "degree",
    reverse = T, 
  )+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) + 
  scale_y_continuous(name = "Intersection size",expand = c(0, 0), limits = c(0, 6))+
  theme(text = element_text(family = "ARIAL", colour = "black", size = 9))+
  xlab('')+
  theme_combmatrix(
    combmatrix.label.text = element_text(family = "ARIAL", colour = "black", size = 9),
  )+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))+
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 9)
  )

ggsave("./upsets_plot/test/test.pdf", plot = p, width = 12)
#### CSF ####
load("~/work/Olink/YKKY0100/Olink_CSF/results_CSF/rendering_cache_Olink/DEG-calculation_721c515ef58d91962b9a95096dd3b72b.RData")
exp_list <- exp_group_object@deg %>% 
  lapply(\(.x){.x@result_df}) %>% 
  lapply(\(.x){
    .x %>% dplyr::filter(Change != "NOT") %>% rownames()
  })

a <- fromList(exp_list)
rownames(a) <- exp_list %>% unlist %>% unique()
p2 <- upset(a,colnames(a))


gene_pathway_membership <- a %>% apply(2, as.logical) %>% t


tidy_pathway_member <- gene_pathway_membership %>%
  as_tibble(rownames = "Pathway") %>%
  gather(Gene, Member, -Pathway) %>%
  filter(Member) %>%
  select(- Member) %>%
  group_by(Gene) %>%
  summarize(Pathways = list(Pathway))

p <- tidy_pathway_member %>%
  ggplot(aes(x = Pathways)) +
  geom_bar() +
  scale_x_upset(
    order_by = "degree",
    reverse = T, 
  )+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) + 
  scale_y_continuous(name = "Intersection size",expand = c(0, 0), limits = c(0, 6))+
  theme(text = element_text(family = "ARIAL", colour = "black", size = 9))+
  xlab('')+
  theme_combmatrix(
    combmatrix.label.text = element_text(family = "ARIAL", colour = "black", size = 9),
  )+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))+
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 9)
  )

ggsave("./upsets_plot/test/CSF_upset.pdf", plot = p, width = 12)


