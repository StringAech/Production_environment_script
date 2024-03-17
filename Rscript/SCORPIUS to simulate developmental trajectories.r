install.packages("SCORPIUS")
library(SCORPIUS)
data(ginhoux)
ginhoux
head(ginhoux$sample_info)
colnames(ginhoux$expression)
rownames(ginhoux$expression)

group_names <- exp_group_object@anno %>% 
  dplyr::select(CancerType)
exp <- merge(group_names, exp_group_object@exp %>% t, by = 0) %>% 
  dplyr::select(-1,-2) %>% 
  as.matrix()
group_names <- group_names$CancerType %>% factor(levels = time.point)
space <- reduce_dimensionality(exp, "spearman")
draw_trajectory_plot(space, group_names, contour = TRUE)
traj <- infer_trajectory(space)
draw_trajectory_plot(space, group_names, traj$path, contour = TRUE)
