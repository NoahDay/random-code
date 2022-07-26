pacman::p_load(GGally,tidymodels,tidyverse,janitor,dplyr,ggplot2)

# Looking at arrests data for this example
data("USArrests")

# Parallel coordinate plot
# Gives us an idea of the variance of each variable (we see assault counts are highly variable compared to the other crimes)
ggparcoord(USArrests, scale = "globalminmax")

# Set up PCA using recipes
US_arrest_PCA <-
  recipe(~., USArrests) %>%
  step_normalize(all_numeric()) %>% 
  step_pca(all_numeric(), num_comp = 4)


# Prep it
arrests_prep <- prep(US_arrest_PCA)


# Loadings plot
# Loadings are the coefficients of the linear combination of the original variables from which the PCs are constructed
pacman::p_load(tidytext)
arrests_prep %>% 
  tidy(n = 2) %>% 
  filter(component %in% c("PC1", "PC2", "PC3", "PC4")) %>% 
  group_by(component) %>% 
  top_n(5) %>% 
  mutate(
    terms_order = reorder_within(terms, abs(value), component)
  ) %>% 
  ggplot(
    aes(value, terms_order, fill = terms) 
  ) + 
  geom_col(show.legend = FALSE) + 
  facet_wrap(~component, scales = "free_y") + 
  scale_y_reordered()



# Scree plot
# This tells us the percentage of the variance that is explained by each PC
sdev <- arrests_prep$steps[[2]]$res$sdev

percent_variation <- sdev^2 / sum(sdev^2)
percent_variation

tibble(
  component = str_c("PC", 1:length(percent_variation)),
  percent_var = percent_variation
) %>%
  mutate(component = fct_inorder(component)) %>%
  ggplot(aes(component, percent_var)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Percent variance explained by each PCA component")


