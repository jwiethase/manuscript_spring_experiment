library(tidyverse)
library(ggnewscale)
library(vegan)

theme_set(theme_bw(base_size = 12))

# Load data =============
treatment_check <- read_csv('data/Treatment_check.csv') %>% 
      filter(Treatment != 'out') 

data_communities <- read_csv('data/biom_pg_L.csv') %>% 
      pivot_longer(cols = 7:43, names_to = 'tank', values_to = 'biomass') %>% 
      group_by(tank) %>% 
      summarise(richness = sum(biomass > 0),
                shannon = diversity(biomass),
                pielou = shannon/log(richness)) %>% 
      rowwise() %>% 
      mutate(Mesocosm = strsplit(tank, '_')[[1]][1],
             Treatment = gsub('_', '', sub("^[^_]+_(.*)_[^_]+$", "\\1", tank)),
             Sampling = as.integer(sub(".*_(.+)$", "\\1", tank)),
             Treatment = ifelse(Treatment == 'controll', 'control', Treatment),
             Temperature = ifelse(Mesocosm %in% LETTERS[1:6], 'cold', 'warm'),
             Mesocosm = ifelse(Treatment == 'control' | Treatment == 'ctrl', 'Control', Mesocosm)) %>% 
      filter(Treatment != 'outside',
             Sampling != 10) %>% 
      select(-tank, -Treatment)

# Check: Do communities in tanks differ? =============
check_data <- treatment_check %>% 
      mutate(Mesocosm = ifelse(Treatment == 'control', 'Control', Mesocosm),
             Mesocosm_ID = paste(Mesocosm, Temperature, sep = '_'),
             isControl = ifelse(Mesocosm == 'Control', 'Control', 'NonControl')) %>% 
      filter(Sampling != 10) %>% 
      left_join(data_communities) %>% 
      pivot_longer(cols = c(Diatom_biovol:dino_perc_of_totcells, richness:pielou), 
                   names_to = 'com_metric', values_to = 'value') %>% 
      filter(com_metric %in% c('cellratio_dd', 'dia_perc_of_totcells', 'richness', 'pielou')) 
check_data$Sampling[check_data$Sampling == 3] <- 6
check_data$Sampling[check_data$Sampling == 1] <- 3
# check_data$Sampling[check_data$Sampling == 10] <- 24

warm_data <- check_data %>%
      filter(Temperature == 'warm')
cold_data <- check_data %>%
      filter(Temperature == 'cold') %>% 
      mutate(Mesocosm_ID = relevel(factor(Mesocosm_ID), ref = 'Control_cold'))

dummy_data <- expand.grid(
      com_metric = c('dia_perc_of_totcells'),
      Temperature = unique(check_data$Temperature),
      Sampling = 3) %>%
      mutate(value = c(0, 100))

t.test(check_data$value[check_data$Treatment != 'control' & check_data$com_metric == 'dia_perc_of_totcells' & check_data$Sampling == 3 & check_data$Temperature == 'warm'], 
       mu = check_data$value[check_data$Treatment == 'control' & check_data$com_metric == 'dia_perc_of_totcells' & check_data$Sampling == 3 & check_data$Temperature == 'warm'])

pvals <- check_data %>%
   group_by(com_metric, Temperature, Sampling) %>%
   summarise(
      mu0 = mean(value[Treatment=='control'], na.rm = TRUE),
      p.value = round(t.test(value[Treatment!='control'], mu = mu0)$p.value, digits = 3),
      y.max = max(value, na.rm = TRUE),
      .groups = 'drop'
   ) %>%
   mutate(
      label = paste0("p = ", signif(p.value, 2))
   )


check_plot <- ggplot() +
      geom_jitter(data = warm_data %>% filter(com_metric %in% c('cellratio_dd', 'dia_perc_of_totcells')),
                  aes(x = Sampling,
                      y = value,
                      fill = Mesocosm_ID,
                      shape = isControl,
                      size = isControl, 
                      alpha = isControl),
                  color = 'black', width = 0.1) +
      scale_fill_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds',
                        guide = guide_legend(override.aes = list(shape = 21, size = 5, color = 'black'))) +
      scale_shape_manual(values = c('Control' = 23, 'NonControl' = 21), 
                         name = element_blank(),
                         guide = guide_legend(override.aes = list(size = 5))) +
      scale_size_manual(values = c('Control' = 5, 'NonControl' = 3)) +
      scale_alpha_manual(values = c('Control' = 1, 'NonControl' = 0.5)) +
      guides(size = 'none', alpha = 'none') + 
      new_scale_fill() + 
      geom_jitter(data = cold_data %>% filter(com_metric %in% c('cellratio_dd', 'dia_perc_of_totcells')),
                  aes(x = Sampling,
                      y = value,
                      fill = Mesocosm_ID,
                      shape = isControl,
                      size = isControl, 
                      alpha = isControl),
                  color = 'black', width = 0.1) +
      scale_fill_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues',
                        guide = guide_legend(override.aes = list(shape = 21, size = 5, color = 'black'))) +
      scale_shape_manual(values = c('Control' = 23, 'NonControl' = 21), 
                         name = element_blank(),
                         guide = guide_legend(override.aes = list(size = 5))) +
      scale_size_manual(values = c('Control' = 5, 'NonControl' = 3)) +
      scale_alpha_manual(values = c('Control' = 1, 'NonControl' = 0.5)) +
      guides(size = 'none', alpha = 'none') +
      geom_blank(data = dummy_data, aes(x = Sampling, y = value)) +
      facet_grid(com_metric ~ Temperature, 
                 scales = 'free', 
                 labeller = as_labeller(c('cellratio_dd' = 'Dia/Din Cellratio',
                                          'dia_perc_of_totcells' = '% Diatoms',
                                          # 'pielou' = "Pielou's Evenness",
                                          # 'richness' = 'Species Richness',
                                          'warm' = '10°C',
                                          'cold' = '5°C'))) +
      ylab(NULL) +
      scale_x_continuous(breaks = seq(1, 6, 1), name = 'Day of study', expand = expansion(0.1)) +
      scale_y_continuous(expand = expansion(0.15)) +
      theme(legend.key.size = unit(1,"line")) +
      geom_text(
         data = pvals %>% filter(Sampling == 3,
                                 com_metric %in% c('cellratio_dd', 'dia_perc_of_totcells')),
         aes(x = Sampling, y = y.max*1.1, label = label), 
         inherit.aes = FALSE, size = 3, nudge_x = 0.5) +
   geom_text(
      data = pvals %>% filter(Sampling == 6,
                              com_metric %in% c('cellratio_dd', 'dia_perc_of_totcells')),
      aes(x = Sampling, y = y.max*1.1, label = label), 
      inherit.aes = FALSE, size = 3, nudge_x = -0.5); check_plot

ggsave('figures/check_plot.jpg', check_plot, width = 7, height = 5, units = 'in')

