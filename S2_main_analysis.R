library(tidyverse)
library(nlme)
library(lme4)
library(ggpubr)
library(DHARMa)
library(ggh4x)
library(patchwork)
library(ggnewscale)
library(writexl)
source('source/misc_functions.R')
theme_set(theme_bw(base_size = 12))

# L&O figure requirements:
# Font size greater than 8 pt at 300+ dpi



# Load data =============
data <- read_csv("data/SpringExp23_predictor_response.csv") %>% 
      filter(Treatment != 'out', !is.na(POC_ug_L)) %>% 
      mutate(Treatment = as.factor(Treatment),
             Temperature = as.factor(Temperature),
             Mesocosm_ID = as.factor(paste(Mesocosm, Temperature, sep = '_')),
             POC_chla = POC_umol_L / Chl_umol_L) 

data_sediment <- read_csv("data/SpringExp23_predictor_response_sedtrap.csv") %>% 
   mutate(Treatment = as.factor(Treatment),
          Temperature = as.factor(Temperature),
          Mesocosm_ID = as.factor(paste(Mesocosm, Temperature, sep = '_')))

data_GHG <- read_csv("data/SpringExp23_predictor_response_GHG.csv") %>% 
   mutate(Treatment = as.factor(Treatment),
          Temperature = as.factor(Temperature),
          Mesocosm_ID = as.factor(paste(Mesocosm, Temperature, sep = '_')),
          CO2_mmol_m2_d = CO2_mol_m2_d*1000) %>% 
   left_join(data %>% select(DOY, Sample_ID, Chl_umol_L)) 

# POC ===================
## Data exploration =============
warm_data <- data %>%
   filter(Temperature == 'warm')
cold_data <- data %>%
   filter(Temperature == 'cold')

warm_curves <- ggplot() +
   geom_point(data = warm_data,
               aes(x = DOY,
                   y = POC_umol_L,
                   fill = Mesocosm_ID),
               color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = warm_data,
              aes(x = DOY,
                  y = POC_umol_L,
                  col = Mesocosm_ID),
              method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab(NULL) +
   ylim(0, 250) + 
   theme(legend.position = 'none'); warm_curves

cold_curves <- ggplot() +
   geom_point(data = cold_data,
              aes(x = DOY,
                  y = POC_umol_L,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = cold_data,
               aes(x = DOY,
                   y = POC_umol_L,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab('POC (µmol L^-1)') +
   ylim(0, 250) + 
   theme(legend.position = 'none'); cold_curves

combined_curves <- ggpubr::ggarrange(cold_curves, warm_curves); combined_curves
ggsave(combined_curves, filename = 'figures/combined_curves_POC.jpg', width = 12)


## Overall temporal trend of POC =============
model <- lmer(POC_umol_L ~ Temperature * poly(DOY, degree = 2) + (1|Mesocosm), 
                               data = data)

simulationOutput <- simulateResiduals(fittedModel = model)
testResiduals(simulationOutput)
plot(simulationOutput)

simulationOutput_2 <- recalculateResiduals(simulationOutput, group = data$DOY) 
testTemporalAutocorrelation(simulationOutput_2, time = unique(data$DOY))

# Very little temporal autocorrelation, meaning that the polynomial term of DOY
# has soaked up the temporal autocorrelation
summary_POC <- as.data.frame(summary(model)[['coefficients']]) %>% 
   mutate(response = 'POC') %>% 
   rownames_to_column('coefficient')
anova_POC <- car::Anova(model) %>% 
   mutate(response = 'POC') %>% 
   rownames_to_column('coefficient')

### Plot model results =============
# Fitted DOY effect
newdata_DOY <- expand.grid(
   DOY = seq(min(data$DOY), max(data$DOY), length.out = 100),
   Temperature = c('cold', 'warm')
)

preds_DOY <- predict(model, newdata = newdata_DOY, se.fit = TRUE, re.form = NA)

newdata_DOY_POC <- newdata_DOY %>%
   mutate(
      fit  = preds_DOY$fit,
      lower = fit - 1.96 * preds_DOY$se,
      upper = fit + 1.96 * preds_DOY$se
   )

jpeg('figures/model_result_POC.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot() +
   geom_point(data = data, aes(x = DOY, y = POC_umol_L, col = Temperature), cex = 2, alpha = 0.3) +
   geom_ribbon(data = newdata_DOY_POC, aes(x = DOY, ymin = lower, ymax = upper, fill = Temperature), alpha = 0.3) +
   geom_line(data = newdata_DOY_POC, aes(x = DOY, y = fit, col = Temperature), size = 1.2) +
   scale_color_manual(values = c('blue', 'red')) +
   scale_fill_manual(values = c('blue', 'red')) +
   labs(x = "DOY",
        y = "Predicted POC (µmol L^-1)")
dev.off()

newdata_DOY_POC$response <- 'Water Column POC'

## Temporal trajectory analysis of POC =============
### Aggregate the data =============
summary_data <- data %>%
   group_by(Temperature, Treatment) %>%
   group_split() %>%
   map_dfr(summarize_tank, var = 'POC_umol_L')

par(mfrow = c(1, 1))
plot_correlation(summary_data, 3:NCOL(summary_data))

jpeg('figures/covar_corr_POC.jpg', res = 200, width = 10, height = 7, units = 'in')
plot_correlation(summary_data, c("peak_raw", "time_to_peak_raw", "cum_pred_fitted"))
dev.off()

### Fit models =============
summary_results <- summary_data %>%
   # Select those metrics that are not highly correlated with each other
   select(Temperature, peak_fit, time_to_peak_fit, cum_pred_fitted) %>%   
   pivot_longer(-c(Temperature), names_to = "metric", values_to = "value") %>%
   group_by(metric) %>%
   summarize(model = list(lm(value ~ Temperature, data = cur_data())), .groups = "drop") %>%
   mutate(summary = map(model, broom::tidy)) %>%
   unnest(summary) %>%
   # Apply Benjamini–Hochberg (BH) p-value adjustment for multiple testing. Slightly less
   # conservative than Bonferroni, maintaining more power
   mutate(adj_p = round(p.adjust(p.value, method = "BH"), 3),
          p.value = round(p.value, 3),
          lower = estimate - 1.96 * std.error,
          upper = estimate + 1.96 * std.error) %>% 
   ungroup()

### Plot the results =============
temp_levels <- sort(unique(summary_data$Temperature))

# For each metric, fit a model and compute predictions (with standard errors)
predictions <- map_dfr(c("cum_pred_fitted", "peak_fit", "time_to_peak_fit"), function(m) {
   mod <- lm(as.formula(paste(m, "~ Temperature")), data = summary_data)
   newdata <- tibble(Temperature = temp_levels)

   pred <- predict(mod, newdata, se.fit = TRUE)
   
   tibble(
      metric = m,
      Temperature = newdata$Temperature,
      estimate = pred$fit,
      lower = pred$fit - 1.96 * pred$se.fit,
      upper = pred$fit + 1.96 * pred$se.fit
   )
})

p_values <- summary_results %>%
   filter(term == "Temperaturewarm") %>%
   select(metric, adj_p) %>% 
   mutate(adj_p = ifelse(adj_p < 0.001, 'p < 0.001', paste0('p = ', adj_p)))

# Plot predictions for each metric using ggplot2
jpeg('figures/temp_result_POC.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot(predictions, aes(x = Temperature, y = estimate, group = metric)) +
   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
   geom_point(pch = 23, size = 5, fill = 'coral') +
   facet_wrap(~ metric, scales = "free_y") +
   labs(
      x = "Temperature",
      y = "Predicted value"
   ) +
   geom_text(
      data = p_values,
      aes(x = 0.5, y = -Inf, label = adj_p),
      hjust = 0, vjust = -1,
      size = 4, color = "black",
      inherit.aes = FALSE
   ) +
   scale_y_continuous(expand = expansion(0.2))
dev.off()

# DOC ===================
## Data exploration =============
data_DOC <- data %>% 
   filter(!is.na(DOC_umol_L))
warm_data_DOC <- data_DOC %>%
   filter(Temperature == 'warm')
cold_data_DOC <- data_DOC %>%
   filter(Temperature == 'cold')

warm_curves_DOC <- ggplot() +
   geom_point(data = warm_data_DOC,
              aes(x = DOY,
                  y = DOC_umol_L,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = warm_data_DOC,
               aes(x = DOY,
                   y = DOC_umol_L,
                   col = Mesocosm_ID),
               method = "lm", se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab(NULL) +
   ylim(400, 700) + 
   theme(legend.position = 'none'); warm_curves_DOC

cold_curves_DOC <- ggplot() +
   geom_point(data = cold_data_DOC,
              aes(x = DOY,
                  y = DOC_umol_L,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = cold_data_DOC,
               aes(x = DOY,
                   y = DOC_umol_L,
                   col = Mesocosm_ID),
               method = "lm", se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab('DOC (µmol L^-1)') +
   ylim(400, 700) + 
   theme(legend.position = 'none'); cold_curves_DOC

combined_curves_DOC <- ggpubr::ggarrange(cold_curves_DOC, warm_curves_DOC)
ggsave(combined_curves_DOC, filename = 'figures/combined_curves_DOC.jpg', width = 12)

## Overall temporal trend of DOC =============
model_DOC <- lmer(log(DOC_umol_L) ~ Temperature * poly(DOY, degree = 2) + (1|Mesocosm), 
                      data = data_DOC)

simulationOutput <- simulateResiduals(fittedModel = model_DOC)
testResiduals(simulationOutput)
plot(simulationOutput)

summary_DOC <- as.data.frame(summary(model_DOC)[['coefficients']]) %>% 
   mutate(response = 'DOC') %>% 
   rownames_to_column('coefficient')
anova_DOC <- car::Anova(model_DOC) %>% 
   mutate(response = 'DOC') %>% 
   rownames_to_column('coefficient')

### Plot model results =============
# Fitted DOY effect
newdata_DOY_DOC <- expand.grid(
   DOY = seq(min(data_DOC$DOY), max(data_DOC$DOY), length.out = 100),
   Temperature = c('cold', 'warm')
)

preds_DOY <- predict(model_DOC, newdata = newdata_DOY_DOC, se.fit = TRUE, re.form = NA)

newdata_DOY_DOC <- newdata_DOY_DOC %>%
   mutate(
      fit    = exp(preds_DOY$fit),
      lower  = exp(preds_DOY$fit - 1.96 * preds_DOY$se.fit),
      upper  = exp(preds_DOY$fit + 1.96 * preds_DOY$se.fit)
   )

jpeg('figures/model_result_DOC.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot() +
   geom_point(data = data_DOC, aes(x = DOY, y = DOC_umol_L, col = Temperature), cex = 2, alpha = 0.3) +
   geom_ribbon(data = newdata_DOY_DOC, aes(x = DOY, ymin = lower, ymax = upper, fill = Temperature), alpha = 0.3) +
   geom_line(data = newdata_DOY_DOC, aes(x = DOY, y = fit, col = Temperature), size = 1.2) +
   scale_color_manual(values = c('blue', 'red')) +
   scale_fill_manual(values = c('blue', 'red')) +
   labs(x = "DOY",
        y = "Predicted DOC (µmol L^-1)")
dev.off()

newdata_DOY_DOC$response <- 'Water Column DOC'

## Temporal trajectory analysis of DOC =============
### Aggregate the data =============
summary_data_DOC <- data_DOC %>%
   group_by(Temperature, Treatment) %>%
   group_split() %>%
   map_dfr(summarize_tank, var = 'DOC_umol_L')

### Fit models =============
summary_results_DOC <- summary_data_DOC %>%
   select(Temperature, peak_fit, time_to_peak_fit, cum_pred_fitted) %>%   
   pivot_longer(-c(Temperature), names_to = "metric", values_to = "value") %>%
   group_by(metric) %>%
   summarize(model = list(lm(value ~ Temperature, data = cur_data())), .groups = "drop") %>%
   mutate(summary = map(model, broom::tidy)) %>%
   unnest(summary) %>% 
   mutate(p.value = round(p.value, 3))

### Plot the results =============
predictions_DOC <- map_dfr(c("cum_pred_fitted", "peak_fit", "time_to_peak_fit"), function(m) {
   mod <- lm(as.formula(paste(m, "~ Temperature")), data = summary_data_DOC)
   newdata <- tibble(Temperature = temp_levels)
   
   pred <- predict(mod, newdata, se.fit = TRUE)
   
   tibble(
      metric = m,
      Temperature = newdata$Temperature,
      estimate = pred$fit,
      lower = pred$fit - 1.96 * pred$se.fit,
      upper = pred$fit + 1.96 * pred$se.fit
   )
})

p_values_DOC <- summary_results_DOC %>%
   filter(term == "Temperaturewarm") %>%
   select(metric, p.value) %>% 
   mutate(p.value = ifelse(p.value < 0.001, 'p < 0.001', paste0('p = ', p.value)))

jpeg('figures/temp_result_DOC.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot(predictions_DOC, aes(x = Temperature, y = estimate, group = metric)) +
   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
   geom_point(pch = 23, size = 5, fill = 'coral') +
   facet_wrap(~ metric, scales = "free_y") +
   labs(
      x = "Temperature",
      y = "Predicted value"
   ) +
   geom_text(
      data = p_values_DOC,
      aes(x = 0.5, y = -Inf, label = p.value),
      hjust = 0, vjust = -1,
      size = 4, color = "black",
      inherit.aes = FALSE
   ) +
   scale_y_continuous(expand = expansion(0.2))
dev.off()


# Sediment POC ===================
## Data exploration =============
warm_data_sed <- data_sediment %>%
   filter(Temperature == 'warm')
cold_data_sed <- data_sediment %>%
   filter(Temperature == 'cold')

warm_curves_sed_POC <- ggplot() +
   geom_point(data = warm_data_sed,
              aes(x = DOY,
                  y = POC_mmol_m2_d,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = warm_data_sed,
               aes(x = DOY,
                   y = POC_mmol_m2_d,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab(NULL) +
   ylim(min(data_sediment$POC_mmol_m2_d), max(data_sediment$POC_mmol_m2_d)) + 
   theme(legend.position = 'none'); warm_curves_sed_POC

cold_curves_sed_POC <- ggplot() +
   geom_point(data = cold_data_sed,
              aes(x = DOY,
                  y = POC_mmol_m2_d,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = cold_data_sed,
               aes(x = DOY,
                   y = POC_mmol_m2_d,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab('POC Sediment [µg m^2 d]') +
   ylim(min(data_sediment$POC_mmol_m2_d), max(data_sediment$POC_mmol_m2_d)) + 
   theme(legend.position = 'none'); cold_curves_sed_POC

combined_curves_sed_POC <- ggpubr::ggarrange(cold_curves_sed_POC, warm_curves_sed_POC)
ggsave(combined_curves_sed_POC, filename = 'figures/combined_curves_sed_POC.jpg', width = 12)

## Overall temporal trend of sediment POC =============
model_sed_POC <- lmer(log(POC_mmol_m2_d) ~ Temperature * poly(DOY, degree = 2) + (1|Mesocosm), 
                  data = data_sediment) 
# Singular fit, this could happen if there was very little variation between Mesocosms.
simulationOutput <- simulateResiduals(fittedModel = model_sed_POC)
testResiduals(simulationOutput)
plot(simulationOutput)

summary_sedPOC <- as.data.frame(summary(model_sed_POC)[['coefficients']]) %>% 
   mutate(response = 'Sediment Trap POC') %>% 
   rownames_to_column('coefficient')
anova_sedPOC <- car::Anova(model_sed_POC) %>% 
   mutate(response = 'Sediment Trap POC') %>% 
   rownames_to_column('coefficient')

### Plot model results =============
# Fitted DOY effect
newdata_DOY_sed <- expand.grid(
   DOY = seq(min(data_sediment$DOY), max(data_sediment$DOY), length.out = 100),
   Temperature = c('cold', 'warm')
)

preds_sedPOC_DOY <- predict(model_sed_POC, newdata = newdata_DOY_sed, se.fit = TRUE, re.form = NA)

newdata_sedPOC_DOY <- newdata_DOY_sed %>%
   mutate(
      fit    = exp(preds_sedPOC_DOY$fit),
      lower  = exp(preds_sedPOC_DOY$fit - 1.96 * preds_sedPOC_DOY$se.fit),
      upper  = exp(preds_sedPOC_DOY$fit + 1.96 * preds_sedPOC_DOY$se.fit)
   )

jpeg('figures/model_result_sedPOC.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot() +
   geom_point(data = data_sediment, aes(x = DOY, y = POC_mmol_m2_d, col = Temperature), cex = 2, alpha = 0.3) +
   geom_ribbon(data = newdata_sedPOC_DOY, aes(x = DOY, ymin = lower, ymax = upper, fill = Temperature), alpha = 0.3) +
   geom_line(data = newdata_sedPOC_DOY, aes(x = DOY, y = fit, col = Temperature), size = 1.2) +
   scale_color_manual(values = c('blue', 'red')) +
   scale_fill_manual(values = c('blue', 'red')) +
   labs(x = "DOY",
        y = "Predicted DOC (µmol L^-1)")
dev.off()

newdata_sedPOC_DOY$response <- 'Sediment Trap POC'

## Temporal trajectory analysis of sediment POC =============
### Aggregate the data =============
summary_data_sedPOC <- data_sediment %>%
   group_by(Temperature, Treatment) %>%
   group_split() %>%
   map_dfr(summarize_tank, var = 'POC_mmol_m2_d')

### Fit models =============
summary_results_sedPOC <- summary_data_sedPOC %>%
   select(Temperature, peak_fit, time_to_peak_fit, cum_pred_fitted) %>%   
   pivot_longer(-c(Temperature), names_to = "metric", values_to = "value") %>%
   group_by(metric) %>%
   summarize(model = list(lm(value ~ Temperature, data = cur_data())), .groups = "drop") %>%
   mutate(summary = map(model, broom::tidy)) %>%
   unnest(summary) %>% 
   mutate(p.value = round(p.value, 3))

### Plot the results =============
predictions_sedPOC <- map_dfr(c("cum_pred_fitted", "peak_fit", "time_to_peak_fit"), function(m) {
   mod <- lm(as.formula(paste(m, "~ Temperature")), data = summary_data_sedPOC)
   newdata <- tibble(Temperature = temp_levels)
   
   pred <- predict(mod, newdata, se.fit = TRUE)
   
   tibble(
      metric = m,
      Temperature = newdata$Temperature,
      estimate = pred$fit,
      lower = pred$fit - 1.96 * pred$se.fit,
      upper = pred$fit + 1.96 * pred$se.fit
   )
})

p_values_sedPOC <- summary_results_sedPOC %>%
   filter(term == "Temperaturewarm") %>%
   select(metric, p.value) %>% 
   mutate(p.value = ifelse(p.value < 0.001, 'p < 0.001', paste0('p = ', p.value)))

jpeg('figures/temp_result_sedPOC.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot(predictions_sedPOC, aes(x = Temperature, y = estimate, group = metric)) +
   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
   geom_point(pch = 23, size = 5, fill = 'coral') +
   facet_wrap(~ metric, scales = "free_y") +
   labs(
      x = "Temperature",
      y = "Predicted value"
   ) +
   geom_text(
      data = p_values_sedPOC,
      aes(x = 0.5, y = -Inf, label = p.value),
      hjust = 0, vjust = -1,
      size = 4, color = "black",
      inherit.aes = FALSE
   ) +
   scale_y_continuous(expand = expansion(0.2))
dev.off()

# Sediment DOC ===================
## Data exploration =============
warm_curves_sed_DOC <- ggplot() +
   geom_point(data = warm_data_sed,
              aes(x = DOY,
                  y = DOC_mmol_m2_d,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = warm_data_sed,
               aes(x = DOY,
                   y = DOC_mmol_m2_d,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab(NULL) +
   ylim(min(data_sediment$DOC_mmol_m2_d), max(data_sediment$DOC_mmol_m2_d)) + 
   theme(legend.position = 'none'); warm_curves_sed_DOC

cold_curves_sed_DOC <- ggplot() +
   geom_point(data = cold_data_sed,
              aes(x = DOY,
                  y = DOC_mmol_m2_d,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = cold_data_sed,
               aes(x = DOY,
                   y = DOC_mmol_m2_d,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab('POC Sediment [µg m^2 d]') +
   # ylim(min(data_sediment$DOC_mmol_m2_d), max(data_sediment$DOC_mmol_m2_d)) + 
   theme(legend.position = 'none'); cold_curves_sed_DOC

combined_curves_sed_DOC <- ggpubr::ggarrange(cold_curves_sed_DOC, warm_curves_sed_DOC)
ggsave(combined_curves_sed_DOC, filename = 'figures/combined_curves_sed_DOC.jpg', width = 12)

## Overall temporal trend of sediment DOC =============
model_sed_DOC <- lmer(DOC_mmol_m2_d ~ Temperature * poly(DOY, degree = 2) + (1|Mesocosm), 
                      data = data_sediment) 
# Singular fit, this could happen if there was very little variation between Mesocosms.
simulationOutput <- simulateResiduals(fittedModel = model_sed_DOC)
testResiduals(simulationOutput)
plot(simulationOutput)

summary_sedDOC <- as.data.frame(summary(model_sed_DOC)[['coefficients']]) %>% 
   mutate(response = 'Sediment Trap DOC') %>% 
   rownames_to_column('coefficient')
anova_sedDOC <- car::Anova(model_sed_DOC) %>% 
   mutate(response = 'Sediment Trap DOC') %>% 
   rownames_to_column('coefficient')


### Plot model results =============
# Fitted DOY effect
preds_sedDOC_DOY <- predict(model_sed_DOC, newdata = newdata_DOY_sed, se.fit = TRUE, re.form = NA)

newdata_sedDOC_DOY <- newdata_DOY_sed %>%
   mutate(
      fit    = (preds_sedDOC_DOY$fit),
      lower  = (preds_sedDOC_DOY$fit - 1.96 * preds_sedDOC_DOY$se.fit),
      upper  = (preds_sedDOC_DOY$fit + 1.96 * preds_sedDOC_DOY$se.fit)
   )

jpeg('figures/model_result_sedDOC.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot() +
   geom_point(data = data_sediment, aes(x = DOY, y = DOC_mmol_m2_d, col = Temperature), cex = 2, alpha = 0.3) +
   geom_ribbon(data = newdata_sedDOC_DOY, aes(x = DOY, ymin = lower, ymax = upper, fill = Temperature), alpha = 0.3) +
   geom_line(data = newdata_sedDOC_DOY, aes(x = DOY, y = fit, col = Temperature), size = 1.2) +
   scale_color_manual(values = c('blue', 'red')) +
   scale_fill_manual(values = c('blue', 'red')) +
   labs(x = "DOY",
        y = "Predicted DOC (µmol L^-1)")
dev.off()

newdata_sedDOC_DOY$response <- 'Sediment Trap DOC'

## Temporal trajectory analysis of sediment DOC =============
### Aggregate the data =============
summary_data_sedDOC <- data_sediment %>%
   group_by(Temperature, Treatment) %>%
   group_split() %>%
   map_dfr(summarize_tank, var = 'DOC_mmol_m2_d')

### Fit models =============
summary_results_sedDOC <- summary_data_sedDOC %>%
   select(Temperature, peak_fit, time_to_peak_fit, cum_pred_fitted) %>%   
   pivot_longer(-c(Temperature), names_to = "metric", values_to = "value") %>%
   group_by(metric) %>%
   summarize(model = list(lm(value ~ Temperature, data = cur_data())), .groups = "drop") %>%
   mutate(summary = map(model, broom::tidy)) %>%
   unnest(summary) %>% 
   mutate(p.value = round(p.value, 3))

### Plot the results =============
predictions_sedDOC <- map_dfr(c("cum_pred_fitted", "peak_fit", "time_to_peak_fit"), function(m) {
   mod <- lm(as.formula(paste(m, "~ Temperature")), data = summary_data_sedDOC)
   newdata <- tibble(Temperature = temp_levels)
   
   pred <- predict(mod, newdata, se.fit = TRUE)
   
   tibble(
      metric = m,
      Temperature = newdata$Temperature,
      estimate = pred$fit,
      lower = pred$fit - 1.96 * pred$se.fit,
      upper = pred$fit + 1.96 * pred$se.fit
   )
})

p_values_sedDOC <- summary_results_sedDOC %>%
   filter(term == "Temperaturewarm") %>%
   select(metric, p.value) %>% 
   mutate(p.value = ifelse(p.value < 0.001, 'p < 0.001', paste0('p = ', p.value)))

jpeg('figures/temp_result_sedDOC.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot(predictions_sedDOC, aes(x = Temperature, y = estimate, group = metric)) +
   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
   geom_point(pch = 23, size = 5, fill = 'coral') +
   facet_wrap(~ metric, scales = "free_y") +
   labs(
      x = "Temperature",
      y = "Predicted value"
   ) +
   geom_text(
      data = p_values_sedDOC,
      aes(x = 0.5, y = -Inf, label = p.value),
      hjust = 0, vjust = -1,
      size = 4, color = "black",
      inherit.aes = FALSE
   ) +
   scale_y_continuous(expand = expansion(0.2))
dev.off()

# GHG ===================
## Data exploration =============
warm_data_GHG <- data_GHG %>%
   filter(Temperature == 'warm')
cold_data_GHG <- data_GHG %>%
   filter(Temperature == 'cold')

warm_curves_GHG <- ggplot() +
   geom_point(data = warm_data_GHG,
              aes(x = DOY,
                  y = CO2_mmol_m2_d,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = warm_data_GHG,
               aes(x = DOY,
                   y = CO2_mmol_m2_d,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab(NULL) +
   ylim(min(data_GHG$CO2_mmol_m2_d), max(data_GHG$CO2_mmol_m2_d)) + 
   theme(legend.position = 'none'); warm_curves_GHG

cold_curves_GHG <- ggplot() +
   geom_point(data = cold_data_GHG,
              aes(x = DOY,
                  y = CO2_mmol_m2_d,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = cold_data_GHG,
               aes(x = DOY,
                   y = CO2_mmol_m2_d,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab('CO2 [mmol m^2 d]') +
   ylim(min(data_GHG$CO2_mmol_m2_d), max(data_GHG$CO2_mmol_m2_d)) + 
   theme(legend.position = 'none'); cold_curves_GHG

combined_curves_GHG <- ggpubr::ggarrange(cold_curves_GHG, warm_curves_GHG)
ggsave(combined_curves_GHG, filename = 'figures/combined_curves_GHG.jpg', width = 12)

## Overall temporal trend of GHG =============
model_GHG <- lmer(CO2_mmol_m2_d ~ Temperature * poly(DOY, degree = 2) + (1|Mesocosm), 
                  data = data_GHG) 

simulationOutput <- simulateResiduals(fittedModel = model_GHG)
testResiduals(simulationOutput)
plot(simulationOutput)

simulationOutput_2 <- recalculateResiduals(simulationOutput, group = data_GHG$DOY) 
testTemporalAutocorrelation(simulationOutput_2, time = unique(data_GHG$DOY))

summary_CO2 <- as.data.frame(summary(model_GHG)[['coefficients']]) %>% 
   mutate(response = 'CO2') %>% 
   rownames_to_column('coefficient')
anova_CO2 <- car::Anova(model_GHG) %>% 
   mutate(response = 'CO2') %>% 
   rownames_to_column('coefficient')

### Plot model results =============
# Fitted DOY effect
newdata_DOY_GHG <- expand.grid(
   DOY = seq(min(data_GHG$DOY), max(data_GHG$DOY), length.out = 100),
   Temperature = c('cold', 'warm')
)

preds_GHG_DOY <- predict(model_GHG, newdata = newdata_DOY_GHG, se.fit = TRUE, re.form = NA)

newdata_GHG_DOY <- newdata_DOY_GHG %>%
   mutate(
      fit    = (preds_GHG_DOY$fit),
      lower  = (preds_GHG_DOY$fit - 1.96 * preds_GHG_DOY$se.fit),
      upper  = (preds_GHG_DOY$fit + 1.96 * preds_GHG_DOY$se.fit)
   )

jpeg('figures/model_result_GHG.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot() +
   geom_point(data = data_GHG, aes(x = DOY, y = CO2_mmol_m2_d, col = Temperature), cex = 2, alpha = 0.3) +
   geom_ribbon(data = newdata_GHG_DOY, aes(x = DOY, ymin = lower, ymax = upper, fill = Temperature), alpha = 0.3) +
   geom_line(data = newdata_GHG_DOY, aes(x = DOY, y = fit, col = Temperature), size = 1.2) +
   scale_color_manual(values = c('blue', 'red')) +
   scale_fill_manual(values = c('blue', 'red')) +
   labs(x = "DOY",
        y = "Predicted CO2 [mmol m^2 d]")
dev.off()

newdata_GHG_DOY$response <- 'CO2'

## Temporal trajectory analysis of GHG =============
### Aggregate the data =============
summary_data_GHG <- data_GHG %>%
   group_by(Temperature, Treatment) %>%
   group_split() %>%
   map_dfr(summarize_tank, var = 'CO2_mmol_m2_d')

### Fit models =============
summary_results_GHG <- summary_data_GHG %>%
   select(Temperature, peak_fit, time_to_peak_fit, cum_pred_fitted) %>%   
   pivot_longer(-c(Temperature), names_to = "metric", values_to = "value") %>%
   group_by(metric) %>%
   summarize(model = list(lm(value ~ Temperature, data = cur_data())), .groups = "drop") %>%
   mutate(summary = map(model, broom::tidy)) %>%
   unnest(summary) %>% 
   mutate(p.value = round(p.value, 3))

### Plot the results =============
predictions_GHG <- map_dfr(c("cum_pred_fitted", "peak_fit", "time_to_peak_fit"), function(m) {
   mod <- lm(as.formula(paste(m, "~ Temperature")), data = summary_data_GHG)
   newdata <- tibble(Temperature = temp_levels)
   
   pred <- predict(mod, newdata, se.fit = TRUE)
   
   tibble(
      metric = m,
      Temperature = newdata$Temperature,
      estimate = pred$fit,
      lower = pred$fit - 1.96 * pred$se.fit,
      upper = pred$fit + 1.96 * pred$se.fit
   )
})

p_values_GHG <- summary_results_GHG %>%
   filter(term == "Temperaturewarm") %>%
   select(metric, p.value) %>% 
   mutate(p.value = ifelse(p.value < 0.001, 'p < 0.001', paste0('p = ', p.value)))

jpeg('figures/temp_result_GHG.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot(predictions_GHG, aes(x = Temperature, y = estimate, group = metric)) +
   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
   geom_point(pch = 23, size = 5, fill = 'coral') +
   facet_wrap(~ metric, scales = "free_y") +
   labs(
      x = "Temperature",
      y = "Predicted value"
   ) +
   geom_text(
      data = p_values_GHG,
      aes(x = 0.5, y = -Inf, label = p.value),
      hjust = 0, vjust = -1,
      size = 4, color = "black",
      inherit.aes = FALSE
   ) +
   scale_y_continuous(expand = expansion(0.2))
dev.off()

# Chlorophyll a ===================
## Data exploration =============
warm_curves_chla <- ggplot() +
   geom_point(data = warm_data,
              aes(x = DOY,
                  y = Chl_umol_L,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = warm_data,
               aes(x = DOY,
                   y = Chl_umol_L,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab(NULL) +
   ylim(0, 0.015) +
   theme(legend.position = 'none'); warm_curves_chla

cold_curves_chla <- ggplot() +
   geom_point(data = cold_data,
              aes(x = DOY,
                  y = Chl_umol_L,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = cold_data,
               aes(x = DOY,
                   y = Chl_umol_L,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab('Chlorophyll a (µmol L^-1)') +
   ylim(0, 0.015) +
   theme(legend.position = 'none'); cold_curves_chla

combined_curves_chla <- ggpubr::ggarrange(cold_curves_chla, warm_curves_chla); combined_curves_chla
ggsave(combined_curves_chla, filename = 'figures/combined_curves_chla.jpg', width = 12)

## Overall temporal trend of Chlorophyll a =============
model_chl_a <- lmer(Chl_umol_L ~ Temperature * poly(DOY, degree = 2) + (1|Mesocosm), 
                    data = data) 
# Singular fit, this could happen if there was very little variation between Mesocosms.
simulationOutput <- simulateResiduals(fittedModel = model_chl_a)
testResiduals(simulationOutput)
plot(simulationOutput)

summary_chl_a <- as.data.frame(summary(model_chl_a)[['coefficients']]) %>% 
   mutate(response = 'Chlorophyll a') %>% 
   rownames_to_column('coefficient')
anova_chl_a <- car::Anova(model_chl_a) %>% 
   mutate(response = 'Chlorophyll a') %>% 
   rownames_to_column('coefficient')

### Plot model results =============
# Fitted DOY effect
newdata_DOY_chl_a <- expand.grid(
   DOY = seq(min(data$DOY), max(data$DOY), length.out = 100),
   Temperature = c('cold', 'warm')
)

preds_chl_a_DOY <- predict(model_chl_a, newdata = newdata_DOY_chl_a, se.fit = TRUE, re.form = NA)

newdata_chl_a_DOY <- newdata_DOY_chl_a %>%
   mutate(
      fit    = (preds_chl_a_DOY$fit),
      lower  = (preds_chl_a_DOY$fit - 1.96 * preds_chl_a_DOY$se.fit),
      upper  = (preds_chl_a_DOY$fit + 1.96 * preds_chl_a_DOY$se.fit)
   )

jpeg('figures/model_result_chl_a.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot() +
   geom_point(data = data, aes(x = DOY, y = Chl_umol_L, col = Temperature), cex = 2, alpha = 0.3) +
   geom_ribbon(data = newdata_chl_a_DOY, aes(x = DOY, ymin = lower, ymax = upper, fill = Temperature), alpha = 0.3) +
   geom_line(data = newdata_chl_a_DOY, aes(x = DOY, y = fit, col = Temperature), size = 1.2) +
   scale_color_manual(values = c('blue', 'red')) +
   scale_fill_manual(values = c('blue', 'red')) +
   labs(x = "DOY",
        y = "Predicted Chlorophyll a (µmol L^-1)")
dev.off()

newdata_chl_a_DOY$response <- 'Chlorophyll a'

## Temporal trajectory analysis of sediment POC =============
### Aggregate the data =============
summary_data_chl_a <- data %>%
   group_by(Temperature, Treatment) %>%
   group_split() %>%
   map_dfr(summarize_tank, var = 'Chl_umol_L')

### Fit models =============
summary_results_chl_a <- summary_data_chl_a %>%
   select(Temperature, peak_fit, time_to_peak_fit, cum_pred_fitted) %>%   
   pivot_longer(-c(Temperature), names_to = "metric", values_to = "value") %>%
   group_by(metric) %>%
   summarize(model = list(lm(value ~ Temperature, data = cur_data())), .groups = "drop") %>%
   mutate(summary = map(model, broom::tidy)) %>%
   unnest(summary) %>% 
   mutate(p.value = round(p.value, 3))

### Plot the results =============
predictions_chl_a <- map_dfr(c("cum_pred_fitted", "peak_fit", "time_to_peak_fit"), function(m) {
   mod <- lm(as.formula(paste(m, "~ Temperature")), data = summary_data_chl_a)
   newdata <- tibble(Temperature = temp_levels)
   
   pred <- predict(mod, newdata, se.fit = TRUE)
   
   tibble(
      metric = m,
      Temperature = newdata$Temperature,
      estimate = pred$fit,
      lower = pred$fit - 1.96 * pred$se.fit,
      upper = pred$fit + 1.96 * pred$se.fit
   )
})

p_values_chl_a <- summary_results_chl_a %>%
   filter(term == "Temperaturewarm") %>%
   select(metric, p.value) %>% 
   mutate(p.value = ifelse(p.value < 0.001, 'p < 0.001', paste0('p = ', p.value)))

jpeg('figures/temp_result_chl_a.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot(predictions_chl_a, aes(x = Temperature, y = estimate, group = metric)) +
   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
   geom_point(pch = 23, size = 5, fill = 'coral') +
   facet_wrap(~ metric, scales = "free_y") +
   labs(
      x = "Temperature",
      y = "Predicted value"
   ) +
   geom_text(
      data = p_values_chl_a,
      aes(x = 0.5, y = -Inf, label = p.value),
      hjust = 0, vjust = -1,
      size = 4, color = "black",
      inherit.aes = FALSE
   ) +
   scale_y_continuous(expand = expansion(0.2))
dev.off()

# POC:Chlorophyll a Ratio ===================
## Data exploration =============
warm_curves_POC_chla <- ggplot() +
   geom_point(data = warm_data,
              aes(x = DOY,
                  y = POC_chla,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = warm_data,
               aes(x = DOY,
                   y = POC_chla,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Warm)', palette = 'Reds', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab(NULL) +
   ylim(0, 4e06) +
   theme(legend.position = 'none'); warm_curves_POC_chla

cold_curves_POC_chla <- ggplot() +
   geom_point(data = cold_data,
              aes(x = DOY,
                  y = POC_chla,
                  fill = Mesocosm_ID),
              color = 'black', pch = 21, alpha = 0.5, cex = 3) +
   geom_smooth(data = cold_data,
               aes(x = DOY,
                   y = POC_chla,
                   col = Mesocosm_ID),
               method = "lm", formula = y ~ poly(x, degree = 2), se = F) +
   scale_fill_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   scale_color_brewer(name = 'Mesocosm ID (Cold)', palette = 'Blues', direction = -1) +
   guides(size = 'none', alpha = 'none') +
   ylab('POC:Chlorophyll a Ratio') +
   ylim(0, 4e06) +
   theme(legend.position = 'none'); cold_curves_POC_chla

combined_curves_POC_chla <- ggpubr::ggarrange(cold_curves_POC_chla, warm_curves_POC_chla); combined_curves_POC_chla
ggsave(combined_curves_POC_chla, filename = 'figures/combined_curves_POC_chla.jpg', width = 12)

## Overall temporal trend of POC:Chlorophyll a Ratio =============
model_POC_chla <- lmer(POC_chla ~ Temperature * poly(DOY, degree = 2) + (1|Mesocosm), 
                       data = data) 
# Singular fit, this could happen if there was very little variation between Mesocosms.
simulationOutput <- simulateResiduals(fittedModel = model_POC_chla)
testResiduals(simulationOutput)
plot(simulationOutput)

summary_POC_chla <- as.data.frame(summary(model_POC_chla)[['coefficients']]) %>% 
   mutate(response = 'POC:Chlorophyll a Ratio') %>% 
   rownames_to_column('coefficient')
anova_POC_chla <- car::Anova(model_POC_chla) %>% 
   mutate(response = 'POC:Chlorophyll a Ratio') %>% 
   rownames_to_column('coefficient')

### Plot model results =============
# Fitted DOY effect
newdata_DOY_POC_chla <- expand.grid(
   DOY = seq(min(data$DOY), max(data$DOY), length.out = 100),
   Temperature = c('cold', 'warm')
)

preds_POC_chla_DOY <- predict(model_POC_chla, newdata = newdata_DOY_POC_chla, se.fit = TRUE, re.form = NA)

newdata_POC_chla_DOY <- newdata_DOY_POC_chla %>%
   mutate(
      fit    = (preds_POC_chla_DOY$fit),
      lower  = (preds_POC_chla_DOY$fit - 1.96 * preds_POC_chla_DOY$se.fit),
      upper  = (preds_POC_chla_DOY$fit + 1.96 * preds_POC_chla_DOY$se.fit)
   )

jpeg('figures/model_result_POC_chla.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot() +
   geom_point(data = data, aes(x = DOY, y = POC_chla, col = Temperature), cex = 2, alpha = 0.3) +
   geom_ribbon(data = newdata_POC_chla_DOY, aes(x = DOY, ymin = lower, ymax = upper, fill = Temperature), alpha = 0.3) +
   geom_line(data = newdata_POC_chla_DOY, aes(x = DOY, y = fit, col = Temperature), size = 1.2) +
   scale_color_manual(values = c('blue', 'red')) +
   scale_fill_manual(values = c('blue', 'red')) +
   labs(x = "DOY",
        y = "Predicted POC:Chlorophyll a Ratio")
dev.off()

newdata_POC_chla_DOY$response <- 'POC:Chlorophyll a Ratio'

## Temporal trajectory analysis of POC:Chlorophyll a Ratio =============
### Aggregate the data =============
summary_data_POC_chla <- data %>%
   group_by(Temperature, Treatment) %>%
   group_split() %>%
   map_dfr(summarize_tank, var = 'POC_chla')

### Fit models =============
summary_results_POC_chla <- summary_data_POC_chla %>%
   select(Temperature, peak_fit, time_to_peak_fit, cum_pred_fitted) %>%   
   pivot_longer(-c(Temperature), names_to = "metric", values_to = "value") %>%
   group_by(metric) %>%
   summarize(model = list(lm(value ~ Temperature, data = cur_data())), .groups = "drop") %>%
   mutate(summary = map(model, broom::tidy)) %>%
   unnest(summary) %>% 
   mutate(p.value = round(p.value, 3))

### Plot the results =============
predictions_POC_chla <- map_dfr(c("cum_pred_fitted", "peak_fit", "time_to_peak_fit"), function(m) {
   mod <- lm(as.formula(paste(m, "~ Temperature")), data = summary_data_POC_chla)
   newdata <- tibble(Temperature = temp_levels)
   
   pred <- predict(mod, newdata, se.fit = TRUE)
   
   tibble(
      metric = m,
      Temperature = newdata$Temperature,
      estimate = pred$fit,
      lower = pred$fit - 1.96 * pred$se.fit,
      upper = pred$fit + 1.96 * pred$se.fit
   )
})

p_values_POC_chla <- summary_results_POC_chla %>%
   filter(term == "Temperaturewarm") %>%
   select(metric, p.value) %>% 
   mutate(p.value = ifelse(p.value < 0.001, 'p < 0.001', paste0('p = ', p.value)))

jpeg('figures/temp_result_POC_chla.jpg', res = 200, width = 10, height = 7, units = 'in')
ggplot(predictions_POC_chla, aes(x = Temperature, y = estimate, group = metric)) +
   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
   geom_point(pch = 23, size = 5, fill = 'coral') +
   facet_wrap(~ metric, scales = "free_y") +
   labs(
      x = "Temperature",
      y = "Predicted value"
   ) +
   geom_text(
      data = p_values_POC_chla,
      aes(x = 0.5, y = -Inf, label = p.value),
      hjust = 0, vjust = -1,
      size = 4, color = "black",
      inherit.aes = FALSE
   ) +
   scale_y_continuous(expand = expansion(0.2))
dev.off()

# Combined plot =============
preds_combined <- rbind(newdata_DOY_POC, 
                        newdata_DOY_DOC, 
                        newdata_sedPOC_DOY,
                        newdata_GHG_DOY,
                        newdata_sedDOC_DOY,
                        newdata_chl_a_DOY,
                        newdata_POC_chla_DOY) %>% 
   mutate(response = factor(response, levels = c(
      "CO2",
      "Water Column POC", 
      "Sediment Trap POC", 
      "Water Column DOC",
      "Sediment Trap DOC",
      "Chlorophyll a",
      "POC:Chlorophyll a Ratio")))

raw_data_combined <- data %>% 
   select(DOY, Mesocosm, Temperature, POC_umol_L, Chl_umol_L, POC_chla) %>% 
   full_join(data_sediment %>% select(DOY, Mesocosm, Temperature, POC_mmol_m2_d, DOC_mmol_m2_d)) %>% 
   full_join(data_GHG %>% select(DOY, Mesocosm, Temperature, CO2_mmol_m2_d)) %>% 
   full_join(data_DOC %>% select(DOY, Mesocosm, Temperature, DOC_umol_L)) %>% 
   pivot_longer(cols = 4:NCOL(.), names_to = 'response', values_to = 'value')

raw_data_combined$response[raw_data_combined$response == 'CO2_mmol_m2_d'] <- 'CO2'
raw_data_combined$response[raw_data_combined$response == 'DOC_mmol_m2_d'] <- 'Sediment Trap DOC'
raw_data_combined$response[raw_data_combined$response == 'POC_mmol_m2_d'] <- 'Sediment Trap POC'
raw_data_combined$response[raw_data_combined$response == 'DOC_umol_L'] <- 'Water Column DOC'
raw_data_combined$response[raw_data_combined$response == 'POC_umol_L'] <- 'Water Column POC'
raw_data_combined$response[raw_data_combined$response == 'Chl_umol_L'] <- 'Chlorophyll a'
raw_data_combined$response[raw_data_combined$response == 'POC_chla'] <- 'POC:Chlorophyll a Ratio'

raw_data_combined <- raw_data_combined %>%
   mutate(facet_col = Temperature,
          response = factor(response, 
                            levels = c("CO2",
                                       "Water Column POC", 
                                       "Sediment Trap POC", 
                                       "Water Column DOC",
                                       "Sediment Trap DOC",
                                       "Chlorophyll a",
                                       "POC:Chlorophyll a Ratio")))

preds_combined <- preds_combined %>%
   mutate(facet_col = "trends")

raw_data_combined$facet_col <- factor(raw_data_combined$facet_col, levels = c("cold", "warm", "trends"))
preds_combined$facet_col <- factor(preds_combined$facet_col, levels = c("cold", "warm", "trends"))

trend_vlines <- preds_combined %>%
   group_by(Temperature, response) %>%
   summarize(vline = ifelse(response == "CO2", DOY[which.min(fit)], DOY[which.max(fit)]), .groups = "drop") %>%
   mutate(facet_col = factor("trends", levels = c("cold", "warm", "trends")))

mesocosm_colors <- c("#102952", "#10397B", "#104FA5", "#748B9A", "#99B1D8",
                     "#B4DEFF", "#732000", "#8C2E0E", "#BD4A29", "#BC4928",
                     "#BD735A", "#E6946A")
temperature_colors <- c('#104FA5', '#BD4A29') 

combined_plot <- ggplot() +
   geom_point(data = raw_data_combined %>% filter(!response %in% c("Chlorophyll a", "POC:Chlorophyll a Ratio")),
              aes(x = DOY, y = value, fill = Mesocosm),
              color = 'black', pch = 21, alpha = 0.5, size = 3, show.legend = F) +
   geom_smooth(data = raw_data_combined %>% filter(!response %in% c("Chlorophyll a", "POC:Chlorophyll a Ratio")),
               aes(x = DOY, y = value, col = Mesocosm),
               method = "lm", formula = y ~ poly(x, degree = 2), se = FALSE, show.legend = F) +
   scale_fill_manual(values = mesocosm_colors) +
   scale_color_manual(values = mesocosm_colors) +
   new_scale_color() +
   new_scale_fill() +
   geom_ribbon(data = preds_combined %>% filter(!response %in% c("Chlorophyll a", "POC:Chlorophyll a Ratio")),
               aes(x = DOY, ymin = lower, ymax = upper, fill = Temperature),
               alpha = 0.3) +
   geom_line(data = preds_combined %>% filter(!response %in% c("Chlorophyll a", "POC:Chlorophyll a Ratio")),
             aes(x = DOY, y = fit, col = Temperature),
             size = 1.2) +
   geom_vline(data = trend_vlines %>% filter(!response %in% c("Chlorophyll a", "POC:Chlorophyll a Ratio")),
              aes(xintercept = vline, col = Temperature),
              linetype = "dashed", size = 1.2, show.legend = F) +
   scale_fill_manual(name = "Temperature", values = temperature_colors, labels = c("5°C", "10°C")) +
   scale_color_manual(name = "Temperature", values = temperature_colors, labels = c("5°C", "10°C")) +
   facet_grid(response ~ facet_col, scales = 'free_y',
              labeller = as_labeller(
                 c(
                    'cold' = "Raw Data at 5°C",
                    'warm' = "Raw Data at 10°C",
                    'trends' = "Fitted Trends (Peaks Highlighted)",
                    "CO2" = "CO\u2082 [mmol m\u207B\u00B2 d\u207B\u00B9]",
                    "Water Column POC" = "POC [µmol L\u207B\u00B9]", 
                    "Sediment Trap POC" = "POC [mmol m\u207B\u00B2 d\u207B\u00B9]", 
                    "Water Column DOC" = "DOC [µmol L\u207B\u00B9]",
                    "Sediment Trap DOC" = "DOC [mmol m\u207B\u00B2 d\u207B\u00B9]"
                 )
              )) +
   scale_x_continuous(breaks = seq(min(raw_data_combined$DOY), max(raw_data_combined$DOY), 5),
                      labels = seq(min(raw_data_combined$DOY), max(raw_data_combined$DOY), 5),
                      expand = expansion(0.1)) +
   scale_y_continuous(expand = expansion(0.1)) +
   ylab(element_blank()) +
   xlab('Day of year') +
   theme(legend.position = "bottom",
         text = element_text(family = "sans")); combined_plot

combined_plot_sub <- ggplot() +
   geom_point(data = raw_data_combined %>% filter(!response %in% c("Water Column DOC", "Sediment Trap DOC")),
              aes(x = DOY, y = value, fill = Mesocosm),
              color = 'black', pch = 21, alpha = 0.5, size = 3, show.legend = F) +
   geom_smooth(data = raw_data_combined %>% filter(!response %in% c("Water Column DOC", "Sediment Trap DOC")),
               aes(x = DOY, y = value, col = Mesocosm),
               method = "lm", formula = y ~ poly(x, degree = 2), se = FALSE, show.legend = F) +
   scale_fill_manual(values = mesocosm_colors) +
   scale_color_manual(values = mesocosm_colors) +
   new_scale_color() +
   new_scale_fill() +
   geom_ribbon(data = preds_combined %>% filter(!response %in% c("Water Column DOC", "Sediment Trap DOC")),
               aes(x = DOY, ymin = lower, ymax = upper, fill = Temperature),
               alpha = 0.3) +
   geom_line(data = preds_combined %>% filter(!response %in% c("Water Column DOC", "Sediment Trap DOC")),
             aes(x = DOY, y = fit, col = Temperature),
             size = 1.2) +
   geom_vline(data = trend_vlines %>% filter(!response %in% c("Water Column DOC", "Sediment Trap DOC")),
              aes(xintercept = vline, col = Temperature),
              linetype = "dashed", size = 1.2, show.legend = F) +
   scale_fill_manual(name = "Temperature", values = temperature_colors, labels = c("5°C", "10°C")) +
   scale_color_manual(name = "Temperature", values = temperature_colors, labels = c("5°C", "10°C")) +
   facet_grid(response ~ facet_col, scales = 'free_y',
              labeller = as_labeller(
                 c(
                    'cold' = "Raw Data at 5°C",
                    'warm' = "Raw Data at 10°C",
                    'trends' = "Fitted Trends (Peaks Highlighted)",
                    "CO2" = "CO\u2082 [mmol m\u207B\u00B2 d\u207B\u00B9]",
                    "Water Column POC" = "POC [µmol L\u207B\u00B9]", 
                    "Sediment Trap POC" = "POC [mmol m\u207B\u00B2 d\u207B\u00B9]", 
                    "Chlorophyll a" = "Chl a [µmol L\u207B\u00B9]",
                    "POC:Chlorophyll a Ratio" = "POC:Chl a"
                 )
              )) +
   scale_x_continuous(breaks = seq(min(raw_data_combined$DOY), max(raw_data_combined$DOY), 5),
                      labels = seq(min(raw_data_combined$DOY), max(raw_data_combined$DOY), 5),
                      expand = expansion(0.1)) +
   scale_y_continuous(expand = expansion(0.1)) +
   ylab(element_blank()) +
   xlab('Day of year') +
   theme(legend.position = "bottom",
         text = element_text(family = "sans")); combined_plot_sub
ggsave(combined_plot_sub, filename = 'figures/combined_plot_sub.jpg', width = 20, height = 25, units = 'cm')


predictions_combined <- rbind(predictions_GHG,
                              predictions,
                              predictions_sedPOC,
                              predictions_DOC,
                              predictions_sedDOC) %>% 
   mutate(response = factor(c(
      rep("CO2", 6),
      rep("Water Column POC", 6),
      rep("Sediment Trap POC", 6),
      rep("Water Column DOC", 6),
      rep("Sediment Trap DOC", 6)),
      levels = c("CO2",
                 "Water Column POC", 
                 "Sediment Trap POC", 
                 "Water Column DOC",
                 "Sediment Trap DOC"))) 

p_values_combined <- rbind(p_values_GHG,
                           p_values %>% 
                              rename(p.value = adj_p),
                           p_values_sedPOC,
                           p_values_DOC,
                           p_values_sedDOC) %>% 
   mutate(response = factor(c(rep("CO2", 3),
                              rep("Water Column POC", 3),
                              rep("Sediment Trap POC", 3),
                              rep("Water Column DOC", 3),
                              rep("Sediment Trap DOC", 3)),
                            levels = c("CO2",
                                       "Water Column POC", 
                                       "Sediment Trap POC", 
                                       "Water Column DOC",
                                       "Sediment Trap DOC"))) 

predictions_plot <- ggplot(predictions_combined, aes(x = Temperature, y = estimate, group = metric, 
                                                     col = Temperature, fill = Temperature)) +
   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
   geom_point(pch = 23, size = 5) +
   geom_text(
      data = p_values_combined,
      aes(x = 0.5, y = -Inf, label = p.value),
      hjust = 0, vjust = -0.7,
      size = 4, color = "black",
      inherit.aes = FALSE
   ) +
   facet_grid2(response ~ metric, scales = "free", independent = "y",
              labeller = as_labeller(
                 c(
                    'cum_pred_fitted' = "Cumulative prediction",
                    'peak_fit' = "Value at Peak",
                    'time_to_peak_fit' = "Days to Peak",
                    "CO2" = "CO\u2082 [mmol m\u207B\u00B2 d\u207B\u00B9]",
                    "Water Column POC" = "POC [µmol L\u207B\u00B9]", 
                    "Sediment Trap POC" = "POC [mmol m\u207B\u00B2 d\u207B\u00B9]", 
                    "Water Column DOC" = "DOC [µmol L\u207B\u00B9]",
                    "Sediment Trap DOC" = "DOC [mmol m\u207B\u00B2 d\u207B\u00B9]"
                 )
              )) +
   labs(
      x = "Temperature",
      y = "Predicted value"
   ) +
   scale_fill_manual(name = "Temperature", values = temperature_colors, labels = c("5°C", "10°C")) +
   scale_color_manual(name = "Temperature", values = temperature_colors, labels = c("5°C", "10°C")) +
   scale_y_continuous(expand = expansion(0.2)) +
   theme(legend.position = "bottom",
         text = element_text(family = "sans")); predictions_plot

ggsave('figures/combined_plot.jpg', width = 20, height = 25, units = 'cm')

full_plot <- combined_plot | predictions_plot + plot_annotation(tag_levels = 'a')
full_plot <- ggarrange(combined_plot, predictions_plot, ncol = 2, labels = c('a', 'b'),
                       font.label = list(size = 18, color = "black", face = "bold", family = NULL))
ggsave('figures/full_plot.jpg', width = 37, height = 25, units = 'cm')

# Model results table =============
model_results <- merge(predictions_combined, p_values_combined) %>% 
   mutate(across(estimate:upper, \(x) round(x, digits = 2)),
          response = as.character(response),
          Temperature = ifelse(Temperature == 'warm', 'Temperature: warm', 'Temperature: cold')) %>% 
   select(Response = metric, `Carbon measurement` = response, Predictor = Temperature, 
          `Median estimate` = estimate, 
          `Lower confidence interval` = lower, 
          `Upper confidence interval` = upper, 
          `P-value` = p.value)

model_results$Response[model_results$Response == 'cum_pred_fitted'] = "Cumulative prediction"
model_results$Response[model_results$Response == 'peak_fit'] = "Value at Peak"
model_results$Response[model_results$Response == 'time_to_peak_fit'] = "Days to Peak"

model_results$`Carbon measurement`[model_results$`Carbon measurement` == "CO2"] = "CO\u2082 [mmol m\u207B\u00B2 d\u207B\u00B9]"
model_results$`Carbon measurement`[model_results$`Carbon measurement` == "Water Column POC"] = "POC [µmol L\u207B\u00B9]"
model_results$`Carbon measurement`[model_results$`Carbon measurement` == "Sediment Trap POC"] = "POC [mmol m\u207B\u00B2 d\u207B\u00B9]"
model_results$`Carbon measurement`[model_results$`Carbon measurement` == "Water Column DOC"] = "DOC [µmol L\u207B\u00B9]"
model_results$`Carbon measurement`[model_results$`Carbon measurement` == "Sediment Trap DOC"] = "DOC [mmol m\u207B\u00B2 d\u207B\u00B9]"


# Export =============
write.csv(predictions_combined, 'model_out/predictions_combined.csv')
write.csv(preds_combined, 'model_out/trend_predictions_combined.csv')
write_xlsx(model_results, "model_out/model_results.xlsx")
