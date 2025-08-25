check_residuals <- function(residuals){
      par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
      qqnorm(residuals, pch = 1, frame = FALSE)
      qqline(residuals, col = "steelblue", lwd = 2)
      hist(residuals)
      mtext(paste0('Shapiro-Wilk test: ', round(shapiro.test(residuals)[[2]], digits = 3)), 
            side = 3, outer = TRUE, line = 0, cex = 1.5) 
      par(mfrow = c(1, 1))
}

get_cum_pred <- function(mod, data, x, baseline = 0) {
   integrate(function(t) predict(mod, setNames(data.frame(t), x)) - baseline,
             min(data[[x]]), max(data[[x]]))$value
}

summarize_tank <- function(df, var) {
   poly_fit_raw <- lm(get(var) ~ poly(DOY, 2, raw = TRUE), data = df)
   poly_fit <- lm(get(var) ~ poly(DOY, 2), data = df)
   fitted_cum_pred <- get_cum_pred(poly_fit_raw, data = df, x = 'DOY')
   peak_raw <- ifelse(grepl('CO2', var), min(df[, var], na.rm = TRUE), max(df[, var], na.rm = TRUE))
   peak_fit <- ifelse(grepl('CO2', var), min(poly_fit_raw$fitted.values), max(poly_fit_raw$fitted.values)) 
   
   time_to_peak_raw <- df$DOY[which.max(df[[var]])] - min(df$DOY)
   time_to_peak_fit <- ifelse(grepl('CO2', var), df$DOY[which.min(poly_fit_raw$fitted.values)] - min(df$DOY), 
                              df$DOY[which.max(poly_fit_raw$fitted.values)] - min(df$DOY))

   if(NROW(df[df$DOY < df$DOY[df[, var] == max(df[, var])], ]) > 2 & NROW(df[df$DOY > df$DOY[df[, var] == max(df[, var])], ]) > 2){
      n <- nrow(df)
      early_slope <- coef(lm(get(var) ~ DOY, data = df[df$DOY < df$DOY[df[, var] == max(df[, var])], ]))[[2]]
      late_slope <- coef(lm(get(var) ~ DOY, data = df[df$DOY > df$DOY[df[, var] == max(df[, var])], ]))[[2]]
   } else {
      early_slope <- NA
      late_slope <- NA
   }
   
   tibble(
      Temperature = unique(df$Temperature),
      Treatment = unique(df$Treatment),
      poly_intercept = coef(poly_fit)[[1]],
      poly_linear = coef(poly_fit)[[2]],
      poly_quadratic = ifelse(length(coef(poly_fit)) > 2, coef(poly_fit)[[3]], NA),
      cum_pred_fitted = fitted_cum_pred,
      peak_raw = peak_raw,
      time_to_peak_raw = time_to_peak_raw,
      peak_fit = peak_fit,
      time_to_peak_fit = time_to_peak_fit,
      mean = mean(df[[var]], na.rm = TRUE),
      median = median(df[[var]], na.rm = TRUE),
      early_slope = early_slope,
      late_slope = late_slope
   )
}

plot_correlation <- function(data, vars) {
   corrplot::corrplot(cor(data %>% select(all_of(vars))), 
                      method = "circle", type = "upper", order = "hclust",
                      addCoef.col = "black", number.cex = 0.65, 
                      tl.col = "black", tl.srt = 45, diag = FALSE)
}