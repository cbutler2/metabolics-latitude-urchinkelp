library(pacman)
p_load(tidyverse, ggplot2, car, glmmTMB, splines, AICcmodavg, DHARMa, metafor)


# Preliminaries -----------------------------------------------------------

##  Basic plot formatting
theme_format <- 
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.title=element_text(size=14),
        axis.text = element_text(size = 12),
        axis.ticks = element_line(colour="black"),
        title = element_text(size = 12),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.key.width=unit(1.2,"cm"),
        strip.text = element_text(size = 10, colour = "black"),
        strip.background = element_rect("white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


palette <- c("#a50f15", "#f16813", "#fde43a", "#6baed6", "#2171b5", "#08306b")


# read in data ------------------------------------------------------------


eck.rates <- readRDS("data/ecklonia_rates.rds")
eck.rates <- eck.rates %>% 
  filter(Contents == 'kelp') %>% 
  droplevels()

glimpse(eck.rates)
summary(eck.rates)


eck.rates.summary <- eck.rates %>%  
  group_by(Site, Temp_treat) %>% 
  summarise(N = length(rate_ww),
            mean = mean(rate_ww),
            sd = sd(rate_ww),
            se = sd/sqrt(N))


# Fit curves to each site ---------------------------------------------------------------

fit1 <- glmmTMB(rate_ww ~ poly(Temp_treat, 2)*Site + (1 | Site:Replicate),
                family = gaussian(),
                data = eck.rates)

res <- simulateResiduals(fit1, n = 1000)
plot(res) # these look fine

summary(fit1)

## predict
newdata <- eck.rates %>%
  group_by(Site, Replicate) %>%
  summarise(
    temp_min = min(Temp_treat, na.rm = TRUE),
    temp_max = max(Temp_treat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Replicate = NA) %>% # including this line is important because it forces predict to ignore the replicates and give population-level response
  rowwise() %>% # this carries out the below code row by row from the above dataframe
  mutate(Temp_treat = list(seq(temp_min, temp_max, by = 0.5))) %>%
  select(Site, Replicate, Temp_treat) %>%
  unnest(Temp_treat)

preds1 <- predict(fit1, newdata = newdata, se = TRUE, type = 'response', re.form = NA)

preds1 <- data.frame(newdata, preds1)


ggplot() +
  geom_ribbon(aes(x = Temp_treat, y = fit, ymin = fit - se.fit, ymax = fit + se.fit, fill = Site), alpha = 0.2, colour = NA, data = preds1) +
  geom_line(aes(x = Temp_treat, y = fit, colour = Site), data = preds1) +
  geom_point(aes(x = Temp_treat, y = mean, colour = Site), data = eck.rates.summary) +
  geom_errorbar(aes(x = Temp_treat, y = mean, ymin = mean - se, ymax = mean + se, colour = Site), data = eck.rates.summary) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_format


# Get bootstrapped parameter estimates ----------------------------------------

# extract optima

optima <- preds1 %>%
  group_by(Site) %>%
  slice_max(fit, n = 1, with_ties = FALSE) %>%
  rename(
    Topt = Temp_treat,
    Rmax = fit
  )

optima


## bootstrap optima estimates to get uncertainty

## function to extract Topt and Rmax

get_optima <- function(mod, newdata) {
  
  preds <- predict(
    mod,
    newdata = newdata,
    type = "response",
    re.form = NA   # population-level curve
  )
  
  preds_df <- data.frame(newdata, preds)
  
  preds_df %>% 
    group_by(Site) %>% 
    slice_max(preds, n = 1, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(Site, Topt = Temp_treat, Rmax = preds)
}


## run bootstrap

set.seed(123)

B <- 500

boot_res <- replicate(
  B,
  {
    y_sim <- simulate(fit1, nsim = 1)[[1]]

    fit1_sim <- update(
      fit1,
      data = transform(eck.rates, rate_ww = y_sim)
    )

    get_optima(fit1_sim, newdata)
  },
  simplify = FALSE
)

boot_df <- bind_rows(boot_res, .id = "boot")

## plot
boot_df <- boot_df %>% 
  pivot_longer(cols = Topt:Rmax, names_to = 'param', values_to = 'value')


boot.summary <- boot_df %>%  
  group_by(Site, param) %>% 
  summarise(N = length(value),
            mean = mean(value),
            sd = sd(value),
            se = sd/sqrt(N))

ggplot() +
  geom_point(aes(x = Site, y = mean), data = boot.summary) +
  geom_errorbar(aes(x = Site, y = mean, ymin = mean - se, ymax = mean + se), data = boot.summary) +
  scale_colour_manual(values = palette) +
  facet_wrap(~param, scale = 'free_y') +
  theme_format


boot_df <- boot_df %>% 
  mutate(Latitude = case_when(Site == 'Sawtell' ~ -30.375, 
                              Site == 'Forster' ~ -32.125, 
                              Site == 'Shellharbour' ~ -34.625, 
                              Site == 'Merimbula' ~ -36.875, 
                              Site == 'Mallacoota' ~ -37.625, 
                              Site == 'Fortescue' ~ -43.125)) %>% 
  select(Site, Latitude, everything())



# Meta-regression of Topt across latitude ----------------------------------------------

topt <- boot_df %>% 
  filter(param == 'Topt') %>% 
  filter(! Site == 'Sawtell') %>%
  droplevels()

## try polynomial model
fit3 <- glm(value ~ Latitude, data = topt)
fit4 <- glm(value ~ Latitude + I(Latitude^2), data = topt)

M.list=list(fit3, fit4)
aictab(M.list) #fit4 is better

simulateResiduals(fit4, plot=TRUE) # need a meta-regression type approach because variance changes between sites

topt_meta <- topt %>%
  group_by(Site, Latitude) %>%
  summarise(yi = mean(value, na.rm = TRUE),
            vi = var(value, na.rm = TRUE),
            se = sqrt(vi),
            n_boot = n(),
            .groups = "drop") # 5 obs long


rma_topt <- rma(yi = yi,
                vi = vi,
                mods = ~ Latitude + I(Latitude^2),
                data = topt_meta,
                method = "REML")


qqnorm(resid(rma_topt))
plot(rma_topt)
summary(rma_topt) # no significant effect of latitude

## predict
newdata <- data.frame(Latitude = seq(min(topt_meta$Latitude), max(topt_meta$Latitude), length.out=100))
newmods <- model.matrix(~ Latitude + I(Latitude^2), data = newdata)[, -1]  # drop intercept

preds <- predict(rma_topt, newmods = newmods)  # gives fit & CI

preds_df <- newdata %>%
  mutate(fit = preds$pred,
         lwr = preds$ci.lb,
         upr = preds$ci.ub)

topt_meta <- topt_meta %>% 
  mutate(weight = weights(rma_topt),
         cex = sqrt(weight/max(weight))*3)

ggplot(topt_meta, aes(x = Latitude, y = yi)) +
  # Bubble points weighted by precision
  geom_point(aes(size = weight), alpha = 0.8) +
  scale_size_continuous(name = "Weight") +
  # Prediction ribbon and line
  geom_ribbon(data = preds_df, aes(x = Latitude, y = fit, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.20, fill = "#4C9F70") +
  geom_line(data = preds_df, aes(x = Latitude, y = fit),
            inherit.aes = FALSE, color = "#2C7FB8", linewidth = 1.2) +
  theme_format



# Meta-regression of rmax across latitude ----------------------------------------------------------------

rmax <- boot_df %>% 
  filter(param == 'Rmax') %>% 
  filter(! Site == 'Sawtell') %>%
  droplevels()

## try polynomial model
fit6 <- glm(value ~ Latitude, data = rmax)
fit7 <- glm(value ~ Latitude + I(Latitude^2), data = rmax)

M.list=list(fit6, fit7)
aictab(M.list) # fit7 is better

simulateResiduals(fit7, plot=TRUE) # need metaregression approach so can weight observations by their variance

rmax_meta <- rmax %>%
  group_by(Site, Latitude) %>%
  summarise(yi = mean(value, na.rm = TRUE),
            vi = var(value, na.rm = TRUE),
            se = sqrt(vi),
            n_boot = n(),
            .groups = "drop")

rma_rmax <- rma(yi = yi, vi = vi, mods = ~ Latitude + I(Latitude^2), data = rmax_meta, method = "REML")


qqnorm(resid(rma_rmax))
plot(rma_rmax)
summary(rma_rmax) # significant effect of latitude

## predict
newdata <- data.frame(Latitude = seq(min(rmax_meta$Latitude), max(rmax_meta$Latitude), length.out=100))
newmods <- model.matrix(~ Latitude + I(Latitude^2), data = newdata)[, -1]  # drop intercept

preds <- predict(rma_rmax, newmods = newmods)  # gives fit & CI

preds_df <- newdata %>%
  mutate(fit = preds$pred,
         lwr = preds$ci.lb,
         upr = preds$ci.ub)

rmax_meta <- rmax_meta %>% 
  mutate(weight = weights(rma_rmax),
         cex = sqrt(weight/max(weight))*3)

ggplot(rmax_meta, aes(x = Latitude, y = yi)) +
  # Bubble points weighted by precision
  geom_point(aes(size = weight), alpha = 0.8) +
  scale_size_continuous(name = "Weight") +
  # Prediction ribbon and line
  geom_ribbon(data = preds_df, aes(x = Latitude, y = fit, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.20, fill = "#4C9F70") +
  geom_line(data = preds_df, aes(x = Latitude, y = fit),
            inherit.aes = FALSE, color = "#2C7FB8", linewidth = 1.2) +
  theme_format




# Sensitivity analysis for topt -------------------------------------------

plot(influence(rma_topt))
topt_inf <- influence(rma_topt)

# site 3 showing higher influence
levels(topt_meta$Site) # merimbula is site 3

plot(influence(rma_topt)$dfbs$Latitude)
# also site 3 maybe stands out (this plot shows the change in each coefficient when a site is left out)

diag_tab <- data.frame(
  site = 1:nrow(topt_meta),
  resid = topt_inf$inf$rstudent,
  dffits = topt_inf$inf$dffits,
  hat   = topt_inf$inf$hat,
  cookd = topt_inf$inf$cook.d,
  covr = topt_inf$inf$cov.r,
  tau2.del = topt_inf$inf$tau2.del,
  QE.del = topt_inf$inf$QE.del,
  weight = topt_inf$inf$weight,
  dfbs_intercept = topt_inf$dfbs$intrcpt,
  dfbs_latitude  = topt_inf$dfbs$Latitude  # change name if your moderator is `lat_abs` or `lat_c`
)

diag_tab <- diag_tab %>% 
  mutate(site = case_when(site == 1 ~ 'Forster', 
                          site == 2 ~ 'Shellharbour', 
                          site == 3 ~ 'Merimbula', 
                          site == 4 ~ 'Mallacoota', 
                          site == 5 ~ 'Fortescue')) %>% 
  pivot_longer(cols = resid:dfbs_latitude, names_to = 'param', values_to = 'value')

## refit without Merimbula
topt <- boot_df %>% 
  filter(param == 'Topt') %>% 
  filter(! Site == 'Sawtell' &
           ! Site == 'Merimbula') %>% 
  droplevels()

topt_meta <- topt %>%
  group_by(Site, Latitude) %>%
  summarise(yi = mean(value, na.rm = TRUE),
            vi = var(value, na.rm = TRUE),
            se = sqrt(vi),
            n_boot = n(),
            .groups = "drop") # 5 obs long


rma_topt <- rma(yi = yi, vi = vi, mods = ~ Latitude + I(Latitude^2), data = topt_meta, method = "REML")

qqnorm(resid(rma_topt))
plot(rma_topt)
summary(rma_topt) # no significant effect of latitude

newdata <- data.frame(Latitude = seq(min(topt_meta$Latitude), max(topt_meta$Latitude), length.out=100))
newmods <- model.matrix(~ Latitude + I(Latitude^2), data = newdata)[, -1]  # drop intercept
toptPreds <- predict(rma_topt, newmods = newmods)  # gives fit & CI

toptPreds_df <- newdata %>%
  mutate(fit = toptPreds$pred,
         lwr = toptPreds$ci.lb,
         upr = toptPreds$ci.ub)

topt_meta <- topt_meta %>% 
  mutate(weight = weights(rma_topt))

ggplot(aes(x = Latitude, y = yi), data = topt_meta) +
  geom_point(aes(size = weight), shape = 21, data = topt_meta) +
  geom_ribbon(aes(x = Latitude, y = fit, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.10, data = toptPreds_df) +
  geom_line(aes(x = Latitude, y = fit),
            inherit.aes = FALSE, data = toptPreds_df, ) +
  scale_size_continuous(name = bquote("Weight (1/SE"^2*")")) +
  theme_format +
  theme(plot.margin = margin(0.1, 1, 1, 1, "cm")) +
  labs(title = 'b) Topt', x = 'Latitude', y = expression(paste("Temperature (",degree,"C)")))


# Sensitivity for rmax ----------------------------------------------------

plot(influence(rma_rmax))
rmax_inf <- influence(rma_rmax)

# site 3 showing higher influence
levels(rmax_meta$Site) # merimbula is site 3

plot(influence(rma_rmax)$dfbs$Latitude)
# also site 3 maybe stands out (this plot shows the change in each coefficient when a site is left out)

diag_tab <- data.frame(
  site = 1:nrow(rmax_meta),
  resid = rmax_inf$inf$rstudent,
  dffits = rmax_inf$inf$dffits,
  hat   = rmax_inf$inf$hat,
  cookd = rmax_inf$inf$cook.d,
  covr = rmax_inf$inf$cov.r,
  tau2.del = rmax_inf$inf$tau2.del,
  QE.del = rmax_inf$inf$QE.del,
  weight = rmax_inf$inf$weight,
  dfbs_intercept = rmax_inf$dfbs$intrcpt,
  dfbs_latitude  = rmax_inf$dfbs$Latitude  # change name if your moderator is `lat_abs` or `lat_c`
)

diag_tab <- diag_tab %>% 
  mutate(site = case_when(site == 1 ~ 'Forster', 
                          site == 2 ~ 'Shellharbour', 
                          site == 3 ~ 'Merimbula', 
                          site == 4 ~ 'Mallacoota', 
                          site == 5 ~ 'Fortescue')) %>% 
  pivot_longer(cols = resid:dfbs_latitude, names_to = 'param', values_to = 'value')

## refit without Merimbula
rmax <- boot_df %>% 
  filter(param == 'Rmax') %>% 
  filter(! Site == 'Sawtell' &
           ! Site == 'Merimbula') %>% 
  droplevels()

rmax_meta <- rmax %>%
  group_by(Site, Latitude) %>%
  summarise(yi = mean(value, na.rm = TRUE),
            vi = var(value, na.rm = TRUE),
            se = sqrt(vi),
            n_boot = n(),
            .groups = "drop") # 5 obs long


rma_rmax <- rma(yi = yi, vi = vi, mods = ~ Latitude + I(Latitude^2), data = rmax_meta, method = "REML")

qqnorm(resid(rma_rmax))
plot(rma_rmax)
summary(rma_rmax) # no significant effect of latitude

newdata <- data.frame(Latitude = seq(min(rmax_meta$Latitude), max(rmax_meta$Latitude), length.out=100))
newmods <- model.matrix(~ Latitude + I(Latitude^2), data = newdata)[, -1]  # drop intercept
rmaxPreds <- predict(rma_rmax, newmods = newmods)  # gives fit & CI

rmaxPreds_df <- newdata %>%
  mutate(fit = rmaxPreds$pred,
         lwr = rmaxPreds$ci.lb,
         upr = rmaxPreds$ci.ub)

rmax_meta <- rmax_meta %>% 
  mutate(weight = weights(rma_rmax))

ggplot(aes(x = Latitude, y = yi), data = rmax_meta) +
  geom_point(aes(size = weight), shape = 21, data = rmax_meta) +
  geom_ribbon(aes(x = Latitude, y = fit, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.10, data = rmaxPreds_df) +
  geom_line(aes(x = Latitude, y = fit),
            inherit.aes = FALSE, data = rmaxPreds_df, ) +
  scale_size_continuous(name = bquote("Weight (1/SE"^2*")")) +
  theme_format +
  theme(plot.margin = margin(0.1, 1, 1, 1, "cm")) +
  labs(title = 'b) Rmax', x = 'Latitude', y = expression(paste("Temperature (",degree,"C)")))


