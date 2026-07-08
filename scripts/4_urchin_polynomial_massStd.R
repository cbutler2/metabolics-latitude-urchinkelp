library(pacman)
p_load(tidyverse, ggplot2, emmeans, car, AICcmodavg, glmmTMB, DHARMa, flextable, broom.mixed)


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


# Centros ------------------------------------------------------------

### Read in data
centro.rates <- readRDS("data/centro_rates.rds")

centro.rates <- centro.rates %>% 
  filter(Contents == 'urchin',
         rate_ww > 0) %>% # this removes the 'zero' data I added for when I declared the urchins dead
  mutate(temp = as.numeric(as.character(Temp_treat))) %>% 
  droplevels()

centro.rates %>%  
  group_by(site, temp) %>% 
  summarise(N = length(rate_ww),
            mean = mean(rate_ww),
            sd = sd(rate_ww),
            se = sd/sqrt(N)) %>% 
  ggplot(aes(x = temp, y = mean, colour = site)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = palette) +
  theme_format

centro.rates.summary <- centro.rates %>%  
  group_by(site, temp) %>% 
  summarise(N = length(rate_ww),
            mean = mean(rate_ww),
            sd = sd(rate_ww),
            se = sd/sqrt(N))

head(centro.rates)

## Fit model
### use 2nd order polynomial because this fits our expectation for a unimodal TPC-like shape
fit.1 <- glmmTMB(rate_ww~poly(temp, 2)*site +(1|site:Replicate), family=Gamma(link="log"), data=centro.rates, REML=TRUE)

###Residuals look good enough
simulateResiduals(fit.1, plot=TRUE)

summary(fit.1)

###Generate type 2 SS Anova for main terms - clear significant effect of the interaction
car::Anova(fit.1)


###Plot
newdata <- centro.rates %>%
  group_by(site, Replicate) %>%
  summarise(
    temp_min = min(temp, na.rm = TRUE),
    temp_max = max(temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>% # this carries out the below code row by row from the above dataframe
  mutate(Replicate = NA) %>% 
  mutate(temp = list(seq(temp_min, temp_max, by = 0.5))) %>%
  select(site, Replicate, temp) %>%
  unnest(temp)


preds1 <- predict(fit.1, newdata = newdata, se = TRUE, type = 'response', re.form = NA)
preds1 <- data.frame(newdata, preds1)

ggplot() +
  geom_ribbon(aes(x = temp, y = fit, ymin = fit - se.fit, ymax = fit + se.fit, fill = site), alpha = 0.2, colour = NA, data = preds1) +
  geom_line(aes(x = temp, y = fit, colour = site), data = preds1) +
  geom_point(aes(x = temp, y = mean, colour = site), data = centro.rates.summary) +
  geom_errorbar(aes(x = temp, y = mean, ymin = mean - se, ymax = mean + se, colour = site), data = centro.rates.summary) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_format

preds1 <- preds1 %>% 
  mutate(site = factor(site, levels = c('Sawtell', 'Forster', 'Shellharbour', 'Merimbula', 'Mallacoota', 'Fortescue')))

### Calculate max/min of shared range
newdata %>% 
  group_by(site) %>%  
  summarise(min = min(temp), max = max(temp)) # shared range across sites is 7.62-29.6

### Calculat emedian
median(seq(7.62, 29.6, by = 0.5)) # median is 18.4

### marginal mean and pairwise comparisons
emm_sites <- emmeans(fit.1,~ site|temp, at = list(temp = c(7.62, 18.4, 29.6)))

pwpm(emm_sites) #matrix shows the EMMs along the diagonal, p-values in the upper triangle, and the differences in the lower triangle
plot(emm_sites, comparisons = TRUE) # if arrows overlap, there is no difference between the pairs


# Helios ------------------------------------------------------------

### Read in data
helio.rates <- readRDS("data/helio_rates.rds")

helio.rates <- helio.rates %>% 
  filter(Contents == 'urchin',
         rate_ww > 0) %>% # this removes the 'zero' data I added for when I declared the urchins dead
  mutate(temp = as.numeric(as.character(Temp_treat))) %>% 
  droplevels()

helio.rates <- helio.rates %>% 
  mutate(site = factor(site, levels = c('Sawtell', 'Forster', 'Shellharbour', 'Merimbula', 'Mallacoota', 'Fortescue')))


helio.rates %>%  
  group_by(site, temp) %>% 
  summarise(N = length(rate_ww),
            mean = mean(rate_ww),
            sd = sd(rate_ww),
            se = sd/sqrt(N)) %>% 
  ggplot(aes(x = temp, y = mean, colour = site)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = palette) +
  theme_format

helio.rates.summary <- helio.rates %>%  
  group_by(site, temp) %>% 
  summarise(N = length(rate_ww),
            mean = mean(rate_ww),
            sd = sd(rate_ww),
            se = sd/sqrt(N))

head(helio.rates)

## Fit model
### use 2nd order polynomial because this fits our expectation for a unimodal TPC-like shape
fit.3 <- glmmTMB(rate_ww~poly(temp, 2)*site +(1|site:Replicate), family=Gamma(link="log"), data=helio.rates, REML=TRUE)

#Residuals look good enough
simulateResiduals(fit.3, plot=TRUE)

summary(fit.3)

#Generate type 2 SS Anova for main terms - clear significant effect of the interaction
car::Anova(fit.3)


## Plot
newdata <- helio.rates %>%
  group_by(site, Replicate) %>%
  summarise(
    temp_min = min(temp, na.rm = TRUE),
    temp_max = max(temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>% # this carries out the below code row by row from the above dataframe
  mutate(Replicate = NA) %>% 
  mutate(temp = list(seq(temp_min, temp_max, by = 0.5))) %>%
  select(site, Replicate, temp) %>%
  unnest(temp)

preds2 <- predict(fit.3, newdata = newdata, se = TRUE, type = 'response', re.form = NA)
preds2 <- data.frame(newdata, preds2)

ggplot() +
  geom_ribbon(aes(x = temp, y = fit, ymin = fit - se.fit, ymax = fit + se.fit, fill = site), alpha = 0.2, colour = NA, data = preds2) +
  geom_line(aes(x = temp, y = fit, colour = site), data = preds2) +
  geom_point(aes(x = temp, y = mean, colour = site), data = helio.rates.summary) +
  geom_errorbar(aes(x = temp, y = mean, ymin = mean - se, ymax = mean + se, colour = site), data = helio.rates.summary) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_format

ggplot() +
  geom_ribbon(aes(x = temp, y = fit, ymin = fit - se.fit, ymax = fit + se.fit, fill = site), alpha = 0.2, colour = NA, data = preds2) +
  geom_line(aes(x = temp, y = fit, colour = site), data = preds2) +
  geom_point(aes(x = temp, y = mean, colour = site), data = helio.rates.summary) +
  geom_errorbar(aes(x = temp, y = mean, ymin = mean - se, ymax = mean + se, colour = site), data = helio.rates.summary) +
  facet_wrap(~site) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_format

preds2 <- preds2 %>% 
  mutate(site = factor(site, levels = c('Sawtell', 'Forster', 'Shellharbour', 'Merimbula', 'Mallacoota', 'Fortescue')))

### Calculate max/min of shared range
newdata %>% 
  group_by(site) %>%  
  summarise(min = min(temp), max = max(temp)) # shared range across sites is 6.42-29.2

### calculate median
newdata %>% 
  filter(site == 'Sawtell') %>% 
  summarise(median = median(temp)) # median is 17.9

## marginal means and pairwise comparisons
emm_sites <- emmeans(fit.3,~ site|temp, at = list(temp = c(6.42, 17.9, 29.2)))

pwpm(emm_sites) #matrix shows the EMMs along the diagonal, p-values in the upper triangle, and the differences in the lower triangle
plot(emm_sites, comparisons = TRUE) # if arrows overlap, there is no difference between the pairs

