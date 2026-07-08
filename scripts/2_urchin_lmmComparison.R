library(pacman)
p_load(tidyverse, ggplot2, car, lmerTest, flextable, broom.mixed)


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


# LMM for difference in resp between species ------------------------------------------

## load helio data
helio.rates <- readRDS("data/helio_rates.rds")

helio.rates <- helio.rates %>% 
  filter(Contents == 'urchin',
         rate_ww > 0) %>%
  mutate(temp = as.numeric(as.character(Temp_treat))) %>% 
  droplevels()


### add latitude
helio.rates <- helio.rates %>% 
  mutate(latitude = case_when(site == 'Sawtell' ~ -30.375, 
                              site == 'Forster' ~ -32.125, 
                              site == 'Shellharbour' ~ -34.625, 
                              site == 'Merimbula' ~ -36.875, 
                              site == 'Mallacoota' ~ -37.625, 
                              site == 'Fortescue' ~ -43.125))

## load centro data
centro.rates <- readRDS("data/centro_rates.rds")

centro.rates <- centro.rates %>% 
  filter(Contents == 'urchin',
         rate_ww > 0) %>%
  mutate(temp = as.numeric(as.character(Temp_treat))) %>% 
  droplevels()

### add latitude
centro.rates <- centro.rates %>% 
  mutate(latitude = case_when(site == 'Sawtell' ~ -30.375, 
                              site == 'Forster' ~ -32.125, 
                              site == 'Shellharbour' ~ -34.625, 
                              site == 'Merimbula' ~ -36.875, 
                              site == 'Mallacoota' ~ -37.625, 
                              site == 'Fortescue' ~ -43.125))


## fit linear mixed model
d1 <- rbind(centro.rates, helio.rates)

fit <- lmer(log(rate_abs) ~ species + (1|latitude:Replicate) + (1|temp), data = d1)

Anova(fit, type = 3, test = "F") 

summary(fit)

opar <- par(mfrow=c(2,2))

plot(fitted(fit), resid(fit),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "black")

qqnorm(resid(fit))
qqline(resid(fit))

plot(fitted(fit), sqrt(abs(resid(fit))),
     xlab = "Fitted Values", ylab = "Sqrt |Residuals|",
     main = "Scale-Location Plot")
abline(h = 0, col = "black")

par(opar)

as_flextable(fit)

### calculate how much less o2 helio consumes relative to centro
(exp(-0.25) - 1) *100 # therefore consumes 22.1% less O2

### calculate mean consumption for each species
d1 %>% 
  group_by(species) %>% 
  summarise(N = length(rate_abs),
            mean = mean(rate_abs),
                 sd = sd(rate_abs),
                 se = sd/sqrt(N))
