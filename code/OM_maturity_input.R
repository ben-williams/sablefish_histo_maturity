# sablefish histological data only
# ben.williams@noaa.gov
# 2020-05

# Notes: dropped early legs from ak summer survey
# see Rodgveller 2018 (poor classification early in season)

# load ----
source(here::here("code/helper.r"))
theme_set(theme_present())

# data ---- 

read_csv("data/WCGBT_ODFW_WDFW_2010_2018_sablefishmaturity_updated022720.csv") %>% 
  dplyr::select(year = Year, date = Date_collected, lat = latitude_dd,
                long = longitude_dd, depth = depth_m, length = Length_cm, 
                age = Age_yrs, weight = Weight_kg, mature = Mature) %>% 
  mutate(Year = factor(year),
         Age = factor(age),
         Mature = factor(mature),
         date = mdy(date),
         day = yday(date)) %>% 
  arrange(year) -> wc 

wc %>% 
  mutate(age = ifelse(age>20, 20, age),
         Age = factor(age), 
         Location = 'WC',
         Breaks = factor(case_when(lat<=36 ~ "R1",
                                   lat>36 & lat <=50 ~ "R2",
                                   lat>50 & long>= -145 ~ "R3",
                                   long< -145 ~ "R4"))) -> wc_plus

read_csv("data/AFSC_winter_2011_histo_v2.csv") %>% 
  dplyr::select(depth = depth_m, 
                weight = TOTALWT,
                length = TOTLENGTH,
                age = AGE,
                lat,
                long,
                mature = Skip_mat) %>% 
  mutate(weight = weight / 1000, 
         length = length / 10, 
         Age = factor(age),
         Mature = factor(mature), 
         year = 2011, 
         Season = 'winter',
         Year = factor(year),
         Area = NA,
         date = NA,
         Atresia = NA,
         julian = NA,
         Location = 'AK') -> ak_2011

read_csv("data/AFSC_winter_2015_histo_v2.csv") %>% 
  dplyr::select(date, 
                mature = Skip_mature,
                length,
                weight = weight_g,
                age = AGE,
                lat = st_lat,
                long = end_long,
                sdepth = st_dep_m,
                edepth = end_dep_m) %>% 
  mutate(date = mdy(date),
         year = year(date),
         month = month(date),
         julian = yday(date),
         weight = weight / 1000, 
         length = length / 10, 
         depth = (sdepth + edepth) / 2,
         Age = factor(age),
         Mature = factor(mature), 
         Season = 'winter',
         Year = factor(year),
         Area = NA,
         Atresia = NA,
         lat = as.numeric(str_sub(lat, 1,2)) + 
           as.numeric(str_sub(lat, 4,5))/60 + 
           as.numeric(str_sub(lat, 7,8))/3600,
         long = -as.numeric(str_sub(long, 1,3)) + 
           as.numeric(str_sub(long, 5,6))/60 + 
           as.numeric(str_sub(long, 8,9))/3600,
         Location = 'AK') %>% 
  dplyr::select(-sdepth, -edepth) -> ak_2015

read_csv("data/AFSC_summer_2015_macro_histo_v2.csv") %>% 
  filter(!is.na(date, leg>5)) %>% 
  dplyr::select(date, 
                mature = histo_mat,
                length,
                weight = total_wt,
                age = Age,
                long = Long,
                lat = Lat,
                code = Depth_stratum_code, 
                macro = `at-sea_mat`) %>% 
  mutate(mature = ifelse(mature>=1, 1, 0),
         date = mdy(date),
         year = year(date),
         month = month(date),
         julian = yday(date),
         weight = weight / 1000, 
         length = length / 10, 
         depth = case_when(code == 1 ~ 50,
                           code == 2 ~ 150,
                           code == 3 ~ 250,
                           code == 4 ~ 350,
                           code == 5 ~ 450,
                           code == 6 ~ 500,
                           code == 7 ~ 700,
                           code == 8 ~ 900),
         Age = factor(age),
         Mature = factor(mature), 
         Season = 'summer',
         Year = factor(year),
         Area = NA,
         Atresia = NA,
         Location = 'AK') %>% 
  dplyr::select(-code) -> ak_2015s

bind_rows(ak_2011, ak_2015, ak_2015s) %>% 
  mutate(day = yday(date),
         Location = 'AK') -> ak

ak %>% 
  mutate(age = ifelse(age>20, 20, age),
         Age  = factor(age),
         Breaks = factor(case_when(lat<=36 ~ "R1",
                                   lat>36 & lat <=50 ~ "R2",
                                   lat>50 & long>= -145 ~ "R3",
                                   long< -145 ~ "R4"))) %>%
  filter(lat<60) -> ak_plus


bind_rows(wc_plus, ak_plus) %>% 
  mutate(Breaks = factor(Breaks)) %>% 
  filter(age>0, !is.na(Breaks)) -> dat

# add relative body condition Kr

lw <- lm(log(weight) ~ log(length) * Age, data = dat)

# removed the large outliers
dat %>% 
  mutate(w = exp(predict(lw, .)) * exp((sigma(lw)^2) / 2),
         Kr = weight / w) %>% 
  filter(Kr<2) -> dat

# eda ----

dat %>% 
  ggplot(aes(depth, fill = Breaks)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~Breaks)

# model ----

fit <- gam(mature ~ s(length, k=4, by = Age) + te(long, lat, depth) + 
             s(age, k=3) + s(Kr, k=3) + Year, 
           data = dat, gamma = 1.4, family = "binomial")

fit1 <- gam(mature ~ s(length, k=4, by = Age) + te(long, lat, depth) + 
             s(age, k=3) + s(Kr, k=3), 
           data = dat, gamma = 1.4, family = "binomial")

fit2 <- gam(mature ~ s(length, k=4, by = Age) + te(long, lat, depth) + 
              s(age, k=3) , 
            data = dat, gamma = 1.4, family = "binomial")

fit3 <- gam(mature ~ s(length, k=4, by = Age) + te(long, lat, depth) + 
             s(Kr, k=3), 
           data = dat, gamma = 1.4, family = "binomial")

fit4 <- gam(mature ~ s(length, k=4, by = Age) + 
              s(age, k=3) + s(Kr, k=3), 
           data = dat, gamma = 1.4, family = "binomial")

fit5 <- gam(mature ~ te(long, lat, depth) + 
             s(age, k=3) + s(Kr, k=3) + Year, 
           data = dat, gamma = 1.4, family = "binomial")

fit6 <- gam(mature ~ s(length, k=4, by = Age) + s(depth,k=3) + 
              s(age, k=3) + s(Kr, k=3), 
            data = dat, gamma = 1.4, family = "binomial")

fit7 <- gam(mature ~ s(length, k=4, by = Age) + te(long, lat) + 
              s(age, k=3) + s(Kr, k=3), 
            data = dat, gamma = 1.4, family = "binomial")

fit8 <- gam(mature ~ s(length, k=4, by = Age) + s(depth,k=3) + Breaks + 
              s(age, k=3) + s(Kr, k=3), 
            data = dat, gamma = 1.4, family = "binomial")

fit9 <- gam(mature ~ s(length, k=4, by = Age) + Breaks + s(age, k=3) + 
              s(Kr, k=3) + te(long, lat), 
            data = dat, gamma = 1.4, family = "binomial")

AIC(fit, fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9) %>% 
  rownames_to_column() %>% 
  mutate(delta = min(AIC) - AIC) %>% 
  arrange(-delta)

summary(fit1)

# check fit 

m <- getViz(fit1)

print(plot(m, allTerms = T), pages = 1) # Calls print.plotGam()

dat %>% 
  group_by(Breaks, age) %>% 
  summarise(long = median(long),
            lat = median(lat),
            length = median(length)) -> lngs

dat %>% 
  group_by(Breaks) %>% 
  summarise(depth = median(depth)) -> dpths

dat %>% 
  group_by(age, Breaks) %>% 
  summarise(long = median(long),
            lat = median(lat),
            length = median(length),
            depth = median(depth, na.rm = T),
            Kr = median(Kr, na.rm = T)) %>% 
  left_join(expand.grid(year = unique(dat$year), 
                        Breaks = factor(c("R1", "R2", "R3", "R4")),
                       age = 1:20), .) %>% 
  group_by(age) %>% 
  mutate(length = ifelse(is.na(length), median(length, na.rm=T), length),
         long = ifelse(is.na(long), median(long, na.rm=T), long),
         lat = ifelse(is.na(lat), median(lat, na.rm=T), lat),
         Year = factor(year),
         Age = factor(age)) %>% 
  ungroup() %>% 
  mutate(fit = predict.gam(fit1, ., type = "response")) %>% 
  ggplot(aes(age, fit, color = Breaks)) + 
  geom_point() +
  geom_line() +
  scale_color_viridis_d() +
  geom_hline(yintercept = c(0.25, 0.5, 0.75, 0.95), lty = 3, alpha = 0.5) +
  ggtitle(bquote("A"[50]== "4.6, 6.9, 8.2")) +
  # coord_cartesian(xlim = c(4,10.2), ylim = c(0.90, 0.96)) +
  scale_y_continuous(breaks =seq(0, 1)) 
  # geom_vline(xintercept = c(10.4))

# perc: 0.25, 0.5, 0.75, 0.95  
# r1:   7.55, 8.3,  8.7, 9.5
# r2:      4, 4.5,  5.1, 6.5
# r3:    4.2, 4.6, 4.97, 7.7
# r4:   6.05, 6.8,  7.7, 10.4


