###############################################################################

# authors: Jenn Lavers, Lise Fournier-Carnoy, Alex Bond
# project: Sable Shearwater, adult CMR survival
# data: LHI SBSH 2010-2025 CMR data

# script objective: clean up raw data for analysis.

###############################################################################

rm(list = ls())

library(tidyverse) # for data manipulation
library(stringr) # for fixing capture histories
library(RMark) # for running rmark
library(zoo)

## Read in raw CMR data -------------------------------------------------------

dat <- read.csv("data/raw/LHI_FFSH_clean_live_recap.csv", skip = 1) %>% 
  mutate(across(where(is.integer), ~ ifelse(is.na(.), 0, 1))) %>% # fill NAs with 0 and non-NAs with 1
  rename_with(~ gsub("^X", "", .), starts_with("X")) %>% 
  mutate(`2019` = ifelse(is.na(`2019`), 0, `2019`)) %>% # no recaptures for 2019 breeding season, fill with 0
  glimpse()

# PLOT: number of birds marked, and number of bird recaptured
dat %>%
  select(`2009`:`2024`) %>%
  # for each bird, find the first year they were marked (first column with a 1)
  mutate(first_year = names(.)[max.col(. != 0, ties.method = "first")]) %>%
  count(first_year) %>%
  mutate(first_year = as.integer(first_year)) %>%
  ggplot(aes(x = first_year, y = n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  scale_x_continuous(breaks = 2009:2024) +
  labs(
    y = "Number of birds marked"
  ) +
  theme_classic()

# remove birds from early years. (there are few birds in these years, so estimating survival in these years is difficult numerically and not very informative)
dat <- dplyr::select(dat, -c(`2009`:`2012`))

# add a 'capture history' (ch) column
dat <- dat %>%
  mutate(ch = apply(select(., first(starts_with("20")):`2024`), 1, paste0, collapse = "")) %>%
  glimpse()

# format chick data so that you're culling chick years
dat <- dat %>%
  mutate(ch = ifelse(Age == "P", 
                     sub("^.*?1.*?(1.*)$", "\\1", ch), 
                     ch)) %>%
  glimpse()

# fill the culled capture histories so everything's the right format
max_length <- max(nchar(dat$ch))  # Find the longest string length

dat <- dat %>%
  mutate(ch = str_pad(ch, width = max_length, side = "left", pad = "0")) %>% # pad the character string with zeroes on the left
  glimpse()

# age class '?' are '1+' birds, so technically adults
dat$Age <- ifelse(dat$Age == "1+", "A", dat$Age)

# remove transients (birds banded as chicks that were only seen once- adults who are only seen once are not considered as transients)
dat <- dat[dat$Age != "P" | str_count(dat$ch, "1") > 1, ]
dat <- dat[str_count(dat$ch, "1") != 0, ]


glimpse(dat)
colSums(dat[3:14]) # make sure the first few years actually have birds in them, otherwise remove those years (line 26)

# save the cleaned data
saveRDS(dat, "data/tidy/LHI_FFSH_capture_histories.rds")


# PLOT 
year_data <- dat %>%
  select(`2013`:`2024`) %>%
  mutate(first_year = names(.)[max.col(. != 0, ties.method = "first")]) %>%
  pivot_longer(`2013`:`2024`, names_to = "Year", values_to = "captured")

marked <- year_data %>%
  filter(captured == 1, Year == first_year) %>%
  count(Year) %>%
  mutate(type = "Banded")

recaptured <- year_data %>%
  filter(captured == 1, Year != first_year) %>%
  count(Year) %>%
  mutate(type = "Recaptured")

recap_summary <- dat %>%
  mutate(total_captures = rowSums(across(`2013`:`2024`)),
         times_recaptured = total_captures - 1) %>%
  summarise(
    total_banded = n(),
    total_recaptures = sum(times_recaptured),
    recaptured_once = sum(times_recaptured == 1),
    recaptured_more_than_once = sum(times_recaptured > 1)
  )
summary_text <- paste0(
  "Total banded: ", recap_summary$total_banded, "\n",
  "Total recaptures: ", recap_summary$total_recaptures, "\n",
  "Recaptured once: ", recap_summary$recaptured_once, "\n",
  "Recaptured >1 time: ", recap_summary$recaptured_more_than_once
)

p <- bind_rows(marked, recaptured) %>%
  mutate(Year = as.integer(Year)) %>%
  complete(Year = 2013:2024, type = c("Banded", "Recaptured"), fill = list(n = 0)) %>%
  ggplot(aes(x = Year, y = n, fill = type)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7), vjust = -0.5, size = 3) +
  scale_x_continuous(breaks = 2013:2024) +
  ylim(c(0, 40)) +
  scale_fill_manual(values = c("Banded" = "darkgray", "Recaptured" = "steelblue")) +
  labs(y = "Number of birds", fill = NULL) +
  annotate("text", x = Inf, y = Inf, label = summary_text,
           hjust = 1.05, vjust = 1.2, size = 3.5, lineheight = 1.1) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("plots/summary_captures_recaptures.png", plot = p, width = 8, height = 5, dpi = 500)


## Environmental data ---------------------------------------------------------

# -- we want to understand whether/what environmental index has an influence over survival.
# -- the wintering period (May-September) is what we will test. We'll use a 3-month average of each variable (July-September)

### El Nino Southern Oscillation (ENSO) ---------------------------------------

enso <- read.csv("data/raw/ENSO values.csv", skip = 1) %>%  # skip the first row (data source info)
  rename(# rename for clarity
         ENSO_jun_jul_aug = JJA,
         #ENSO_jul_aug_sep = JAS,
         #ENSO_aug_sep_oct = ASO
         ) %>% 
  glimpse()

# Select the desired columns
enso <- enso[, c("time", 
                 "ENSO_jun_jul_aug"
                 #"ENSO_jul_aug_sep" 
                 #"ENSO_aug_sep_oct"
                 )]
enso <- enso[enso$time >= 2013 & enso$time <= 2024,]

print(enso) # check the format - years have to be called 'time' for rmark to recognise it.

# save the tidy environmental dat
saveRDS(enso, "data/tidy/ENSO_tidy.rds")


### Arctic Oscillation --------------------------------------------------------

AO <- read.csv("data/raw/AO values.csv", skip = 1) %>%  # skip the first row (data source info)
  glimpse()

# Select the desired columns
ao_clean <- AO %>%
  mutate(
    AO_jun_jul_aug = rowMeans(select(., June, July, August), na.rm = TRUE)
    #AO_jul_aug_sep = rowMeans(select(., July, August, September), na.rm = TRUE)
  ) %>%
  dplyr::select(time, 
                AO_jun_jul_aug
                #AO_jul_aug_sep
  ) %>% 
  dplyr::filter(time >= 2013 & time <= 2024) %>% 
  glimpse()

# save the tidy environmental dat
saveRDS(ao_clean, "data/tidy/AO_tidy.rds")


### Wintering temperature -----------------------------------------------------

sst_w <- read.csv("data/raw/SST japan values.csv", skip = 1) %>%
  glimpse()

sst_w <- data.frame(time = sst_w[,1], sst_w_mean = rowMeans(sst_w[,-1]))

plot(sst_w$sst_w_mean ~ sst_w$time)

saveRDS(sst_w, "data/tidy/temp_w_tidy.rds")



### END ###