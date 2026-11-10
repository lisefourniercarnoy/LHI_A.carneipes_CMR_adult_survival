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

# remove transients (birds banded as chicks that were only seen once- adults who are only seen once are not considered as transients)
chicks <- dat |> dplyr::filter(Age == "P" & str_count(ch, "1") > 1)
adults <- dat |> dplyr::filter(Age != "P")
dat <- rbind(chicks, adults)

# format chick data so that you're culling chick years
dat <- dat %>%
  mutate(ch = ifelse(Age == "P", 
                     sub("^.*?1.*?(1.*)$", "\\1", ch), 
                     ch)) %>%
  glimpse()

# zero the chick-year capture
cols   <- as.character(2013:2024)
chick_rows <- which(dat$Age == "P")

for (ROW in chick_rows) {
  COL <- 1
  while (COL <= length(cols) && dat[ROW, cols[COL]] != 1) {
    COL <- COL + 1
  }
  if (COL <= length(cols)) {
    dat[ROW, cols[COL]] <- 0
  }
}

# fill the culled capture histories so everything's the right format
max_length <- max(nchar(dat$ch))  # Find the longest string length

dat <- dat %>%
  mutate(ch = str_pad(ch, width = max_length, side = "left", pad = "0")) %>% # pad the character string with zeroes on the left
  glimpse()

# age class '?' are '1+' birds, so technically adults
dat$Age <- ifelse(dat$Age == "1+", "A", dat$Age)


glimpse(dat)
colSums(dat[3:14]) # make sure the first few years actually have birds in them, otherwise remove those years (line 26)

# save the cleaned data
saveRDS(dat, "data/tidy/LHI_FFSH_capture_histories.rds")


### Figure 1. summary of mark-recapture sample --------------------------------

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
  ylim(c(0, max(c(marked$n, recaptured$n))+7)) +
  scale_fill_manual(values = c("Banded" = "darkgray", "Recaptured" = "steelblue")) +
  labs(y = "Number of birds", fill = NULL) +
  annotate("text", x = Inf, y = Inf, label = summary_text,
           hjust = 1.05, vjust = 1.2, size = 3.5, lineheight = 1.1) +
  theme_classic() +
  theme(legend.position = "bottom")
p
ggsave("plots/fig_1_summary_captures_recaptures.png", plot = p, width = 8, height = 5, dpi = 500)


## Environmental data ---------------------------------------------------------

# -- we want to understand whether/what environmental index has an influence over survival.
# -- the wintering period (May-September) is what we will test. We'll use a 3-month average of each variable (June-August)

### El Nino Southern Oscillation (ENSO) ---------------------------------------

enso <- read.csv("data/raw/ENSO values.csv", skip = 1) %>%  # skip the first row (data source info)
  rename(# rename for clarity
         ENSO_jun_jul_aug = JJA,
         ENSO_oct_nov_dec = OND,
         ENSO_feb_mar_apr = FMA
         ) %>% 
  glimpse()

# Select the desired columns
enso <- enso[, c("time", 
                 "ENSO_jun_jul_aug",
                 "ENSO_oct_nov_dec" ,
                 "ENSO_feb_mar_apr"
                 )]
enso <- enso[enso$time >= 2013 & enso$time <= 2024,]

print(enso) # check the format - years have to be called 'time' for rmark to recognise it.

# save the tidy environmental dat
saveRDS(enso, "data/tidy/ENSO_tidy.rds")


### Arctic Oscillation (AO) ---------------------------------------------------

AO <- read.csv("data/raw/AO values.csv", skip = 1) %>%  # skip the first row (data source info)
  glimpse()

# Select the desired columns
ao_clean <- AO %>%
  mutate(
    AO_jun_jul_aug = rowMeans(select(., June, July, August), na.rm = TRUE),
    AO_oct_nov_dec = rowMeans(select(., October, November, December), na.rm = TRUE),
    AO_feb_mar_apr = rowMeans(select(., February, March, April), na.rm = TRUE),
    
  ) %>%
  dplyr::select(time, 
                AO_jun_jul_aug,
                AO_oct_nov_dec,
                AO_feb_mar_apr
  ) %>% 
  dplyr::filter(time >= 2013 & time <= 2024) %>% 
  glimpse()

# save the tidy environmental dat
saveRDS(ao_clean, "data/tidy/AO_tidy.rds")

### Sea surface temperature (SST) ---------------------------------------------

#### Wintering SST ------------------------------------------------------------

sst_winter <- read.csv("data/raw/SST japan values.csv", skip = 1) %>%
  glimpse()

sst_winter <- data.frame(time = sst_winter[,1], sst_winter = rowMeans(sst_winter[,-1])) |> dplyr::filter(time >= 2013)

par(mfrow = c(3,1)); plot(sst_winter$sst_winter ~ sst_winter$time)


#### Breeding SST -------------------------------------------------------------

## sst is the only environmental thing that is missing data to make an 'in-between breeding seasons' metric (july-sept.) so we are lumping june-oct for this one only.

sst <- read.csv("data/raw/IMOS_-_Australian_National_Mooring_Network_(ANMN)_-_CTD_Profiles.csv", skip = 29) %>%
  mutate(datetime = ymd_hms(TIME)) %>% 
  dplyr::select(c(datetime, 
                  LATITUDE, LONGITUDE, 
                  DEPTH,
                  TEMP, TEMP_quality_control)) %>% # select relevant columns
  dplyr::filter(DEPTH > 10 & DEPTH < 15) %>% 
  dplyr::filter(TEMP_quality_control == 1 | TEMP_quality_control == 2) %>%  # quality control codes 1 = good data, 2 = prbably good data, 4 is bad data
  group_by(datetime) %>% 
  mutate(temp_10_15m = mean(TEMP)) %>% 
  glimpse()
unique(sst$datetime)


# Extract year and month
sst_breeding <- sst %>%
  mutate(
    time = year(datetime),
    month = month(datetime)
  ) %>%
  group_by(time) %>%
  summarise(

    sst_oct_nov_dec = mean(temp_10_15m[month %in% c(10:12)], na.rm = TRUE),
    sst_feb_mar_apr = mean(temp_10_15m[month %in% c(2:4)], na.rm = TRUE),

    .groups = "drop"
  )

sst_breeding <- sst_breeding %>% dplyr::filter(time >= 2013, 
                                               time < 2025)
print(sst_final)

plot(sst_breeding$sst_oct_nov_dec ~ sst_breeding$time)
plot(sst_breeding$sst_feb_mar_apr ~ sst_breeding$time)


# merge and save
sst_final <- merge(sst_winter, sst_breeding)
saveRDS(sst_final, "data/tidy/temp_tidy.rds")


### Pacific Decadal Oscillation (PDO) -----------------------------------------

pdo <- readLines("data/raw/PDO values.txt")
pdo_data <- pdo[-c(1,2)] # Remove the first two header lines

pdo_split <- strsplit(pdo_data, "\\s+") # Split each line by whitespace
pdo_matrix <- do.call(rbind, lapply(pdo_split, function(x) x[nzchar(x)]))
pdo <- as.data.frame(pdo_matrix, stringsAsFactors = FALSE) # Convert to data frame 
colnames(pdo) <- c("Year", month.abb) # name columns
pdo[] <- lapply(pdo, function(x) ifelse(suppressWarnings(!is.na(as.numeric(x))), as.numeric(x), x)) # Convert numeric columns

pdo_clean <- pdo %>%
  mutate(
    time = Year,
    PDO_jul_aug_sep = rowMeans(select(., Jul, Aug, Sep), na.rm = TRUE)
  ) %>%
  dplyr::select(time, PDO_jul_aug_sep) %>% 
  dplyr::filter(time >= 2013 & time <= 2024) %>% 
  glimpse()

print(pdo_clean) # resulting data has breeding season PDO values, non-breeding season, as well as 1 year lag of each

saveRDS(pdo_clean, "data/tidy/pdo_tidy.rds")



### END ###