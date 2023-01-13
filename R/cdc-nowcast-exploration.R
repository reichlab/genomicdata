## explore CDC variant nowcast data
## downloaded from here: https://data.cdc.gov/Laboratory-Surveillance/SARS-CoV-2-Variant-Proportions/jr58-6ysp

library(tidyverse)
library(lubridate)

v <- read_csv("data/20230113-SARS-CoV-2_Variant_Proportions.csv") |> 
  mutate(week_ending = mdy(substr(week_ending,start = 0, stop=10)),
         published_date = mdy(substr(published_date,start = 0, stop=10)),
         share_lo = as.numeric(share_lo))

## all estimates for USA, faceted by publication week
v |> 
  filter(usa_or_hhsregion == "USA") |> 
  ggplot(aes(x = week_ending, color = variant)) +
  geom_line(aes(y=share)) +
  facet_wrap(.~published_date, dir="v") +
  geom_vline(aes(xintercept = published_date), linetype=2) +
  theme(legend.position = "none")

## all estimates for USA, faceted by variant
v |> 
  group_by(variant) |> 
  mutate(large_variant = any(share > 0.2)) |> 
  ungroup() |> 
  filter(usa_or_hhsregion == "USA", large_variant) |> 
  ggplot(aes(x = week_ending, color = published_date)) +
  geom_line(aes(y=share, group=published_date)) +
  facet_wrap(.~variant, dir="v") +
  theme(legend.position = "none")

## all estimates for region 1, faceted by variant
v |> 
  group_by(variant) |> 
  mutate(large_variant = any(share > 0.2)) |> 
  ungroup() |> 
  filter(usa_or_hhsregion == "1", large_variant) |> 
  ggplot(aes(x = week_ending, color = published_date)) +
  geom_line(aes(y=share, group=published_date)) +
  facet_wrap(.~variant, dir="v") +
  theme(legend.position = "none")
