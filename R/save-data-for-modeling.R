## clean genomic data
library(dplyr)
library(readr)
library(MMWRweek)
library(ggplot2)

#file of national level human genomic data samples. Filtered in the file "filtering genomic file.RMD"
national_level_raw <- read_csv("practicing_stan/national_level.csv") |> 
  filter(!is.na(date)) %>% #filter out rows with no date 
  mutate(year = MMWRweek(date)$MMWRyear) |> ## using this year function as it handles the dates at shoulders better
  mutate(epiweek = MMWRweek(date)$MMWRweek) 

d <- national_level_raw |> 
  group_by(year, epiweek, Nextstrain_clade) |> 
  summarize(clade_samples = n()) |> 
  ungroup() |> 
  group_by(year, epiweek) |> 
  mutate(total_samples = sum(clade_samples)) |> 
  ungroup() |> 
  mutate(epiweek_year = paste(year, formatC(epiweek, width=2, flag="0"), sep="_"),
         epidate = MMWRweek2Date(year, epiweek),
         clade_pct = clade_samples/total_samples)

write_csv(d, file="data/genomic_data_for_modeling.csv")

# Visualize data
ggplot(d, aes(x=epidate, y=clade_pct, color=Nextstrain_clade))+
  geom_line()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


