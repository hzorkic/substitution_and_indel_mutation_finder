library(tidyverse)
library(csv)
library(dplyr)

GISAID <- read.csv("GISAID.csv")


df = data.frame(GISAID) 

count(df, c("type", "mutations"))


head(GISAID)

