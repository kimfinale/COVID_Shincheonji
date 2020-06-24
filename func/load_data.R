library(readxl)
library(tidyverse)
library(lubridate)
POP_DAEGU <- 2.465*1e6
POP_SHINCHEONJI_DAEGU <- 9334
AVG_HOUSEHOLD_SIZE_DAEGU <- 2.5
# 대구광역시 전체 243만 1,523명
# '20.04, KOSIS (행정안전부, 주민등록인구현황)
# d <- read_xlsx( "data/covid_daegu.xlsx", "incidence")
# DAT <- d[ , 1:5 ]
# DAT$date <- as.Date(d$date)
# nr <- nrow(DAT)
# names( DAT ) <- c( "date", "daily_total", "daily_shincheonji", "total", "shincheonji" )
# 
# DAT$daily_total[2:nr] <- DAT$total[2:nr] - DAT$total[1:(nr-1)]
# DAT$daily_shincheonji[2:nr] <- DAT$shincheonji[2:nr] - DAT$shincheonji[1:(nr-1)]
# 
# DAT$daily_low <- DAT$daily_total - DAT$daily_shincheonji
# DAT$low <- DAT$total - DAT$shincheonji
# 
# id <- which( DAT$daily_low < 0 )
# # 
# # DAT[ id, c("daily_shincheonji","daily_low")] <- NA
# d <- DAT
# DAT[ id, "daily_total"] <- d[ id, "daily_shincheonji"]
# DAT[ id, "daily_shincheonji"] <- d[ id, "daily_total"]
# DAT[ id, "daily_low"] <-  - d[ id, "daily_low"]
# 
# 
# dates_low_abnoraml <- c( as.Date("2020-02-22"), as.Date("2020-02-24"), as.Date("2020-02-25"), as.Date("2020-03-01"), as.Date("2020-03-06") ) 
# 
# DAT[ DAT$date %in% dates_low_abnoraml, c("daily_low") ] <- NA
# write_csv( DAT, "outputs/DAT_Temp.csv" )
# DAT <- read_csv( "outputs/DAT_Temp.csv" )
# 
# DAT$date <- mdy( DAT$date )

DAT <- read_csv( "outputs/DAT_Imputed.csv" )

