library(readxl)
library(tidyverse)
library(lubridate)
POP_DAEGU <- 2.465*1e6
POP_SHINCHEONJI_DAEGU <- 9334
AVG_HOUSEHOLD_SIZE_DAEGU <- 2.5
# 대구광역시 전체 243만 1,523명
# '20.04, KOSIS (행정안전부, 주민등록인구현황)
DAT <- read_csv( "outputs/DAT_Imputed.csv" )

