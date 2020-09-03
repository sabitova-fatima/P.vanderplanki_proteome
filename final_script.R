setwd("C:/Users/User/Desktop/patients")

library(readxl)
library(qdapTools)
library(stringr)
library(readr)
library(xlsx)
library(lubridate)
library(dplyr)

temp <- list.files(pattern="*.docx")
patients <- data.frame("шифр" = 1:length(temp))
patients_in <- sapply(temp, read_docx)

################## ПОСТУПЛЕНИЕ #####################

for(i in 1:length(temp)){
patients$шифр[i] <- str_sub(temp[i], start = 1, end = -6)
patients$рост[i] <- str_extract(patients_in[[i]][which(str_detect(patients_in[[i]], "Рост"))[1]], "[:digit:][:digit:][:digit:]")
patients$вес[i] <- str_extract(patients_in[[i]][which(str_detect(patients_in[[i]], "Вес"))[1]], "[:digit:][:digit:]\\s")
patients$сахарный_диабет[i] <- str_sub(patients_in[[i]][which(str_detect(patients_in[[i]], "Сахарный"))[1]], start = 17, end = -1)
patients$Систолическое_давление[i] <- str_sub(str_extract(patients_in[[i]][which(str_detect(patients_in[[i]], "Систолическое"))[1]], "Систолическое\\sдавление..............."),  start = 24, end = 27)
patients$Диастолическое_давление[i] <- str_sub(str_extract((patients_in[[i]][which(str_detect(patients_in[[i]], "Диастолическое"))[1]]), "Диастолическое\\sдавление:.............."), start = 26, end = 28)
patients$ЧСС[i] <- str_sub(str_extract((patients_in[[i]][which(str_detect(patients_in[[i]], "ЧСС"))[1]]), "ЧСС:\\s........"), start = 6, end = 8)
patients$общее_состояние[i] <- str_sub(patients_in[[i]][which(str_detect(patients_in[[i]], "Общее состояние"))[1]], start = 18, end = 35)
patients$сознание[i] <- str_sub(str_extract(patients_in[[i]][which(str_detect(patients_in[[i]], "Общее состояние"))[1]], "Сознание.................."), start = 10, end = 16)
patients$ЧДД[i] <- str_sub(str_extract(patients_in[[i]][which(str_detect(patients_in[[i]], "ЧДД"))[1]], "ЧДД: ......."), start = 5, end = 8)
patients$spo2[i] <- str_sub(str_extract(patients_in[[i]][which(str_detect(patients_in[[i]], "SPO2"))[1]], "SPO2:....."), start = 6, end = -3)
patients$NEWS[i] <- str_extract(patients_in[[i]][which(str_detect(patients_in[[i]], "NEWS"))[1]], "NEWS.....")
patients$дата_рождения[i] <- str_sub(patients_in[[i]][which(str_detect(patients_in[[i]], "Дата рождения"))[1]], start = 16, end = 25) 
patients$возраст[i] <- str_sub(patients_in[[i]][which(str_detect(patients_in[[i]], "Дата рождения"))[1]], start = 27, end = 29) 
patients$дата_поступления[i] <- str_sub(patients_in[[i]][which(str_detect(patients_in[[i]], "Дата поступления"))[1]], start = 31, end = -8) 
patients$время_поступления[i] <- str_sub(patients_in[[i]][which(str_detect(patients_in[[i]], "Дата поступления"))[1]], start = 41, end = -3) 
  
}

################## КОЛИЧЕСТВО ДНЕЙ #####################

regex_date <- paste0("^[:digit:][:digit:].[:digit:][:digit:].[:digit:][:digit:][:digit:][:digit:]")

patients_days <- data.frame(шифр = NA, first = 1:length(patients_in), last = NA, n_days = NA)

for(n in 1:length(patients_in)){
  
  p <- patients_in[[n]]
  
  dates_index <- which(str_detect(patients_in[[n]], regex_date))
  
  patients_days$шифр[n] <- str_sub(temp[n], start = 1, end = 13)
  patients_days$first[n] <- str_sub(p[dates_index[1]], start = 1, end = 10)
  patients_days$last[n] <- str_sub(p[dates_index][length(dates_index)], start = 1, end = 10)
  patients_days$n_days[n] <- as.integer(dmy(str_sub(p[dates_index][length(dates_index)], start = 1, end = 10)) - dmy(str_sub(p[dates_index[1]], start = 1, end = 10)))
}


patients_days_for_print <- patients_days
colnames(patients_days_for_print) <- c("шифр","первый_день", "посл_день","всего_дней")


################## СТАЦИОНАР #####################

ALL_DATA <- data.frame(шифр = NA, дата = NA, номер_суток = NA, время = NA, тип = NA, общее_состояние = NA, сознание = NA, ЧДД = NA, Spo2 = NA, Сист_АД = NA, Диаст_АД = NA, ЧСС = NA, NEWS = NA)

# divider <- data.frame(шифр = "0", дата = "0", номер_суток = "0", время = "0", тип = "0", общее_состояние = "0", сознание = "0", ЧДД = "0", Spo2 = "0", Сист_АД = "0", Диаст_АД = "0", ЧСС = "0", NEWS = "0")

for(j in 1:length(temp)) {
  current_patient <- patients_in[[j]]
  
  key_words <- paste0(c("ПЕРВИЧНЫЙ ОСМОТР", "ОБХОД В РЕАНИМАЦИИ", "ОСМОТР В РЕАНИМАЦИИ", 
                        "ДНЕВНИК", "ОСМОТР ДЕЖУРНОГО ВРАЧА", "СОВМЕСТНЫЙ", "ОБХОД В РЕАНИМАЦИИ",
                        "ОСМОТР ЗАВЕДУЮЩЕГО ОТДЕЛЕНИЕМ"),collapse = '|')
  
  indexes <- which(str_detect(current_patient, key_words))
  
  to_delete <- c(10000)
  
  for(n in 1:(length(indexes)-1)){
    if(indexes[n] + 1 == indexes[n+1]){
      to_delete <- c(to_delete, n+1)}
  }
  indexes <- indexes[-to_delete]
  
  temporary <- data.frame(шифр = 1:(length(indexes)-1), дата = NA, номер_суток = NA, время = NA, тип = NA, общее_состояние = NA, сознание = NA, ЧДД = NA, Spo2 = NA, Сист_АД = NA, Диаст_АД = NA, ЧСС = NA, NEWS = NA)
  
  for(i in 1:(length(indexes)-1)){
    current_check <- current_patient[indexes[i]:indexes[i+1]]
    
    temporary$шифр <- str_sub(temp[j], start = 1, end = -6)
    temporary$дата[i] <- current_check[1] %>% str_sub(start = 1, end = 10)
    temporary$время[i] <- current_check[1] %>% str_sub(start = 12, end = 16)
    temporary$тип[i] <- current_check[1] %>% str_sub(start = 18, end = -1)
    temporary$номер_суток[i] <- as.integer(dmy(temporary$дата[i]) - dmy(patients_days$first[j]))
    temporary$ЧСС[i] <- str_sub(str_extract(current_check[which(str_detect(current_check, "ЧСС"))], "ЧСС:...."), start = 6, end = 8)[1]
    
    temporary$температура[i] <- str_extract(current_check[which(str_detect(current_check, "3[:digit:],[:digit:]"))][1], "3[:digit:],[:digit:]")
    temporary$общее_состояние[i] <- str_sub(current_check[which(str_detect(current_check, "Общее состояние"))[1]], start = 18, end = 32)
    temporary$сознание[i] <- str_sub(str_extract(current_check[which(str_detect(current_check, "Общее состояние"))[1]], "Сознание.................."), start = 10, end = 15)
    temporary$ЧДД[i] <- str_sub(str_extract(current_check[which(str_detect(current_check, "ЧДД"))[1]], "ЧДД: ......."), start = 5, end = 8)
    temporary$Spo2[i] <- str_sub(str_extract(current_check[which(str_detect(current_check, "SPO2"))[1]], "SPO2:....."), start = 6, end = -3)
    temporary$Сист_АД[i] <- str_sub(str_extract(current_check[which(str_detect(current_check, "Систолическое"))[1]], "Систолическое\\sдавление..............."),  start = 24, end = 27)
    temporary$Диаст_АД[i] <- str_sub(str_extract((current_check[which(str_detect(current_check, "Диастолическое"))[1]]), "Диастолическое\\sдавление:.............."), start = 26, end = 28)
    temporary$NEWS[i] <- str_sub(str_extract(current_check[which(str_detect(current_check, "NEWS"))[1]], "NEWS....."), start = 7, end = 8)
    
  }
  ALL_DATA <- bind_rows(ALL_DATA, temporary) 
}

ALL_DATA <- ALL_DATA[-1,]


################## СОХРАНЕНИЕ #####################
write.xlsx(patients, "поступление_6.xlsx")
write.xlsx(patients_days_for_print, "количество_дней_6.xlsx")
write.xlsx(ALL_DATA, "стационар6.xlsx")

