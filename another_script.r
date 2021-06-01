setwd("C:/Users/User/Desktop/proteome coursework")

library(limma)
library(tidyverse)
library(readxl)
library(ggplot2)
library(esquisse)

####### Import #########

yusurika <- read_excel("19c15301_Proteome_yusurika_data_2019 working.xlsx", sheet = 2) %>%
  filter(is.na(`average ratio`) != TRUE)

new <- read_excel("New assembly annotation.xlsx")

yusurika %>%
  left_join(new, by = c('Description' = 'Transcript')) -> The_Proteome

####### Cleaning the data #########

The_Proteome <- select(The_Proteome, - Protein, - Gene) 

# Переименовываем столбцы
mutate(The_Proteome, T0_1 = `Abundances Scaled F1 Sample yusurika_T0`, 
       T0_2 = `Abundances Scaled F2 Sample yusurika_T0`,
       T0_3 = `Abundances Scaled F3 Sample yusurika_T0`,
       T0_4 = `Abundances Scaled F4 Sample yusurika_T0`, 
       T24_1 = `Abundances Scaled F5 Sample yusurika_T24`,
       T24_2 = `Abundances Scaled F6 Sample yusurika_T24`,
       T24_3 = `Abundances Scaled F7 Sample yusurika_T24`,
       T24_4 = `Abundances Scaled F8 Sample yusurika_T24`) %>%
  select(-c(`Abundances Scaled F1 Sample yusurika_T0`, `Abundances Scaled F2 Sample yusurika_T0`,
            `Abundances Scaled F3 Sample yusurika_T0`,
            `Abundances Scaled F4 Sample yusurika_T0`, 
            `Abundances Scaled F5 Sample yusurika_T24`,
            `Abundances Scaled F6 Sample yusurika_T24`,
            `Abundances Scaled F7 Sample yusurika_T24`,
            `Abundances Scaled F8 Sample yusurika_T24`)) -> The_Proteome

# Удаляем пропущенные значения
sum(is.na(The_Proteome$T0_1)) # 64 

#Вместо удаления всех белков, где есть хотя бы один NA в T0 - T24, заменяем эти NA на среднее значение других 3 проб

#T0_1
for (i in 1:length(The_Proteome$T0_1)) {
  if (is.na(The_Proteome$T0_1[i])) {
    The_Proteome$T0_1[i] = mean(c(The_Proteome$T0_2[i], The_Proteome$T0_3[i], The_Proteome$T0_4[i]), na.rm = TRUE)
  }
}

#T0_2
for (i in 1:length(The_Proteome$T0_2)) {
  if (is.na(The_Proteome$T0_2[i])) {
    The_Proteome$T0_2[i] = mean(c(The_Proteome$T0_1[i], The_Proteome$T0_3[i], The_Proteome$T0_4[i]), na.rm = TRUE)
  }
}

#T0_3
for (i in 1:length(The_Proteome$T0_3)) {
  if (is.na(The_Proteome$T0_3[i])) {
    The_Proteome$T0_3[i] = mean(c(The_Proteome$T0_2[i], The_Proteome$T0_1[i], The_Proteome$T0_4[i]), na.rm = TRUE)
  }
}

#T0_4
for (i in 1:length(The_Proteome$T0_4)) {
  if (is.na(The_Proteome$T0_4[i])) {
    The_Proteome$T0_4[i] = mean(c(The_Proteome$T0_2[i], The_Proteome$T0_3[i], The_Proteome$T0_1[i]), na.rm = TRUE)
  }
}

#T24_1
for (i in 1:length(The_Proteome$T24_1)) {
  if (is.na(The_Proteome$T24_1[i])) {
    The_Proteome$T24_1[i] = mean(c(The_Proteome$T24_2[i], The_Proteome$T24_3[i], The_Proteome$T24_4[i]), na.rm = TRUE)
  }
}

#T24_2
for (i in 1:length(The_Proteome$T24_2)) {
  if (is.na(The_Proteome$T24_2[i])) {
    The_Proteome$T24_2[i] = mean(c(The_Proteome$T24_1[i], The_Proteome$T24_3[i], The_Proteome$T24_4[i]), na.rm = TRUE)
  }
}

#T24_3
for (i in 1:length(The_Proteome$T24_3)) {
  if (is.na(The_Proteome$T24_3[i])) {
    The_Proteome$T24_3[i] = mean(c(The_Proteome$T24_2[i], The_Proteome$T24_1[i], The_Proteome$T24_4[i]), na.rm = TRUE)
  }
}

#T24_4
for (i in 1:length(The_Proteome$T24_4)) {
  if (is.na(The_Proteome$T24_4[i])) {
    The_Proteome$T24_4[i] = mean(c(The_Proteome$T24_2[i], The_Proteome$T24_3[i], The_Proteome$T24_1[i]), na.rm = TRUE)
  }
}

sum(is.na(The_Proteome$T0_1)) # 0 - NA не осталось

####### Stats #########

glimpse(The_Proteome)
glimpse(yusurika)

design <- model.matrix(~ 0+factor(c(1,1,1,1,2,2,2,2)))
colnames(design) <- c("T0", "T24")
fit <- lmFit(The_Proteome[21:28], design)
contrast.matrix <- makeContrasts(T24-T0, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit.data <- as.data.frame(fit2)

glimpse(fit.data)

volcanoplot(fit2, pch=16, cex=0.5, coef = 1L, style = "p-value", col = heat.colors(10))

# fit.data - dataframe со статистическими данными. Т.к. гены расположены в той же последовательности, добавляем названия генов простым присоединением. 
# nrow(fit.data) == nrow(The_Proteome)  ----   TRUE (4428)

The_Proteome <- bind_cols(The_Proteome, fit.data)
glimpse(The_Proteome) # теперь это данные о белках + статистика 

# добавляем столбцы для отбора белков по p-value
The_Proteome$Good_p_value <- The_Proteome$F.p.value < 0.05
The_Proteome$p_value <- The_Proteome$Good_p_value

for (i in 1:length(The_Proteome$p_value)){
  if (The_Proteome$p_value[i] == TRUE){
    The_Proteome$p_value[i] <- "< 0.05"
  }else{
    The_Proteome$p_value[i] <- "> 0.05"
  }
}

# и столбец для удобства построения графиков
The_Proteome$minus_log10pvalue <- -log10(The_Proteome$F.p.value)


####### Exploratory analysis #########

#Смотрим, сколько p-value ниже 0,05
nrow(filter(The_Proteome, Good_p_value == TRUE)) #2233
nrow(filter(The_Proteome, Good_p_value == FALSE)) #2195

#то же на графике
p <- data.frame('P_value' = c('less_than_0,05', 'more_than_0,05'), 
                'Count' = c(2233, 2195))

lbls <- c('P_value < 0,05 (good)', 'P_value > 0,05 (bad)')
pie(p$Count, labels = lbls, main="Pie Chart of P-values")

# и на еще одном графике
ggplot(The_Proteome) +
  aes(x = p_value, fill = p_value) +
  geom_bar() +
  labs(title = "P-value distribution in all proteins") +
  theme_minimal()


# фильтруем, оставляя только белки с p value < 0.05
The_Proteome_f <- filter(The_Proteome, F.p.value < 0.05)
glimpse(The_Proteome)
nrow(The_Proteome_f) #2233 белка с хорошим P-value осталось


#смотрим распределение fold change
ggplot(The_Proteome_f) +
  aes(x = coefficients) +
  geom_histogram(bins = 90L, fill = "#0c4c8a") +
  theme_minimal()

#График со всеми 4428 белками
ggplot(The_Proteome, aes(coefficients, min_log10pvalue)) +
  geom_point() 

# Volcano plot со всеми белками
ggplot(The_Proteome) +
  aes(x = coefficients, y = minus_log10pvalue) +
  geom_point(size = 1L, colour = "#cb181d") +
  labs(x = "log2 Fold Change", y = "-log10 P-value", title = "All proteins volcano plot") +
  theme_minimal()

# Насколько fold change зависит от длины? 
ggplot(The_Proteome, aes(coefficients, Length))+
  geom_point()

# А от количества пептидов? 
ggplot(The_Proteome, aes(coefficients, `Number of Peptides by Search Engine Sequest HT`))+
  geom_point() # зависимость такая же 

# чтобы подтвердить, что все норм, должен получиться график логарифма
ggplot(The_Proteome, aes(coefficients, `average ratio`))+
  geom_point() # да, все работает нормально

# хромосомки :)
ggplot(The_Proteome, aes(`Scaffold Id`))+
  geom_bar(aes(fill=`Scaffold Id`))+
  scale_colour_brewer(palette = "Set2")# что за странные значения?

ggplot(The_Proteome, aes(x=Description, y=coefficients)) + 
  geom_point()  +
  coord_flip()

####### Full analysis. Часть 1. Анализ функций #########

###### Анализ по Interpro ########

The_Proteome %>%
  filter(Good_p_value == TRUE)%>%
  filter(Interpro != "NA")%>%
  filter(Interpro != "")%>%
  filter(Interpro != " ")%>%
  select(Description, coefficients, Interpro) %>%
  na.omit() -> Interpro

t <- as.data.frame(str_split_fixed(Interpro$Interpro, "\n", 22))

Interpro <- bind_cols(Interpro, t)
Interpro <- select(Interpro, -'Interpro')

# Я убила 2 часа на попытку сделать все одним циклом, но увы
# мне жаль

V1 <- Interpro[,c(1,2,3)]
V1 <- select(V1, Description, coefficients, functions = V1)
V2 <- Interpro[,c(1,2,4)]
V2 <- select(V2, Description, coefficients, functions = V2)
V3 <- Interpro[,c(1,2,5)]
V3 <- select(V3, Description, coefficients, functions = V3)
V4 <- Interpro[,c(1,2,6)]
V4 <- select(V4, Description, coefficients, functions = V4)
V5 <- Interpro[,c(1,2,7)]
V5 <- select(V5, Description, coefficients, functions = V5)
V6 <- Interpro[,c(1,2,8)]
V6 <- select(V6, Description, coefficients, functions = V6)
V7 <- Interpro[,c(1,2,9)]
V7 <- select(V7, Description, coefficients, functions = V7)
V8 <- Interpro[,c(1,2,10)]
V8 <- select(V8, Description, coefficients, functions = V8)
V9 <- Interpro[,c(1,2,11)]
V9 <- select(V9, Description, coefficients, functions = V9)
V10 <- Interpro[,c(1,2,12)]
V10 <- select(V10, Description, coefficients, functions = V10)
V11 <- Interpro[,c(1,2,13)]
V11 <- select(V11, Description, coefficients, functions = V11) 
V12 <- Interpro[,c(1,2,14)]
V12 <- select(V12, Description, coefficients, functions = V12)
V13 <- Interpro[,c(1,2,15)]
V13 <- select(V13, Description, coefficients, functions = V13)
V14 <- Interpro[,c(1,2,16)]
V14 <- select(V14, Description, coefficients, functions = V14)
V15 <- Interpro[,c(1,2,17)]
V15 <- select(V15, Description, coefficients, functions = V15)
V16 <- Interpro[,c(1,2,18)]
V16 <- select(V16, Description, coefficients, functions = V16)
V17 <- Interpro[,c(1,2,19)]
V17 <- select(V17, Description, coefficients, functions = V17)
V18 <- Interpro[,c(1,2,20)]
V18 <- select(V18, Description, coefficients, functions = V18)
V19 <- Interpro[,c(1,2,21)]
V19 <- select(V19, Description, coefficients, functions = V19)
V20 <- Interpro[,c(1,2,22)]
V20 <- select(V20, Description, coefficients, functions = V20)
V21 <- Interpro[,c(1,2,23)]
V21 <- select(V21, Description, coefficients, functions = V21)
V22 <- Interpro[,c(1,2,24)]
V22 <- select(V22, Description, coefficients, functions = V22) 


Interpro <- bind_rows(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, 
                      V13, V14, V15, V16, V17, V18, V19, V20, V21, V22) 

Interpro <- filter(Interpro, functions != "")
nrow(Interpro) # 7892
sum(is.na(Interpro)) # 0

# превращаем вот это
# IPR021418: THO complex, subunitTHOC2, C-terminal (871:1172)\r\r"
# В это 
# "IPR021418: THO complex, subunitTHOC2, C-terminal "
# да, очень странным способом. Опрделенно, stringr умеет лучше, но я не умею (4 утра)

t <- as.data.frame(str_split_fixed(Interpro$functions, "\\(", 2))
t <- select(t, -"V2")
Interpro <- bind_cols(Interpro, t)
Interpro <- select(Interpro, -'functions', functions = V1) 

# ура! Мы получили красивый датасет с функциями и белками, такой же, как в 'Molecular function'
#$ Description  <chr> "g1000.t1", "g10026.t1"
#$ coefficients <dbl> -14.4500, 44.4250, -11.7250, 
#$ functions    <fct> "IPR021418: THO complex, subunitTHOC2, C-terminal "

# Mean fold change
mean_fc <- Interpro %>% 
  group_by(functions) %>% summarise(mean = mean(coefficients), count = n())

# Готовим данные для красивых графиков
mean_fc <- mean_fc %>% arrange(mean) 
mean_fc <- mean_fc[order(mean_fc$mean), ] 
mean_fc$functions <- factor(mean_fc$functions, levels = mean_fc$functions)
# делаем z-преобразование
mean_fc$mean_z <- round((mean_fc$mean - mean(mean_fc$mean))/sd(mean_fc$mean), 2)

mean_fc$type <- ifelse(mean_fc$mean_z < 0, "Downregulated", "Upregulated")

# Diverging Lollipop Chart
ggplot(mean_fc[1:100, ], aes(x = functions, y = mean_z, label = mean_z))+ 
  geom_point(stat='identity', fill="black", size=5)  +
  geom_segment(aes(y = 0, 
                   x = functions, 
                   yend = mean_z, 
                   xend = functions), 
               color = "black") +
  geom_text(color="white", size=2)+ 
  labs(subtitle="Первые 100 наблюдений", 
       title= "Средний fold change белков для каждой функции") + 
  ylim(-5, 2.5) +
  coord_flip()

nrow(mean_fc) # 3465 неужели так много разных функций?! 

###### Анализ по Molecular function ########

##### делаем простой dafaframe с description, coefficiens и GO Molecular function 
# это не так просто, потому что функции засунуты все в одну ячейку, поэтому кода будет много

The_Proteome %>%
  filter(Good_p_value == TRUE)%>%
  filter(`Molecular Function` != "NA")%>%
  select(Description, coefficients, `Molecular Function`) %>%
  na.omit() %>% rename( molecular_function = `Molecular Function`) -> mol_func

t <- as.data.frame(str_split_fixed(mol_func$molecular_function, "/", 11))

molecular_function <- bind_cols(mol_func, t)
molecular_function <- select(molecular_function, -'molecular_function')

# В оригинальной таблице функции каждого белка были указаны в одной ячейке через /
# В каждой ячейке было от 1 до 10 функций белка в виде GO:0000166 nucleotide binding
# Я разделила ячейки на 10 разных столбцов, которые снова привела в один столбец

# да, я опять не пишу циклы. Каюсь

V1 <- molecular_function[,1:3]
V1 <- rename(V1, functions = V1)
V2 <- molecular_function[,c(1,2,4)]
V2 <- V2[str_which(molecular_function$V2, 'GO:'),]
V2 <- rename(V2, functions = V2)
V3 <- molecular_function[,c(1,2,5)]
V3 <- V3[str_which(molecular_function$V3, 'GO:'),]
V3 <- rename(V3, functions = V3)
V4 <- molecular_function[,c(1,2,6)]
V4 <- V4[str_which(molecular_function$V4, 'GO:'),]
V4 <- rename(V4, functions = V4)
V5 <- molecular_function[,c(1,2,7)]
V5 <- V5[str_which(molecular_function$V5, 'GO:'),]
V5 <- rename(V5, functions = V5)
V6 <- molecular_function[,c(1,2,8)]
V6 <- V6[str_which(molecular_function$V6, 'GO:'),]
V6 <- rename(V6, functions = V6)
V7 <- molecular_function[,c(1,2,9)]
V7 <- V7[str_which(molecular_function$V7, 'GO:'),]
V7 <- rename(V7, functions = V7)
V8 <- molecular_function[,c(1,2,10)]
V8 <- V8[str_which(molecular_function$V8, 'GO:'),]
V8 <- rename(V8, functions = V8)
V9 <- molecular_function[,c(1,2,11)]
V9 <- V9[str_which(molecular_function$V9, 'GO:'),]
V9 <- rename(V9, functions = V9)
V10 <- molecular_function[,c(1,2,12)]
V10 <- V10[str_which(molecular_function$V10, 'GO:'),]
V10 <- rename(V10, functions = V10)
V11 <- molecular_function[,c(1,2,13)]
V11 <- V11[str_which(molecular_function$V11, 'GO:'),]
V11 <- rename(V11, functions = V11) 
V11 <- slice(V11, c(2,3)) # удалила странную строку 

Simple_proteome <- bind_rows(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11) %>% arrange(functions)

# Ура! Мы получили Simple_proteome - датасет с теми самыми тремя колонками, пойдем его анализировать!

## здесь молек функция + среднее foldchange по всем белкам, которые выполняют эту функцию, всего 620 функций

mean_fold_change <- Simple_proteome %>% 
  transmute(func = functions, coef = coefficients) %>%
  group_by(func) %>% summarise(mean = mean(coef))


# Готовим данные для красивых графиков
mean_fold_change_arr <- mean_fold_change %>% arrange(mean) 
mean_fold_change_arr <- mean_fold_change_arr[order(mean_fold_change_arr$mean), ] 
mean_fold_change_arr$func <- factor(mean_fold_change_arr$func, levels = mean_fold_change_arr$func)
# делаем z-преобразование
mean_fold_change_arr$mean_z <- round((mean_fold_change_arr$mean - mean(mean_fold_change_arr$mean))/sd(mean_fold_change_arr$mean), 2)

mean_fold_change_arr$type <- ifelse(mean_fold_change_arr$mean_z < 0, "Downregulated", "Upregulated")

# крошечные боксплоты
ggplot(mean_fold_change_arr[1:100,]) +
  aes(x = func, y = mean) +
  geom_boxplot(fill = "#0c4c8a") +
  theme_minimal()+
  labs(subtitle="Первые 100 наблюдений", 
       title= "Средний fold change белков для каждой функции") + 
  coord_flip()


# Diverging Lollipop Chart
ggplot(mean_fold_change_arr[1:100, ], aes(x = func, y = mean_z, label = mean_z))+ 
  geom_point(stat='identity', fill="black", size=6)  +
  geom_segment(aes(y = 0, 
                   x = func, 
                   yend = mean_z, 
                   xend = func), 
               color = "black") +
  geom_text(color="white", size=2)+ 
  labs(subtitle="Первые 100 наблюдений", 
       title= "Средний fold change белков для каждой функции") + 
  ylim(-2.5, 2.5) +
  coord_flip()
# увы, наблюдений слишком много, чтобы вместить их в один график

###

# посмотрим, какие функции выполняют самые сильно изменившиеся белки
upreg <- filter(Simple_proteome, coefficients > 100)
downreg <- filter(Simple_proteome, coefficients < -100)


u <- ggplot(upreg, aes(coefficients))+
  geom_histogram(aes(fill=functions), 
                 bins=5, 
                 col="black", 
                 size=.1) +
  labs(title="Upregulated proteins", 
       subtitle="log2FoldChange > 100")  

d <- ggplot(downreg, aes(coefficients))+
  geom_histogram(aes(fill=functions), 
                 bins=5, 
                 col="black", 
                 size=.1) +
  labs(title="Downregulated proteins", 
       subtitle="log2FoldChange > 100")  

###### Анализ по Cellular Component ########

The_Proteome %>%
  filter(Good_p_value == TRUE)%>%
  filter(`Cellular Component` != "NA")%>%
  select(Description, coefficients, `Cellular Component`) %>%
  na.omit() %>% rename(cellular_comp = `Cellular Component`) -> cell_comp

t <- as.data.frame(str_split_fixed(cell_comp$cellular_comp, "/", 3))

cellular_component <- bind_cols(cell_comp, t)
cellular_component <- select(cellular_component, -'cellular_comp') #560 строк

# В оригинальной таблице функции каждого белка были указаны в одной ячейке через /
# В каждой ячейке было от 1 до 10 функций белка в виде GO:0000166 nucleotide binding
# Я разделила ячейки на 10 разных столбцов, которые снова привела в один столбец

V1 <- cellular_component[,1:3]
V1 <- rename(V1, functions = V1)

V2 <- cellular_component[,c(1,2,4)]
V2 <- V2[str_which(cellular_component$V2, 'GO:'),]
V2 <- rename(V2, functions = V2)

V3 <- cellular_component[,c(1,2,5)]
V3 <- V3[str_which(cellular_component$V3, 'GO:'),]
V3 <- rename(V3, functions = V3)

Simple_proteome_cell_comp <- bind_rows(V1, V2, V3) %>% arrange(functions)

mean_fold_change_cc <- Simple_proteome_cell_comp %>% 
  transmute(func = functions, coef = coefficients) %>%
  group_by(func) %>% summarise(mean = mean(coef))

# Готовим данные для красивых графиков
mean_fold_change_cc_arr <- mean_fold_change_cc %>% arrange(mean) 
mean_fold_change_cc_arr <- mean_fold_change_cc_arr[order(mean_fold_change_cc_arr$mean), ] 
mean_fold_change_cc_arr$func <- factor(mean_fold_change_cc_arr$func, levels = mean_fold_change_cc_arr$func)
# делаем z-преобразование
mean_fold_change_cc_arr$mean_z <- round((mean_fold_change_cc_arr$mean - mean(mean_fold_change_cc_arr$mean))/sd(mean_fold_change_cc_arr$mean), 2)

mean_fold_change_cc_arr$type <- ifelse(mean_fold_change_cc_arr$mean_z < 0, "Downregulated", "Upregulated")

# Diverging Lollipop Chart
ggplot(mean_fold_change_cc_arr[1:100, ], aes(x = func, y = mean_z, label = mean_z))+ 
  geom_point(stat='identity', fill="black", size=6)  +
  geom_segment(aes(y = 0, 
                   x = func, 
                   yend = mean_z, 
                   xend = func), 
               color = "black") +
  geom_text(color="white", size=2)+ 
  labs(subtitle="Первые 100 наблюдений", 
       title= "Средний fold change белков для каждой функции") + 
  ylim(-2.5, 2.5) +
  coord_flip()

# посмотрим, какие функции выполняют самые сильно изменившиеся белки
upreg <- filter(Simple_proteome_cell_comp, coefficients > 100)
downreg <- filter(Simple_proteome_cell_comp, coefficients < -100)


u <- ggplot(upreg, aes(coefficients))+
  geom_histogram(aes(fill=functions), 
                 bins=5, 
                 col="black", 
                 size=.1) +
  labs(title="Upregulated proteins", 
       subtitle="log2FoldChange > 100")  

d <- ggplot(downreg, aes(coefficients))+
  geom_histogram(aes(fill=functions), 
                 bins=5, 
                 col="black", 
                 size=.1) +
  labs(title="Downregulated proteins", 
       subtitle="log2FoldChange > 100") 


###### Анализ по Biological process ########

The_Proteome %>%
  filter(Good_p_value == TRUE)%>%
  filter(`Biological Process` != "NA")%>%
  select(Description, coefficients, `Biological Process`) %>%
  na.omit() %>% rename(biological_process = `Biological Process`) -> biol_proc

t <- as.data.frame(str_split_fixed(biol_proc$biological_process, "/", 6))

biological_process <- bind_cols(biol_proc, t)
biological_process <- select(biological_process, -'biological_process') #902 строки

# В оригинальной таблице функции каждого белка были указаны в одной ячейке через /
# В каждой ячейке было от 1 до 10 функций белка в виде GO:0000166 nucleotide binding
# Я разделила ячейки на 10 разных столбцов, которые снова привела в один столбец

V1 <- biological_process[,1:3]
V1 <- rename(V1, functions = V1)
V2 <- biological_process[,c(1,2,4)]
V2 <- V2[str_which(biological_process$V2, 'GO:'),]
V2 <- rename(V2, functions = V2)
V3 <- biological_process[,c(1,2,5)]
V3 <- V3[str_which(biological_process$V3, 'GO:'),]
V3 <- rename(V3, functions = V3)
V4 <- biological_process[,c(1,2,6)]
V4 <- V4[str_which(biological_process$V4, 'GO:'),]
V4 <- rename(V4, functions = V4)
V5 <- biological_process[,c(1,2,7)]
V5 <- V5[str_which(biological_process$V5, 'GO:'),]
V5 <- rename(V5, functions = V5)
V6 <- biological_process[,c(1,2,8)]
V6 <- V6[str_which(biological_process$V6, 'GO:'),]
V6 <- rename(V6, functions = V6)

Simple_biol_proc <- bind_rows(V1, V2, V3, V4, V5, V6) %>% arrange(functions)
# 521 unique functions

# Средний fold change по функциям
mean_fc_bp <- Simple_biol_proc %>% 
  transmute(functions, coefficients, chr) %>%
  group_by(functions) %>% summarise(mean = mean(coefficients))

# Готовим данные для красивых графиков
mean_fc_bp_arr <- mean_fc_bp %>% arrange(mean) 
mean_fc_bp_arr <- mean_fc_bp_arr[order(mean_fc_bp_arr$mean), ] 
mean_fc_bp_arr$functions <- factor(mean_fc_bp_arr$functions, levels = mean_fc_bp_arr$functions)
mean_fc_bp_arr$mean_z <- round((mean_fc_bp_arr$mean - mean(mean_fc_bp_arr$mean))/sd(mean_fc_bp_arr$mean), 2) # делаем z-преобразование
mean_fc_bp_arr$type <- ifelse(mean_fc_bp_arr$mean_z < 0, "Downregulated", "Upregulated")

# Diverging Lollipop Chart
ggplot(mean_fc_bp_arr[1:100, ], aes(x = functions, y = mean_z, label = mean_z))+ 
  geom_point(stat='identity', fill="black", size=6)  +
  geom_segment(aes(y = 0, 
                   x = functions, 
                   yend = mean_z, 
                   xend = functions), 
               color = "black") +
  geom_text(color="white", size=2)+ 
  labs(subtitle="Первые 100 наблюдений", 
       title= "Средний fold change белков для каждой функции") + 
  ylim(-2.5, 2.5) +
  coord_flip()

# посмотрим, какие функции выполняют самые сильно изменившиеся белки
upreg_100 <- filter(Simple_biol_proc, coefficients > 100) # 29
downreg_100 <- filter(Simple_biol_proc, coefficients < -100) # 25

upreg_50 <- filter(Simple_biol_proc, coefficients > 50) # 106
downreg_50 <- filter(Simple_biol_proc, coefficients < -50) # 127

u_100 <- ggplot(upreg_100, aes(coefficients))+
  geom_histogram(aes(fill=functions), 
                 bins=5, 
                 col="black", 
                 size=.1) +
  labs(title="Upregulated proteins", 
       subtitle="log2FoldChange > 100")  

d_100 <- ggplot(downreg_100, aes(coefficients))+
  geom_histogram(aes(fill=functions), 
                 bins=5, 
                 col="black", 
                 size=.1) +
  labs(title="Downregulated proteins", 
       subtitle="log2FoldChange > 100") 

u_50 <- ggplot(upreg_50, aes(coefficients))+
  geom_histogram(aes(fill=functions), 
                 bins=5, 
                 col="black", 
                 size=.1) +
  labs(title="Upregulated proteins", 
       subtitle="log2FoldChange > 100")  

d_50 <- ggplot(downreg_50, aes(coefficients))+
  geom_histogram(aes(fill=functions), 
                 bins=5, 
                 col="black", 
                 size=.1) +
  labs(title="Downregulated proteins", 
       subtitle="log2FoldChange > 100") 


# Посмотрим по хромосомам
The_Proteome %>%
  filter(Good_p_value == TRUE)%>%
  filter(`Biological Process` != "NA")%>%
  filter(`Scaffold Id` %in% c('chr_1', 'chr_2', 'chr_3', 'chr_4'))%>%
  select(Description, chr =`Scaffold Id`) -> The_cut_proteome_bf

Simple_biol_proc %>% left_join(The_cut_proteome_bf) -> Simple_biol_proc

chr_1_bp <- filter(Simple_biol_proc, chr == "chr_1") # 454 rows
chr_2_bp <- filter(Simple_biol_proc, chr == "chr_2") # 518 rows
chr_3_bp <- filter(Simple_biol_proc, chr == "chr_3") # 614 rows
chr_4_bp <- filter(Simple_biol_proc, chr == "chr_4") # 137 rows

# Хромосома 1
mean_fc_bp_1 <- chr_1_bp %>% 
  transmute(functions, coefficients, chr) %>%
  group_by(functions) %>% summarise(mean = mean(coefficients))

# Готовим данные для красивых графиков
mean_fc_bp_1 <- mean_fc_bp_1 %>% arrange(mean) 
mean_fc_bp_1 <- mean_fc_bp_1[order(mean_fc_bp_1$mean), ] 
mean_fc_bp_1$functions <- factor(mean_fc_bp_1$functions, levels = mean_fc_bp_1$functions)
mean_fc_bp_1$mean_z <- round((mean_fc_bp_1$mean - mean(mean_fc_bp_1$mean))/sd(mean_fc_bp_1$mean), 2) # делаем z-преобразование
mean_fc_bp_1$type <- ifelse(mean_fc_bp_1$mean_z < 0, "Downregulated", "Upregulated")


# Diverging Lollipop Chart
ggplot(mean_fc_bp_1[1:150, ], aes(x = functions, y = mean_z, label = mean_z))+ 
  geom_point(stat='identity', fill="black", size=5)  +
  geom_segment(aes(y = 0, 
                   x = functions, 
                   yend = mean_z, 
                   xend = functions), 
               color = "black") +
  geom_text(color="white", size=2)+ 
  labs(subtitle="Первые 100 наблюдений", 
       title= "Средний fold change белков для каждой функции") + 
  ylim(-3, 2.5) +
  coord_flip()

# Хромосома 2
mean_fc_bp_2 <- chr_2_bp %>% 
  transmute(functions, coefficients, chr) %>%
  group_by(functions) %>% summarise(mean = mean(coefficients))

# Готовим данные для красивых графиков
mean_fc_bp_2 <- mean_fc_bp_2 %>% arrange(mean) 
mean_fc_bp_2 <- mean_fc_bp_2[order(mean_fc_bp_2$mean), ] 
mean_fc_bp_2$functions <- factor(mean_fc_bp_2$functions, levels = mean_fc_bp_2$functions)
mean_fc_bp_2$mean_z <- round((mean_fc_bp_2$mean - mean(mean_fc_bp_2$mean))/sd(mean_fc_bp_2$mean), 2) # делаем z-преобразование
mean_fc_bp_2$type <- ifelse(mean_fc_bp_2$mean_z < 0, "Downregulated", "Upregulated")


# Diverging Lollipop Chart
ggplot(mean_fc_bp_2[1:150, ], aes(x = functions, y = mean_z, label = mean_z))+ 
  geom_point(stat='identity', fill="black", size=5)  +
  geom_segment(aes(y = 0, 
                   x = functions, 
                   yend = mean_z, 
                   xend = functions), 
               color = "black") +
  geom_text(color="white", size=2)+ 
  labs(subtitle="Первые 100 наблюдений", 
       title= "Средний fold change белков для каждой функции") + 
  ylim(-3, 2.5) +
  coord_flip()

# Хромосома 3
mean_fc_bp_3 <- chr_3_bp %>% 
  transmute(functions, coefficients, chr) %>%
  group_by(functions) %>% summarise(mean = mean(coefficients))

# Готовим данные для красивых графиков
mean_fc_bp_3 <- mean_fc_bp_3 %>% arrange(mean) 
mean_fc_bp_3 <- mean_fc_bp_3[order(mean_fc_bp_3$mean), ] 
mean_fc_bp_3$functions <- factor(mean_fc_bp_3$functions, levels = mean_fc_bp_3$functions)
mean_fc_bp_3$mean_z <- round((mean_fc_bp_3$mean - mean(mean_fc_bp_3$mean))/sd(mean_fc_bp_3$mean), 2) # делаем z-преобразование
mean_fc_bp_3$type <- ifelse(mean_fc_bp_3$mean_z < 0, "Downregulated", "Upregulated")


# Diverging Lollipop Chart
ggplot(mean_fc_bp_3[1:150, ], aes(x = functions, y = mean_z, label = mean_z))+ 
  geom_point(stat='identity', fill="black", size=5)  +
  geom_segment(aes(y = 0, 
                   x = functions, 
                   yend = mean_z, 
                   xend = functions), 
               color = "black") +
  geom_text(color="white", size=2)+ 
  labs(subtitle="Первые 100 наблюдений", 
       title= "Средний fold change белков для каждой функции") + 
  ylim(-3, 2.5) +
  coord_flip()

# Хромосома 4
mean_fc_bp_4 <- chr_4_bp %>% 
  transmute(functions, coefficients, chr) %>%
  group_by(functions) %>% summarise(mean = mean(coefficients))

# Готовим данные для красивых графиков
mean_fc_bp_4 <- mean_fc_bp_4 %>% arrange(mean) 
mean_fc_bp_4 <- mean_fc_bp_4[order(mean_fc_bp_4$mean), ] 
mean_fc_bp_4$functions <- factor(mean_fc_bp_4$functions, levels = mean_fc_bp_4$functions)
mean_fc_bp_4$mean_z <- round((mean_fc_bp_4$mean - mean(mean_fc_bp_4$mean))/sd(mean_fc_bp_4$mean), 2) # делаем z-преобразование
mean_fc_bp_4$type <- ifelse(mean_fc_bp_4$mean_z < 0, "Downregulated", "Upregulated")


# Diverging Lollipop Chart
ggplot(mean_fc_bp_4, aes(x = functions, y = mean_z, label = mean_z))+ 
  geom_point(stat='identity', fill="black", size=8)  +
  geom_segment(aes(y = 0, 
                   x = functions, 
                   yend = mean_z, 
                   xend = functions), 
               color = "black") +
  geom_text(color="white", size=2)+ 
  labs(subtitle="Первые 100 наблюдений", 
       title= "Средний fold change белков для каждой функции") + 
  ylim(-3, 2.5) +
  coord_flip()

####### Full analysis. Часть 2. Отдельно по важным функциям #########

####### SELECTING #######

sum(str_detect(The_Proteome$Interpro, "LEA")) # 3 
sum(str_detect(The_Proteome$Interpro, "HSP")) # 28 
sum(str_detect(The_Proteome$Interpro, "heat")) # 7
sum(str_detect(The_Proteome$Interpro, "transporter")) # 89 
sum(str_detect(The_Proteome$Interpro, "Protease")) # 4  
sum(str_detect(The_Proteome$Interpro, "protease")) # 23 
sum(str_detect(The_Proteome$Interpro, "Protease Inhibitor")) # 0 REMOVE
sum(str_detect(The_Proteome$Interpro, "Protease inhibitor")) # 1 REMOVE
sum(str_detect(The_Proteome$Interpro, "protease inhibitor")) # 0 REMOVE
sum(str_detect(The_Proteome$Interpro, "Apoptosis")) # 2 
sum(str_detect(The_Proteome$Interpro, "apoptosis")) # 4 
sum(str_detect(The_Proteome$Interpro, "Transcription Factor")) # 2 
sum(str_detect(The_Proteome$Interpro, "Transcription factor")) # 40
sum(str_detect(The_Proteome$Interpro, "transcription factor")) # 16
sum(str_detect(The_Proteome$Interpro, "Protein Kinase")) # 0 
sum(str_detect(The_Proteome$Interpro, "Protein kinase")) # 127
sum(str_detect(The_Proteome$Interpro, "protein kinase")) # 83
sum(str_detect(The_Proteome$Interpro, "Ubiquitin")) # 109
sum(str_detect(The_Proteome$Interpro, "ubiquitin")) # 30
sum(str_detect(The_Proteome$Interpro, "DNA repair")) # 6
sum(str_detect(The_Proteome$Interpro, "Signal Transduction")) # 0
sum(str_detect(The_Proteome$Interpro, "Signal transduction")) # 0
sum(str_detect(The_Proteome$Interpro, "signal transduction")) # 0
sum(str_detect(The_Proteome$Interpro, "Ribosomal Protein")) # 2
sum(str_detect(The_Proteome$Interpro, "Ribosomal")) # 196
sum(str_detect(The_Proteome$Interpro, "Actin")) # 33
sum(str_detect(The_Proteome$Interpro, "actin")) # 66
sum(str_detect(The_Proteome$Interpro, "Globin")) # 1
sum(str_detect(The_Proteome$Interpro, "globin")) # 0
sum(str_detect(The_Proteome$Interpro, "Cytochrome Oxidase")) # 0
sum(str_detect(The_Proteome$Interpro, "Cytochrome oxidase")) # 1
sum(str_detect(The_Proteome$Interpro, "cytochrome oxidase")) # 0
sum(str_detect(The_Proteome$Interpro, "NADH dehydrogenase")) # 18

# 1. HSP
# 2. Transporters
# 3. Protease
# 4. Apoptosis
# 5. Transcriptional factor
# 6. Ubiquitin
# 7. DNA repair
# 8. Ribosomal Proteins
# 9. Actin
# 10. NADH dehydrogenase
# 11. Protein kinase
# 12. Antioxidants

####### Make CSV files to analyze them in python. See proteome_plots.pdf #######

HSP <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "HSP")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(HSP, file = "HSP.csv")
HSP$func <- rep("HSP", length(HSP$Description))

transporter <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "transporter")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(transporter, file = "transporter.csv")
transporter$func <- rep("transporter", length(transporter$Description))

protease <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "protease")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(protease, file = "protease.csv")
protease$func <- rep("protease", length(protease$Description))

apoptosis <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "apoptosis")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(apoptosis, file = "apoptosis.csv")
apoptosis$func <- rep("apoptosis", length(apoptosis$Description))

Transcription_factor <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "Transcription factor")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(Transcription_factor, file = "Transcription_factor.csv")
Transcription_factor$func <- rep("Transcription_factor", length(Transcription_factor$Description))

Ubiquitin <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "Ubiquitin")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(Ubiquitin, file = "Ubiquitin.csv")
Ubiquitin$func <- rep("Ubiquitin", length(Ubiquitin$Description))

DNA_repair <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "DNA repair")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(DNA_repair, file = "DNA_repair.csv")
DNA_repair$func <- rep("DNA_repair", length(DNA_repair$Description))

Ribosomal <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "Ribosomal")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(Ribosomal, file = "Ribosomal.csv")
Ribosomal$func <- rep("Ribosomal", length(Ribosomal$Description))

Actin <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "Actin")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(Actin, file = "Actin.csv")
Actin$func <- rep("Actin", length(Actin$Description))

Protein_kinase <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "Protein kinase")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(Protein_kinase, file = "Protein_kinase.csv")
Protein_kinase$func <- rep("Protein_kinase", length(Protein_kinase$Description))

NADH_dehydrogenase <- filter(The_Proteome, str_detect(The_Proteome$Interpro, "NADH dehydrogenase")) %>% 
  unite(Description, Description, Good_p_value, remove = FALSE)
write.csv(NADH_dehydrogenase, file = "NADH_dehydrogenase.csv")
NADH_dehydrogenase$func <- rep("NADH_dehydrogenase", length(NADH_dehydrogenase$Description))

All <- bind_rows(HSP, transporter, protease, apoptosis, Transcription_factor, Ubiquitin, DNA_repair, Ribosomal, Protein_kinase, NADH_dehydrogenase)

write_csv(All, "All.csv") # для анализа в питоне

mean_all <- All %>% 
  transmute(func, coefficients) %>%
  group_by(func) %>% summarise(mean = mean(coefficients))

write_csv(mean_all, "mean_all.csv") # для анализа в питоне

####### cryptogenes #######

read_excel('cryptogenes.xlsx', sheet = 1) %>%
  select(genes = `New ID 3.0`, functions = "Query") %>%  
  left_join(The_Proteome, by = c("genes" = "Pv.09 Assembly Counterpart")) %>% 
  unite(Description, Description, functions, remove = FALSE) %>% 
  filter(Good_p_value == T) -> cryptogenes

write.csv(cryptogenes, "cryptogenes.csv")

# График - в питоне, сс proteome_plots.pdf






