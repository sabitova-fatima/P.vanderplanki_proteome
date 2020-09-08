setwd("C:/Users/User/Desktop/NEW_COURSEWORK")

####### Adding packages #########

library(limma)
library(tidyverse)
library(readxl)
library(esquisse)
library(phyloseq)

####### Importing and analyzing #########

The_Proteome <- read_excel("fixed_proteome.xlsx") 
The_Proteome$coefficients <- log2(The_Proteome$fc)
write.csv(The_Proteome, "The_Proteome.csv")

####### Exploratory data analysis #########

#Смотрим, сколько p-value ниже 0,05
nrow(filter(The_Proteome, Good_p_value == TRUE)) #2233
nrow(filter(The_Proteome, Good_p_value == FALSE)) #2195

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
ggplot(The_Proteome) +
  aes(x = coefficients) +
  geom_histogram(bins = 90L, fill = "#0c4c8a") +
  theme_minimal()+
  ggtitle('Распределение сhange всех белков')+
  labs(x = "log2 Fold Change", y = "Количество")

ggplot(The_Proteome, aes(y = coefficients)) +
  geom_boxplot() +
  theme_minimal()+
  ggtitle('Распределение log2 fold change')+
  labs(x = "log2 Fold Change", y = "Количество")+
  coord_flip()

# Volcano plot со всеми белками
ggplot(The_Proteome) +
  aes(x = coefficients, y = minus_log10pvalue, color = Good_p_value) +
  geom_point(size = 1L) +
  labs(x = "log2 Fold Change", y = "-log10 P-value", title = "Volcano plot всех белков") +
  theme_minimal()+
  scale_colour_discrete(name  ="P-value",
                        breaks=c("FALSE", "TRUE"),
                        labels=c("> 0.05", "< 0.05"))

# Насколько fold change зависит от длины? 
ggplot(The_Proteome, aes(coefficients, Length))+
  geom_point(color = "#0c4c8a")+
  labs(x = "log2 Fold Change", y = "Длина белка в п.н.")+
  ggtitle('Зависимость изменения концентрации белка от его длины')+
  theme_minimal()

# А от количества пептидов? 
ggplot(The_Proteome, aes(coefficients, `Number of Peptides by Search Engine Sequest HT`))+
  geom_point() +
  labs(x = "log2 Fold Change", y = "-log10 P-value")+
  ggtitle('Зависимость изменения концентрации от количества пептидов')# зависимость такая же 

# чтобы подтвердить, что все норм, должен получиться график логарифма
ggplot(The_Proteome, aes(coefficients, `average ratio`))+
  geom_point() # да, все работает нормально

# хромосомки :)
The_Proteome_n <- filter(The_Proteome, `Scaffold Id` %in% c("chr_1", "chr_2", "chr_3", "chr_4"))

ggplot(The_Proteome_n, aes(`Scaffold Id`))+
  geom_bar(aes(fill=`Scaffold Id`))+
  scale_colour_brewer(palette = "Set2")+
  labs(x = "Хромосомы", y = "Количество обнаруженных белков")+
  ggtitle('Распределение обнаруженных белков по хромосомам')+
  theme_minimal()

####### Full analysis. Часть 1. Анализ функций #########

###### Анализ по Interpro ########

Interpro <- read.csv("Interpro.csv")

#$ Description  <chr> "g1000.t1", "g10026.t1"
#$ coefficients <dbl> -14.4500, 44.4250, -11.7250, 
#$ functions    <fct> "IPR021418: THO complex, subunitTHOC2, C-terminal "

# Mean fold change
mean_fc <- Interpro %>% 
  group_by(functions) %>% summarise(mean = mean(coefficients), count = n()) %>% arrange(mean) 

mean_fc <- filter(mean_fc, count > 1)

mean_fc <- mean_fc[order(mean_fc$mean), ] 
mean_fc$functions <- factor(mean_fc$functions, levels = mean_fc$functions)
# делаем z-преобразование
mean_fc$mean_z <- round((mean_fc$mean - mean(mean_fc$mean))/sd(mean_fc$mean), 2)
mean_fc$type <- ifelse(mean_fc$mean_z < 0, "Downregulated", "Upregulated")

# Diverging Lollipop Chart - переделать в питоне
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

nrow(mean_fc) # 3465 

# посмотрим, какие функции выполняют самые сильно изменившиеся белки
upreg <- filter(Interpro, coefficients > 100)
downreg <- filter(Interpro, coefficients < -100)

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

###### Анализ по Molecular function ########
Simple_proteome <- read.csv("Molecular_function.csv")

## здесь молек функция + среднее foldchange по всем белкам, которые выполняют эту функцию, всего 620 функций

mean_fold_change <- Simple_proteome %>% 
  transmute(func = functions, coef = coefficients) %>%
  group_by(func) %>% summarise(mean = mean(coef), count = n())

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

# посмотрим, какие функции выполняют самые сильно изменившиеся белки
upreg <- filter(Simple_proteome, coefficients > 100)
downreg <- filter(Simple_proteome, coefficients < -100)


u <- ggplot(upreg, aes(coefficients))+
  geom_histogram(aes(fill=functions), 
                 bins=5, 
                 col="black", 
                 size=.1) +
  labs(title="Upregulated proteins by their molecular functions", 
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
  group_by(func) %>% summarise(mean = mean(coef), count = n())

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
upreg <- filter(Simple_proteome_cell_comp, coefficients > 50)
downreg <- filter(Simple_proteome_cell_comp, coefficients < -50)


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
  transmute(functions, coefficients) %>%
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

####### Full analysis. Часть 2. Отдельно по важным функциям #########

####### 11 maybe-realated to anhydrobiosis functions #######

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
# 12. Antioxidants (later)

# Making CSV files to analyze them in python. See proteome_plots.pdf 

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

####### cryptogenes #######

read_excel('cryptogenes.xlsx', sheet = 1) %>%
  select(genes = `New ID 3.0`, Query = "Query") %>%  
  left_join(The_Proteome, by = c("genes" = "Pv.09 Assembly Counterpart")) %>% 
  filter(Good_p_value == T) -> cryptogenes

read_excel('cryptogenes.xlsx', sheet = 1) %>%
  select(genes = `New ID 3.0`, Query = "Query") %>%  
  left_join(The_Proteome, by = c("genes" = "Pv.09 Assembly Counterpart"))-> cryptogenes1

cryptogenes$func <- rep("cryptogene", length(cryptogenes$Description))

write.csv(cryptogenes, "cryptogenes.csv")

mean(cryptogenes$coefficients) # 99.77849
mean(The_Proteome$coefficients) #-3.775512

# The_Proteome_f - with good p-value (< 0,05)
mean(The_Proteome_f$coefficients) #-6.20365

Non_cry <- setdiff(select(mutate(The_Proteome, func = "cryptogene"), -`Pv.09 Assembly Counterpart`), select(cryptogenes, -Query, -genes)) %>% mutate(func = "non_cryptogene")

# :')
Please_dont_cry <- bind_rows(Non_cry, select(cryptogenes, -Query, -genes)) %>% 
  filter(coefficients > -600)

ggplot(Please_dont_cry, aes(coefficients, minus_log10pvalue, color = func))+
  geom_point()+
  scale_color_manual(values = c("#FC4E07", "#999999"))+
  xlab("-log10 p-value") + 
  ylab("log2 Fold Change") +
  ggtitle("Cryptogene proteins")

####### Антиоксиданты #########
# glutathione peroxidase
# thioredoxin
# oxidative stress
# superoxide dismutase
# glutathione S-transferase

sum(str_detect(The_Proteome$`Molecular Function`, "glutathione peroxidase"), na.rm = T) #1207 1208 1209
sum(str_detect(The_Proteome$`Molecular Function`, "thioredoxin"), na.rm = T) #1188 2827
sum(str_detect(The_Proteome$`Molecular Function`, "oxidative stress"), na.rm = T) # 0
sum(str_detect(The_Proteome$`Molecular Function`, "superoxide dismutase"), na.rm = T) # 4 1677 1678 3775 4046
sum(str_detect(The_Proteome$`Molecular Function`, "glutathione S-transferase"), na.rm = T) # 0

sum(str_detect(The_Proteome$Interpro, "glutathione peroxidase")) # 0
sum(str_detect(The_Proteome$Interpro, "thioredoxin")) #546 1771 2771 2899
sum(str_detect(The_Proteome$Interpro, "oxidative stress")) # 0
sum(str_detect(The_Proteome$Interpro, "superoxide dismutase")) # 1677 1678 2900 3775 4046
sum(str_detect(The_Proteome$Interpro, "glutathione S-transferase")) # 2 468 2693

antioxidants <- The_Proteome[c(1207, 1208, 1209, 1188, 2827, 1677, 1678, 3775, 4046, 546, 1771, 2771, 2899, 2900, 468, 2693),] 

antioxidants %>% unite(Description, Description, Good_p_value, remove = FALSE) -> antioxidants
antioxidants$func <- rep("antioxidants", length(antioxidants$Description))

write.csv(antioxidants, file = "antioxidants.csv")

Non_antiox <- setdiff(mutate(The_Proteome, func = "antioxidants"), antioxidants) %>% mutate(func = "non_antioxidants")

Compare <- bind_rows(Non_antiox, antioxidants) %>% 
  filter(coefficients > -600)

ggplot(Compare, aes(coefficients, minus_log10pvalue, color = func))+
  geom_point()+
  scale_color_manual(values = c("#FC4E07", "#999999"))+
  xlab("-log10 p-value") + 
  ylab("log2 Fold Change") +
  ggtitle("Antioxidants")

###### все  вместе #######

Special <- bind_rows(HSP, transporter, protease, apoptosis, Transcription_factor, Ubiquitin, DNA_repair, Ribosomal, Protein_kinase, NADH_dehydrogenase, cryptogenes, antioxidants)

The_Proteome_t <- unite(The_Proteome, Description, Description, Good_p_value, remove = FALSE)

Non_special <- setdiff(mutate(The_Proteome_t, func = "special"), select(Special, - Query, - genes)) %>% mutate(func = "non_special")

Compare1 <- bind_rows(Non_special, Special) %>% 
  filter(coefficients > -600)

ggplot(Compare1, aes(coefficients, minus_log10pvalue, color = func))+
  geom_point(stat="identity")+
  xlab("-log10 p-value") + 
  ylab("log2 Fold Change") +
  ggtitle("Proteins special and non-special for dehydration")

mean_all <- All %>% 
  transmute(func, coefficients) %>%
  group_by(func) %>% summarise(mean = mean(coefficients), count = n())

median_all <- All %>% 
  transmute(func, coefficients) %>%
  group_by(func) %>% summarise(median = median(coefficients), count = n())

write_csv(All, "All.csv") # для анализа в питоне
write_csv(mean_all, "mean_all.csv") # для анализа в питоне

###### сравниваем особенные гены с обычными #######

The_Proteome_f$func <- ""
Non_special <- setdiff(The_Proteome_f, All1)

Special <- separate(All, Description, c("Description", "m"), sep = "_") %>% select(-m, -func, -genes, -functions) %>% filter(F.p.value < 0.05)

mean_sp <- mean(Special$coefficients) # -3.104199
mean_non_sp <- mean(Non_special$coefficients) # 1.175383
mean_cry <- mean(cryptogenes)

######## РЕЗУЛЬТАТЫ  ########

# all the proteins
mean(The_Proteome$coefficients) #-0.04589744

# The proteins with good p-value (< 0,05)
mean(The_Proteome_f$coefficients) #-0.07248773

# cryptogenes
mean(cryptogenes$coefficients) # 10.85544

# antioxidants
mean(antioxidants$coefficients) # 1.056678


esquisse::esquisser(The_Proteome)
