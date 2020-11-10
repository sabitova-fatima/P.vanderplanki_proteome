# Polypedilum vanderplanki proteome exploratory analysis

## ЭТАПЫ РАБОТЫ ##

## 1. Очистка результатов

1. Удаление NA
2. Удаление сомнительных результатов
3. Удаление ненужных для анализа столбцов

## 2. Статистическая обработка

Для анализа использовался пакет «limma» в среде R. В частности, была
применена функция, которая проводит t-test путем эмпирической байесовской модерации стандартных отклонений. Значения дифференциальной экспрессии при применении этого метода приведены к логарифмическому виду, что облегчает их анализ: значения ниже нуля показывают снижение концентрации, а значения выше – его повышение. За статистически достоверный результат были приняты данные с p-value ниже 0.05.

## 3. Анализ данных
Объект - клеточная линия pv11, полученная от эмбрионов комара P. vanderplanki. Обработана трегалозойЮ контрольная группа - клетки без обработки.


## 4. Визуализация

Графики сделаны в среде R с помощью пакета ggplot2


