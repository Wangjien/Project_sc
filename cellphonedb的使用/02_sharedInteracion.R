library(ktplots)
library(ggplot2)
library(dplyr)
library(patchwork)

pvals = read.delim('pvalues.csv', sep = ',', header=T, check.names=F)
tmp = pvals %>% select('interacting_pair',12:last_col())
# 转换为长数据框
tmp.long =  tmp %>% 
        tidyr::pivot_longer(cols=-1,
               names_to = "group",
               values_to = "value")
# 去除不显著的组合
tmp.long$sig = tmp.long$value < 0.05
tmp.long.sig = tmp.long %>% filter(value < 0.05)
tmp.long.sig = tmp.long.sig %>% mutate(
        inter1 = str_split_fixed('\\|')[1]
)

for(i in 1:20){
        tmp = as.data.frame(split_data[[i]])
        colnames(tmp) = c('id','filename','md5','size','state')
        print(tmp)
        filename = paste0('data','_',i,'.txt')
        print(filename)
        write.table(tmp,file = filename, row.names=F,sep = '\t', quote=F)
}