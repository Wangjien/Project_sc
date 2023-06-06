library(ggplot2)
library(dplyr)
library(patchwork)

# 读入文件
df = read.table('zhang.txt', sep = '\t', header=T)

# r$> head(df)
#      time1 MolarMass.1    time2           LS    time3          RI
# 1 37.51058   359420315 34.02582 -0.000025600 34.02723 0.000569725
# 2 37.52741   352332147 34.04266  0.000009270 34.04302 0.000582702
# 3 37.54425   345304832 34.05949  0.000043300 34.05881 0.000645462
# 4 37.56108   338640082 34.07633  0.000078500 34.07460 0.000644269
# 5 37.57791   332398066 34.09316  0.000115006 34.09039 0.000663840
# 6 37.59475   326056835 34.11000  0.000152457 34.10618 0.000699025
# 宽数据转换为长数据
p1 = ggplot(data = df, aes(x = time1, y = MolarMass.1)) + 
    geom_line()+
    labs(title = 'Molar Mass')+
    theme_bw()

p2 = ggplot(data = df, aes(x = time2, y = LS)) + 
    geom_line(color = 'red')+
    labs(title='LS')+
    theme_bw()

p3 = ggplot(data = df, aes(x = time3, y = RI)) + 
    geom_line(color='blue')+
    labs(title  = 'RI')+
    theme_bw()

p = ggplot(data = df, aes(x = time, y =value))+
    geom_line() +
    facet_wrapper(.~group)   

p + scale_color_manual(values = c('#ff3434','#b6b6ff','black')) + 
    labs(x = 'Time (min)', y = expression(log10^('Molar Mass + 1')*'(g/mol)'))+ 
    theme(axis.title = element_text(size = 25, color = 'black'),
        axis.title.y = element_text(angle = 90),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_blank()
        ) + theme_bw()
