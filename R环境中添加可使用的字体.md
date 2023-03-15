## R环境中添加可使用的字体

```R
# 安装 extrafont 包
install.packages("extrafont")
# 加载 extrafont 包
library(extrafont)
# 查看可以使用的字体
fonttable()
# 添加字体
font_import(paths = "/usr/share/fonts/truetype/msttcorefonts/", pattern = "Times_New_Roman.ttf")
font_import(paths = "/usr/share/fonts/truetype/msttcorefonts/", pattern = "Verdana.ttf")
font_import(paths = "/usr/share/fonts/truetype/msttcorefonts/", pattern = "Courier_New.ttf")
font_import(paths = "/usr/share/fonts/truetype/msttcorefonts/", pattern = "Georgia.ttf")
font_import(paths = "/usr/share/fonts/truetype/msttcorefonts/", pattern = "Comic_Sans_MS_Bold.ttf")
font_import(paths = "/usr/share/fonts/truetype/msttcorefonts/", pattern = "Webdings.ttf")

# 使用字体
df <- data.frame(
    x = 1,
    y = 1:4,
    label = c('Arial','Times New Roman', 'Courier','Verdana')
)
png('./查看不同的字体.png', height = 2000, width = 2000, res= 300)
ggplot(data = df,aes(x =x , y =y ))+
    geom_text(aes(label=label),family = c('arial','times new roman', 'courier','verdana'),size = 10)
dev.off()
```

[![pp1gEDg.md.png](https://s1.ax1x.com/2023/03/15/pp1gEDg.md.png)](https://imgse.com/i/pp1gEDg)