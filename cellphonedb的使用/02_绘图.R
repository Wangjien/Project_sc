# //--------------------------------------------------------------
# // 绘制cellphonedb结果
# // wangjien 20230620
# //--------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(patchwork)
library(dpyr)
library(tidyverse)
library(tidyr)

# 绘制interacting_pair number
plotInteractingPair1 = function(pvales){
    pvales_sub = pvales[,12:dim(pvales)[2]]
    statdf = as.data.frame(colSums(pvales_sub < 0.05))
    colnames(statdf) = 'names'
    statdf$indexb=str_replace(rownames(statdf),"^.*\\|","")
    statdf$indexa=str_replace(rownames(statdf),"\\|.*$","")
    rankname=sort(unique(statdf$indexa))
    #转成因子类型，画图时，图形将按照预先设置的顺序排列
    statdf$indexa=factor(statdf$indexa,levels = rankname)
    statdf$indexb=factor(statdf$indexb,levels = rankname)
    p = statdf%>%ggplot(aes(x=indexa,y=indexb,fill=number))+geom_tile(color="white")+
        geom_text(aes(label = round(number, 2)), color = "black") +
        scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,20))+
        scale_x_discrete("cluster 1 produces molecule 1")+
        scale_y_discrete("cluster 2 produces molecule 2")+
        theme_minimal()+
        theme(
            axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
            panel.grid = element_blank()
        )
    return(p)
}

# 绘制对称的interacting_pair热图，A-B = 10 B-A = 20 A-B = B-A = 10+20
plotInteractingPair2 = function(pvales){
    pvales_sub = pvales[,12:dim(pvales)[2]]
    statdf = as.data.frame(colSums(pvales_sub < 0.05))
    colnames(statdf) = 'names'
    statdf$indexb=str_replace(rownames(statdf),"^.*\\|","")
    statdf$indexa=str_replace(rownames(statdf),"\\|.*$","")
    rankname=sort(unique(statdf$indexa))    

    statdf$total_number=0
    for (i in 1:dim(statdf)[1]) {
    tmp_indexb=statdf[i,"indexb"]
    tmp_indexa=statdf[i,"indexa"]
        if (tmp_indexa == tmp_indexb) {
            statdf[i,"total_number"] = statdf[i,"number"]
        } else {
            statdf[i,"total_number"] = statdf[statdf$indexb==tmp_indexb & statdf$indexa==tmp_indexa,"number"]+
            statdf[statdf$indexa==tmp_indexb & statdf$indexb==tmp_indexa,"number"]
        }
    }
    # set factor
    statdf$indexa=factor(statdf$indexa,levels = rankname)
    statdf$indexb=factor(statdf$indexb,levels = rankname)
    # plot
    p = statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
        scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"))+
        geom_text(aes(label = round(total_number, 2)), color = "black") + 
        scale_x_discrete("cluster 1")+
        scale_y_discrete("cluster 2")+
        theme_minimal()+
        theme(
            axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
            panel.grid = element_blank()
        )
    return(p)
}



