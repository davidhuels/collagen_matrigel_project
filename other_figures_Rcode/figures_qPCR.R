library(tidyverse)
library(here)
library(ggridges)
library(RColorBrewer)
lam<-read.csv(here("manuscript_figures/Figure_5/20201122_lamininIF_development_more_samples.csv"), header = TRUE, stringsAsFactors = FALSE)

tail(lam)
lam$time<-gsub(lam$tissue, pattern = "_[0-9]", replacement = "")
unique(lam$time)
lam$time<-factor(lam$time, levels = c("e14", "e18", "p7", "adult"))
head(lam)

lam%>%
  filter(time != "e14")%>%
  ggplot(aes(x=time, y=lam_intensity))+
    geom_boxplot(fill="lightgrey")+
    geom_jitter(width = 0.1, size=2.5)+
    ylab("Laminin intensity across crypt region")+
    ylim(0,260)+
    scale_x_discrete(breaks=c("e18", "p7", "adult"),
                   labels=c("E18", "P7", "Adult"))+
    theme_classic()+
    theme(axis.title.x = element_blank())

# use ggridges for distribution plotting
lam%>%
  filter(time != "e14")%>%
  ggplot(aes(x=lam_intensity, y=time))+
    geom_density_ridges()


epimes<-read.csv(here("06112020_gene expression laminins _P7_epi_mes.csv"))
head(epimes)
epimes$group<-gsub(epimes$tissue, pattern = "[0-9]", replacement = "")
epimes%>%
  ggplot(aes(x=gene, y=expression, colour=group))+
    geom_point()+
    scale_y_continuous(trans = 'log10')

epimes%>%
  group_by(gene, group)%>%
  summarise(average = mean(expression, na.rm = TRUE))%>%
  group_by(gene)%>%
  mutate(rel.expression = average/ average[group == "m"])%>%
  ggplot(aes(x=gene, y=rel.expression, fill=group))+
  geom_col(position = position_dodge())+
  coord_flip()+
  theme_bw()

#plot indvidual rel. expression and add bar chart later

epimes%>%
  group_by(gene)%>%
  mutate(av.m = mean(expression[group == "m"]),
         rel.expression.to.m = expression/av.m)%>%
  ggplot(aes(x=gene, y=rel.expression.to.m, colour=group, fill=group))+
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean")+
  geom_point(aes(x=gene, y=rel.expression.to.m, colour=group),
             position = position_dodge(width = 1), colour = "black")+
  coord_flip()+
  theme_bw()


colours_tSNE<-c("#d82e2e", "#884664", "#4876a4", "#498f8a", "#53a361",
                "#d97141", "#f79a4b", "#f5c643", "#bf9041", "#a25b3c", "#e684b5")
length(colours_tSNE)
barplot(1:11,col=colours_tSNE)
epimes%>%
  group_by(gene)%>%
  mutate(av.m = mean(expression[group == "m"]),
         rel.expression.to.m = expression/av.m)%>%
  ggplot(aes(x=gene, y=rel.expression.to.m, colour=group, fill=group))+
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean", 
           fill=rep(x=c("black","darkgrey"), times=9))+
  geom_point(aes(x=gene, y=rel.expression.to.m, colour=group),
             position = position_dodge(width = 1), colour = "black")+
  coord_flip()+
  theme_bw()


display.brewer.pal(n = 8, name = 'Greys')
display.brewer.pal(n=9, name = "Greys")
brewer.pal(n=9, name = "Greys")[c(4,7)]
head(epimes)
class(epimes$gene)
levels(epimes$gene)
epimes$gene<-factor(x=epimes$gene, levels = c("Col4a1", "Lama1", "Lama2", "Lama4", "Lama5",
                                              "Lamb1", "Lamb2", "Lamc1", "Intga6"))

epimes%>%
  group_by(gene)%>%
  mutate(av.m = mean(expression[group == "m"]),
         rel.expression.to.m = expression/av.m)%>%
  ggplot(aes(x=gene, y=rel.expression.to.m, group=group, fill=group))+
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean", 
           fill="grey")+
  geom_point(aes(x=gene, y=rel.expression.to.m, group=group),
             position = position_dodge(width = 1), colour = "black")+
  ylab("expression rel. to mesenchyme")+
  theme_bw()+
  theme(axis.title.x = element_blank())




epimes%>%
  group_by(gene)%>%
  mutate(av.m = mean(expression[group == "m"]),
         rel.expression.to.m = expression/av.m)%>%
  ggplot(aes(x=gene, y=rel.expression.to.m, group=group, fill=group))+
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean", 
           fill=c(rep(x="#BDBDBD", times=9), rep(x="#525252", times=9)))+
  geom_point(aes(x=gene, y=rel.expression.to.m, group=group),
             position = position_dodge(width = 1), colour = "black", size = 2)+
  ylab("expression rel. to mesenchyme")+
  theme_bw()+
  theme(axis.title.x = element_blank())


