library(tidyverse)
library(here)
library(viridis)

facs<-read.csv(here("Ly6a_FACS_all_samples.csv"))
facs
facs_df<-data.frame(Ly6a_perc = c(facs$Matrigel_Ly6aPerc,facs$collagen_Ly6aPerc), 
                     sample=c(rep("matrigel", times=length(facs$Matrigel_Ly6aPerc)),
                              rep("collagen", times=length(facs$collagen_Ly6aPerc))))

facs_df<-facs_df[complete.cases(facs_df),]
facs_df
class(facs_df$Ly6a_perc[1])

facs_df$sample<-factor(facs_df$sample,
                          levels=c("matrigel", "collagen"))
ggplot(facs_df, aes(x=sample, y=Ly6a_perc))+
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean", fill="lightgrey", color="black")+
  geom_jitter(width = 0.2, size=4, alpha=0.5)+
  ylab("Ly6a positive cells [%]")+
  theme_classic()+
  theme(axis.title.x = element_blank())


organoids<-read.csv(here("Ly6a_sorting_organoid_growth.csv"), header = TRUE)
organoids$percent_outgrowth<-(organoids$organoids/organoids$cellnumber)*100

organoids$source_passage<-as.character(organoids$source_passage)

# save as pdf with dim 5 x 8 inch
barplot(1:10, col = viridis(10))
viridis(10)[c(1,5,10)]
organoids%>%
  ggplot(aes(x=Ly6a_stain, percent_outgrowth, color=source_passage))+
    geom_bar(position = "dodge", stat = "summary", fun.y = "mean", fill="lightgrey", color="black")+
    geom_point(size=5)+
    geom_point(shape = 1, size = 5,colour = "black")+
    scale_color_manual(values = viridis(10)[c(1,5,10)])+
    ylim(0,2.5)+
    theme_classic()+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
