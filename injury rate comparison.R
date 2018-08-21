#injury rate comparison graph####
library(ggplot2)
umk_rate <- read.csv("/Users/sai/OneDrive/Summer Project/data/GBD/_std_rate_peryear.csv")
umk_rate <- umk_rate[,colnames(umk_rate) %in% c('sex_name','location_name','year','val')]
SA_rate <- read.csv("/Users/sai/OneDrive/Summer Project/data/GBD/rate_injury_SA_2007-15.csv")
SA_rate <- SA_rate[,colnames(SA_rate) %in% c('sex_name','location_name','year','val')]
SA_rate$val <- SA_rate$val/10
both_rate <- rbind(umk_rate,SA_rate)
if(TRUE) both_rate <- both_rate[both_rate$year!=2015,]

png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_std_inj_rate_peryear_compare.png",sep = ""), width = 800, height = 760)
ggplot(both_rate, aes(x=year, y=val, col=sex_name, linetype=location_name))+
  geom_line()+
  ggtitle("Standardized injury-related mortality rate in \n uMkhanyakude vs national average")+
  ylab("Age-standardized injury-related mortality per 10,000 PY")+
  theme(text = element_text(size=16))
dev.off()
#6.3 (4.0-10.1), 25.8 (19.8-33.9) and 14.9 (11.9-18.9) 
# ggplot(global_rate)+
#   geom_line(aes(x=year, y= val, col=sex_name)) +
#   geom_hline(yintercept=c(25.8,14.9,6.3), color=c('blue','red','green'),linetype="dashed")
# 
# ggplot(SA_rate)+
#   geom_line(aes(x=year, y= val, col=sex_name)) +
#   geom_hline(yintercept=c(25.8,14.9,6.3), color=c('blue','red','green'),linetype="dashed")
# 
# 
# #injury rate comparison graph####
# library(ggplot2)
# global_rate <- read.csv("/Users/sai/OneDrive/Summer Project/data/GBD/rate_injury_global_2007-15.csv")
# SA_rate <- read.csv("/Users/sai/OneDrive/Summer Project/data/GBD/rate_injury_SA_2007-15.csv")
# global_rate$val <- global_rate$val/10
# SA_rate$val <- SA_rate$val/10
# both_rate <- rbind(global_rate,SA_rate)
# 
# #6.3 (4.0-10.1), 25.8 (19.8-33.9) and 14.9 (11.9-18.9) 
# ggplot(global_rate)+
#   geom_line(aes(x=year, y= val, col=sex_name)) +
#   geom_hline(yintercept=c(25.8,14.9,6.3), color=c('blue','red','green'),linetype="dashed")
# 
# ggplot(SA_rate)+
#   geom_line(aes(x=year, y= val, col=sex_name)) +
#   geom_hline(yintercept=c(25.8,14.9,6.3), color=c('blue','red','green'),linetype="dashed")
# 
# 
#   ggplot(both_rate, aes(x=year, y=val, col=sex_name))+
#     geom_line()+facet_wrap(~location_name)+
#     geom_hline(yintercept=rep(c(25.8,14.9,6.3),2), color=rep(c('blue','red','green'),2),linetype="dashed") 
#     
#     geom_hline(yintercept=c(25.8,14.9,6.3), color=c('blue','red','green'),linetype="dashed")
#     
#     ggplot(both_rate, aes(x=year, y=val, col=sex_name, linetype=location_name))+
#       geom_line()