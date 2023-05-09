library(spectrolab)
library(ggplot2)
library(patchwork)

setwd("C:/Users/querc/Dropbox/PostdocProjects/LeafSenescence")

###########################################
## read in data

spectra_tagged<-readRDS("processed_data/spectra_tagged.rds")

###########################################
## calculate some basic indices

## change to English for plots
levels(meta(spectra_tagged)$site_id)[levels(meta(spectra_tagged)$site_id)=="MSB_Tourbiere"]<-"MSB_Bog"

meta(spectra_tagged)$ARI<-1/spectra_tagged[,550]-1/spectra_tagged[,700]
meta(spectra_tagged)$PSRI<-(spectra_tagged[,680]-spectra_tagged[,500])/spectra_tagged[,750]
meta(spectra_tagged)$GMI<-(spectra_tagged[,750]-spectra_tagged[,705])/(spectra_tagged[,750]+spectra_tagged[,705])
meta(spectra_tagged)$mND705<-(spectra_tagged[,750]-spectra_tagged[,705])/(spectra_tagged[,750]+spectra_tagged[,705]-2*spectra_tagged[,445])

colors<-c("Acer rubrum Linnaeus"="firebrick1",
          "Betula populifolia Marshall"="gold1")

ARI_plot<-ggplot(meta(spectra_tagged),
                 aes(x=date,y=ARI,color=species_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  scale_color_manual(values=colors)+
  theme_bw()+theme(text=element_text(size=20),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = c(0.25,0.7),
                   legend.key.width = unit(1.5,"cm"),
                   legend.background = element_rect(fill="transparent"))+
  labs(x="Date",color="Species",shape="Site",linetype="Site",
       y="Anthocyanin Reflectance Index")

ggplot(meta(spectra_tagged),
       aes(x=date,y=PSRI,color=species_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  scale_color_manual(values=colors)+
  theme_bw()+theme(text=element_text(size=20),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = c(0.25,0.7),
                   legend.key.width = unit(1.5,"cm"),
                   legend.background = element_rect(fill="transparent"))+
  labs(x="Date",color="Species",shape="Site",linetype="Site",
       y="Plant Senescence Reflectance Index")

GMI_plot<-ggplot(meta(spectra_tagged),
                 aes(x=date,y=GMI,color=species_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  scale_color_manual(values=colors)+
  theme_bw()+theme(text=element_text(size=20),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = c(0.25,0.3),
                   legend.key.width = unit(1.5,"cm"),
                   legend.background = element_rect(fill="transparent"))+
  labs(x="Date",color="Species",shape="Site",linetype="Site",
       y="Chlorophyll Content Index")

ggplot(meta(spectra_tagged),
       aes(x=date,y=mND705,color=species_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  theme_bw()+theme(text=element_text(size=20),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = c(0.2,0.3))+
  labs(x="Date",color="Species",shape="Site",linetype="Site")

pdf("manuscript/ARI.pdf",width=10,height=7)
ARI_plot
dev.off()

pdf("manuscript/GMI.pdf",width=10,height=7)
GMI_plot
dev.off()

####################################
## summarize spectra

# spectra_Hermine_Acer<-spectra_tagged[meta(spectra_tagged)$species_id=="Acer rubrum Linnaeus" & meta(spectra_tagged)$site_id=="MSB_Hermine",]
# spectra_Hermine_Betula<-spectra_tagged[meta(spectra_tagged)$species_id=="Betula populifolia Marshall" & meta(spectra_tagged)$site_id=="MSB_Hermine",]
# spectra_Tourbiere_Acer<-spectra_tagged[meta(spectra_tagged)$species_id=="Acer rubrum Linnaeus" & meta(spectra_tagged)$site_id=="MSB_Tourbiere",]
# spectra_Tourbiere_Betula<-spectra_tagged[meta(spectra_tagged)$species_id=="Betula populifolia Marshall" & meta(spectra_tagged)$site_id=="MSB_Tourbiere",]
# 
# spectraHA_agg<-aggregate(as.matrix(spectra_Hermine_Acer),
#                          by=list(meta(spectra_Hermine_Acer)$date),FUN=median)
# spectraHA_agg_long<-melt(spectraHA_agg,id.vars = "Group.1")
# spectraHA_agg_long$variable<-as.numeric(as.character(spectraHA_agg_long$variable))
# spectraHA_agg_long$Group.1<-as.factor(spectraHA_agg_long$Group.1)
# 
# spectraHA<-ggplot(spectraHA_agg_long,aes(x=variable,y=value,color=Group.1))+
#   geom_line(size=1)+
#   theme_bw()+
#   coord_cartesian(ylim=c(0,0.55))+
#   labs(title = "Hermine Acer")
# 
# spectraHB_agg<-aggregate(as.matrix(spectra_Hermine_Betula),
#                          by=list(meta(spectra_Hermine_Betula)$date),FUN=median)
# spectraHB_agg_long<-melt(spectraHB_agg,id.vars = "Group.1")
# spectraHB_agg_long$variable<-as.numeric(as.character(spectraHB_agg_long$variable))
# spectraHB_agg_long$Group.1<-as.factor(spectraHB_agg_long$Group.1)
# 
# spectraHB<-ggplot(spectraHB_agg_long,aes(x=variable,y=value,color=Group.1))+
#   geom_line(size=1)+
#   theme_bw()+
#   coord_cartesian(ylim=c(0,0.55))+
#   labs(title = "Hermine Betula")
# 
# spectraTA_agg<-aggregate(as.matrix(spectra_Tourbiere_Acer),
#                          by=list(meta(spectra_Tourbiere_Acer)$date),FUN=median)
# spectraTA_agg_long<-melt(spectraTA_agg,id.vars = "Group.1")
# spectraTA_agg_long$variable<-as.numeric(as.character(spectraTA_agg_long$variable))
# spectraTA_agg_long$Group.1<-as.factor(spectraTA_agg_long$Group.1)
# 
# spectraTA<-ggplot(spectraTA_agg_long,aes(x=variable,y=value,color=Group.1))+
#   geom_line(size=1)+
#   theme_bw()+
#   coord_cartesian(ylim=c(0,0.55))+
#   labs(title = "Tourbiere Acer")
# 
# spectraTB_agg<-aggregate(as.matrix(spectra_Tourbiere_Betula),
#                          by=list(meta(spectra_Tourbiere_Betula)$date),FUN=median)
# spectraTB_agg_long<-melt(spectraTB_agg,id.vars = "Group.1")
# spectraTB_agg_long$variable<-as.numeric(as.character(spectraTB_agg_long$variable))
# spectraTB_agg_long$Group.1<-as.factor(spectraTB_agg_long$Group.1)
# 
# spectraTB<-ggplot(spectraTB_agg_long,aes(x=variable,y=value,color=Group.1))+
#   geom_line(size=1)+
#   theme_bw()+
#   coord_cartesian(ylim=c(0,0.55))+
#   labs(title = "Tourbiere Betula")
# 
# (spectraTA/spectraTB)|(spectraHA/spectraHB)

spectra_tagged_long<-melt(as.data.frame(spectra_tagged),
                           id.vars = c("sample_name",colnames(meta(spectra_tagged))))
spectra_tagged_long$variable<-as.numeric(as.character(spectra_tagged_long$variable))

pdf("manuscript/spectra_panels.pdf",height=20,width=16)
ggplot(spectra_tagged_long,aes(x=variable,y=value))+
  stat_summary(geom="line",fun="median",size=1.5,aes(color=site_id))+
  facet_grid(week~species_id)+
  theme_bw()+theme(text=element_text(size=20))
dev.off()