library(reshape2)
library(ggplot2)
library(lme4)
library(lmerTest)

setwd("C:/Users/querc/Dropbox/PostdocProjects/LeafSenescence")

###################################
## plot senescence dates

sCombined<-read.csv("senescence_combined.csv")
## change to English for figures
sCombined$site_id[which(sCombined$site_id=="MSB_Tourbiere")]<-"MSB_Bog"
sCombined_long<-melt(sCombined,id.vars = c("plant_id","species_id","site_id","leaf_id","sen_date"))
colnames(sCombined_long)<-c("plant_id","species_id","site_id","leaf_id",
                           "sen_date","date","presence")
sCombined_long$date<-gsub(pattern = "X",replacement="",x = sCombined_long$date)
sCombined_long$date<-gsub(pattern = "\\.",replacement="/",x = sCombined_long$date)

sCombined$sen_date<-as.Date(sCombined$sen_date,"%m/%d/%Y")
sCombined$sen_JD<-as.numeric(format(sCombined$sen_date,"%j"))
sCombined_long$date<-as.Date(sCombined_long$date,"%m/%d/%Y")
sCombined_long$sen_date<-as.Date(sCombined_long$sen_date,"%m/%d/%Y")

colors<-c("Acer rubrum Linnaeus"="firebrick1",
          "Betula populifolia Marshall"="gold1")

frac_leaves<-ggplot(sCombined_long,aes(x=date,y=presence,
                                       color=species_id,shape=site_id,
                                       linetype=site_id))+
  geom_point(stat="summary",fun="mean",size=4)+
  geom_line(stat="summary",fun="mean",size=2)+
  scale_color_manual(values=colors)+
  theme_bw()+theme(text=element_text(size=20),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = c(0.25,0.4),
                   legend.key.width = unit(1.5,"cm"),
                   legend.background = element_rect(fill="transparent"))+
  labs(y="Fraction of leaves remaining",x="Date",
       color="Species",shape="Site",linetype="Site")+
  ylim(c(0,1))

pdf("manuscript/LeafFraction.pdf",width=10,height=7)
frac_leaves
dev.off()

sen_JD_mixed<-lmer(sen_JD~species_id*site_id+(1|plant_id),data=sCombined)
summary(sen_JD_mixed)
difflsmeans(sen_JD_mixed)

sen_JD_anova<-aov(sen_JD~species_id*site_id+Error(plant_id),data=sCombined)
summary(sen_JD_anova)

avg_sen<-aggregate(sCombined[,c("sen_JD")],
                   by=list(sCombined$plant_id),
                   FUN=mean,na.rm=T)

################################
## calculate N resorption efficiency

plants<-read.csv("summary/plants.csv")
gestionnaire<-read.csv("summary/2021-LeafSenescence_gestionnaire - Labwork management sheet.csv")
plants$avg_sen<-avg_sen$x[match(plants$plant_id,avg_sen$Group.1)]

## change to English for figures
plants$site_id[plants$site_id=="MSB_Tourbiere"]<-"MSB_Bog"

## initial N per mass
N_init<-read.csv("traits/CN/c_n_initial.csv")
N_init_match<-sapply(plants$plant_id,function(x) grep(pattern = x,x = N_init$analysis_remarks))
plants$n_mass_init<-N_init$n_perc[N_init_match]

## initial LMA
LMA<-read.csv("traits/LMA/2021-LeafSenescence_Leaf_disks_water_samples_for_LMA - Sheet1.csv")
LMA$plant_id<-gestionnaire$Plant_id[match(LMA$Bulk_id,gestionnaire$bulk_id)]
LMA_init<-LMA[LMA$date_sampled %in% c("2021-08-17","2021-08-19"),]
LMA_init_match<-sapply(plants$plant_id,function(x) grep(pattern = x,x = LMA_init$plant_id))
plants$LMA_init<-LMA$leaf_mass_per_area_g_m2[LMA_init_match]
plants$n_area_init<-(plants$n_mass_init/100)*plants$LMA_init

## final N per mass
N_final<-read.csv("traits/CN/c_n_final.csv")
N_final_match<-sapply(plants$plant_id,function(x) grep(pattern = x,x = N_final$analysis_remarks))
plants$n_mass_final<-N_final$n_perc[N_final_match]

## final LMA
LMA_final<-read.csv("traits/LMA/2021-LeafSenescence_gestionnaire - Senesced leaf management.csv")
LMA_final$leaf_mass_per_area_g_m2<-LMA_final$LMA_dry_mass_g/(6*(7/2)^2*pi/1000000)
LMA_final_agg<-aggregate(LMA_final$leaf_mass_per_area_g_m2,
                         by=list(LMA_final$plant_id),
                         FUN=mean,na.rm=T)
colnames(LMA_final_agg)<-c("plant_id","leaf_mass_per_area_g_m2")
LMA_final_match<-sapply(plants$plant_id,function(x) grep(pattern = x,x = LMA_final_agg$plant_id))
plants$LMA_final<-LMA_final_agg$leaf_mass_per_area_g_m2[LMA_final_match]
plants$n_area_final<-(plants$n_mass_final/100)*plants$LMA_final

plants$n_resorp<-(plants$n_area_init-plants$n_area_final)/plants$n_area_init*100

n_resorp_plot<-ggplot(plants,aes(y=n_resorp,x=site_id,color=scientific_name))+
  geom_violin(size=2)+
  labs(x="Scientific name",y="N resorption efficiency (%)",color="Site")+
  theme_bw()+
  scale_color_manual(values=colors)+
  theme(text=element_text(size=20),
        legend.position=c(0.7,0.2),
        legend.background = element_rect(fill="transparent"))

ggplot(plants,aes(y=n_area_init,x=site_id,color=scientific_name))+
  geom_violin(size=2)+
  labs(x="Scientific name",
       y=expression(paste("Initial N (g ",m^-2,")")),color="Site")+
  theme_bw()+
  scale_color_manual(values=colors)+
  theme(text=element_text(size=20),
        legend.position=c(0.3,0.85),
        legend.background = element_rect(fill="transparent"))

n_area_final_plot<-ggplot(plants,aes(y=n_area_final,x=site_id,color=scientific_name))+
  geom_violin(size=2)+
  labs(x="Scientific name",
       y=expression(paste("Final N (g ",m^-2,")")),color="Site")+
  theme_bw()+
  scale_color_manual(values=colors)+
  theme(text=element_text(size=20),
        legend.position=c(0.7,0.8),
        legend.background = element_rect(fill="transparent"))

ggplot(plants,aes(y=n_mass_init,x=site_id,color=scientific_name))+
  geom_violin(size=2)+
  labs(x="Scientific name",
       y="Initial N (%)",color="Site")+
  theme_bw()+
  scale_color_manual(values=colors)+
  theme(text=element_text(size=20),
        legend.position=c(0.7,0.25),
        legend.background = element_rect(fill="transparent"))

ggplot(plants,aes(y=n_mass_final,x=site_id,color=scientific_name))+
  geom_violin(size=2)+
  labs(x="Scientific name",
       y="Final N (%)",color="Site")+
  theme_bw()+
  scale_color_manual(values=colors)+
  theme(text=element_text(size=20),
        legend.position=c(0.3,0.8),
        legend.background = element_rect(fill="transparent"))

ggplot(plants,aes(y=n_resorp,x=avg_sen,shape=site_id,color=scientific_name))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  labs(x="Average abscission date",
       y="N resorption efficiency (%)",
       color="Species",
       shape="Site")+
  theme_bw()+
  theme(text=element_text(size=20))

pdf("manuscript/n_resorp_plot.pdf",width=8,height=7)
n_resorp_plot
dev.off()

pdf("manuscript/n_area_final_plot.pdf",width=8,height=7)
n_area_final_plot
dev.off()

###############################
## read in carbon fractions from initial bulks

carbon_fractions<-read.csv("traits/CarbonFractions/carbon_fractions_bags.csv")
chem_samples<-read.csv("traits/CarbonFractions/leaf_chemistry_samples.csv")
bulk_samples<-read.csv("summary/bulk_leaf_samples.csv")

chem_samples$plant_id<-bulk_samples$plant_id[match(chem_samples$sample_id,bulk_samples$sample_id)]
carbon_fractions$plant_id<-chem_samples$plant_id[match(carbon_fractions$bottle_id,chem_samples$bottle_id)]
carbon_fractions<-carbon_fractions[!is.na(carbon_fractions$plant_id),]

plants$soluble_perc<-carbon_fractions$soluble_perc[match(plants$plant_id,carbon_fractions$plant_id)]
plants$hemicellulose_perc<-carbon_fractions$hemicellulose_perc[match(plants$plant_id,carbon_fractions$plant_id)]
plants$cellulose_perc<-carbon_fractions$cellulose_perc[match(plants$plant_id,carbon_fractions$plant_id)]
plants$lignin_perc<-carbon_fractions$lignin_perc[match(plants$plant_id,carbon_fractions$plant_id)]

ggplot(carbon_fractions,
       aes(y=hemicellulose_perc,
           x=scientific_name,
           color=site_id))+
  geom_violin(size=2)+
  labs(x="Scientific name",
       y="Hemicellulose (%)",color="Site")+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position=c(0.2,0.8))
