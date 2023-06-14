library(spectrolab)
library(ggplot2)
library(patchwork)
library(pls)
library(caret)

setwd("C:/Users/querc/Dropbox/PostdocProjects/LeafSenescence/")

#######################################
## read in spectral and trait data from bulk samples

## spectral data
spectra_bulk<-readRDS("processed_data/spectra_bulk.rds")

## aggregate by bulk id
spectra_bulk_agg<-aggregate(spectra_bulk,by=meta(spectra_bulk)$bulk_id,
                            FUN=try_keep_txt(mean))

## LMA and water-related traits
LMA_samples<-read.csv("traits/LMA/2021-LeafSenescence_Leaf_disks_water_samples_for_LMA - Sheet1.csv")

LMA_samples$equivalent_water_thickness_cm[which(is.na(LMA_samples$leafdisks_fresh_mass_g))]<-NA
LMA_samples$actual_leaf_dry_matter_content_perc[which(is.na(LMA_samples$leafdisks_fresh_mass_g))]<-NA
LMA_samples$actual_leaf_dry_matter_content_perc<-as.numeric(LMA_samples$actual_leaf_dry_matter_content_perc)

LMA_match<-match(meta(spectra_bulk_agg)$bulk_id,LMA_samples$Bulk_id)
meta(spectra_bulk_agg)$LMA<-LMA_samples$leaf_mass_per_area_g_m2[LMA_match]
meta(spectra_bulk_agg)$LDMC<-LMA_samples$actual_leaf_dry_matter_content_perc[LMA_match]*10
meta(spectra_bulk_agg)$EWT<-LMA_samples$equivalent_water_thickness_cm[LMA_match]*10

## pigments
pigments<-read.csv("traits/ChlCar/pigments_extracts.csv")
pigments<-pigments[-which(pigments$quality_flag_extract=="bad"),]
pigments<-pigments[which(pigments$sample_type=="sample"),]

pigments_match<-match(meta(spectra_bulk_agg)$bulk_id,pigments$vial_id)
meta(spectra_bulk_agg)$chla_mass<-pigments$chla_mg_g_disk_mass[pigments_match]
meta(spectra_bulk_agg)$chlb_mass<-pigments$chlb_mg_g_disk_mass[pigments_match]

## we'll use only total chlorophyll rather than a and b
meta(spectra_bulk_agg)$chl_mass<-meta(spectra_bulk_agg)$chla_mass+meta(spectra_bulk_agg)$chlb_mass
meta(spectra_bulk_agg)$car_mass<-pigments$carot_mg_g_disk_mass[pigments_match]

meta(spectra_bulk_agg)$chl_area<-meta(spectra_bulk_agg)$chl_mass*meta(spectra_bulk_agg)$LMA
meta(spectra_bulk_agg)$car_area<-meta(spectra_bulk_agg)$car_mass*meta(spectra_bulk_agg)$LMA

meta(spectra_bulk_agg)$sp_site_time<-paste(meta(spectra_bulk_agg)$site_id,
                                           meta(spectra_bulk_agg)$species_id,
                                           meta(spectra_bulk_agg)$week,
                                           sep="_")

###################################
## testing / training split

## stratified by species x site x time combination
train.sample <- createDataPartition(
  y = meta(spectra_bulk_agg)$sp_site_time,
  p = .75,
  list = FALSE
)

spectra_train<-spectra_bulk_agg[train.sample]
spectra_test<-spectra_bulk_agg[-train.sample]

LMA_model<-plsr(meta(spectra_train)$LMA~as.matrix(spectra_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA <- selectNcomp(LMA_model, method = "onesigma", plot = FALSE)
LMA_valid <- which(!is.na(meta(spectra_train)$LMA))
LMA_pred<-data.frame(bulk_id=meta(spectra_train)$bulk_id[LMA_valid],
                     species_id=meta(spectra_train)$species_id[LMA_valid],
                     site_id=meta(spectra_train)$site_id[LMA_valid],
                     date=meta(spectra_train)$date[LMA_valid],
                     plant_id=meta(spectra_train)$plant_id[LMA_valid],
                     measured=meta(spectra_train)$LMA[LMA_valid],
                     val_pred=LMA_model$validation$pred[,,ncomp_LMA])

LDMC_model<-plsr(meta(spectra_train)$LDMC~as.matrix(spectra_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LDMC <- selectNcomp(LDMC_model, method = "onesigma", plot = FALSE)
LDMC_valid <- which(!is.na(meta(spectra_train)$LDMC))
LDMC_pred<-data.frame(bulk_id=meta(spectra_train)$bulk_id[LDMC_valid],
                      species_id=meta(spectra_train)$species_id[LDMC_valid],
                      site_id=meta(spectra_train)$site_id[LDMC_valid],
                      date=meta(spectra_train)$date[LDMC_valid],
                      plant_id=meta(spectra_train)$plant_id[LDMC_valid],
                      measured=meta(spectra_train)$LDMC[LDMC_valid],
                      val_pred=LDMC_model$validation$pred[,,ncomp_LDMC])

EWT_model<-plsr(meta(spectra_train)$EWT~as.matrix(spectra_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT <- selectNcomp(EWT_model, method = "onesigma", plot = FALSE)
EWT_valid <- which(!is.na(meta(spectra_train)$EWT))
EWT_pred<-data.frame(bulk_id=meta(spectra_train)$bulk_id[EWT_valid],
                     species_id=meta(spectra_train)$species_id[EWT_valid],
                     site_id=meta(spectra_train)$site_id[EWT_valid],
                     date=meta(spectra_train)$date[EWT_valid],
                     plant_id=meta(spectra_train)$plant_id[EWT_valid],
                     measured=meta(spectra_train)$EWT[EWT_valid],
                     val_pred=EWT_model$validation$pred[,,ncomp_EWT])

chl_area_model<-plsr(sqrt(meta(spectra_train)$chl_area)~as.matrix(spectra_train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chl_area <- selectNcomp(chl_area_model, method = "onesigma", plot = FALSE)
chl_area_valid <- which(!is.na(meta(spectra_train)$chl_area))
chl_area_pred<-data.frame(bulk_id=meta(spectra_train)$bulk_id[chl_area_valid],
                          species_id=meta(spectra_train)$species_id[chl_area_valid],
                          site_id=meta(spectra_train)$site_id[chl_area_valid],
                          date=meta(spectra_train)$date[chl_area_valid],
                          plant_id=meta(spectra_train)$plant_id[chl_area_valid],
                          measured=meta(spectra_train)$chl_area[chl_area_valid],
                          val_pred=chl_area_model$validation$pred[,,ncomp_chl_area]^2)

car_area_model<-plsr(sqrt(meta(spectra_train)$car_area)~as.matrix(spectra_train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_area <- selectNcomp(car_area_model, method = "onesigma", plot = FALSE)
car_area_valid <- which(!is.na(meta(spectra_train)$car_area))
car_area_pred<-data.frame(bulk_id=meta(spectra_train)$bulk_id[car_area_valid],
                          species_id=meta(spectra_train)$species_id[car_area_valid],
                          site_id=meta(spectra_train)$site_id[car_area_valid],
                          date=meta(spectra_train)$date[car_area_valid],
                          plant_id=meta(spectra_train)$plant_id[car_area_valid],
                          measured=meta(spectra_train)$car_area[car_area_valid],
                          val_pred=car_area_model$validation$pred[,,ncomp_car_area]^2)

## apply to testing data to estimate model performance
LMA_val_pred<-data.frame(bulk=meta(spectra_test)$bulk_id,
                         val_pred=predict(LMA_model,newdata=as.matrix(spectra_test),ncomp=ncomp_LMA)[,,1],
                         measured=meta(spectra_test)$LMA)

LDMC_val_pred<-data.frame(bulk=meta(spectra_test)$bulk_id,
                         val_pred=predict(LDMC_model,newdata=as.matrix(spectra_test),ncomp=ncomp_LDMC)[,,1],
                         measured=meta(spectra_test)$LDMC)

EWT_val_pred<-data.frame(bulk=meta(spectra_test)$bulk_id,
                         val_pred=predict(EWT_model,newdata=as.matrix(spectra_test),ncomp=ncomp_EWT)[,,1],
                         measured=meta(spectra_test)$EWT)

chl_area_val_pred<-data.frame(bulk=meta(spectra_test)$bulk_id,
                                 val_pred=predict(chl_area_model,newdata=as.matrix(spectra_test),ncomp=ncomp_chl_area)[,,1]^2,
                                 measured=meta(spectra_test)$chl_area)

car_area_val_pred<-data.frame(bulk=meta(spectra_test)$bulk_id,
                              val_pred=predict(car_area_model,newdata=as.matrix(spectra_test),ncomp=ncomp_car_area)[,,1]^2,
                              measured=meta(spectra_test)$car_area)

####################################
## apply to the non-bulk samples

spectra_tagged<-readRDS("processed_data/spectra_tagged.rds")

## change to English for plots
levels(meta(spectra_tagged)$site_id)[levels(meta(spectra_tagged)$site_id)=="MSB_Tourbiere"]<-"MSB_Bog"

meta(spectra_tagged)$LMA<-predict(LMA_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_LMA)[,,1]
meta(spectra_tagged)$LDMC<-predict(LDMC_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_LDMC)[,,1]
meta(spectra_tagged)$EWT<-predict(EWT_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_EWT)[,,1]
meta(spectra_tagged)$chl_area<-predict(chl_area_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_EWT)[,,1]^2
meta(spectra_tagged)$car_area<-predict(car_area_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_EWT)[,,1]^2

## an anthocyanin index, until we have actual anthocyanin data
meta(spectra_tagged)$ARI<-(1/spectra_tagged[,550]-1/spectra_tagged[,700])/100
meta(spectra_tagged)$ARI2<-(1/rowMeans(as.matrix(spectra_tagged[,540:560]))-1/rowMeans(as.matrix(spectra_tagged[,690:710])))/100

meta(spectra_tagged)$mARI<-(1/spectra_tagged[,550]-1/spectra_tagged[,700])*spectra_tagged[,780]
meta(spectra_tagged)$mARI2<-(1/rowMeans(as.matrix(spectra_tagged[,540:560]))-1/rowMeans(as.matrix(spectra_tagged[,690:710])))*rowMeans(as.matrix(spectra_tagged[,760:800]))

## based on ARI regression model from Gitelseon et al. 2001
meta(spectra_tagged)$anth_est<- (meta(spectra_tagged)$ARI-0.0047)/0.0038
## based on ARI regression model from Gitelson et al. 2009
## not suitable for samples with very high ARI
meta(spectra_tagged)$anth_est2<- -log(1-meta(spectra_tagged)$ARI2/0.22)/0.029

colors<-c("Acer rubrum Linnaeus"="firebrick1",
          "Betula populifolia Marshall"="gold1")

ggplot(meta(spectra_tagged),aes(x=date,y=LMA,color=species_id,fill=site_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  coord_cartesian(ylim=c(0,70))+
  scale_color_manual(values=colors)+
  theme_bw()+theme(text=element_text(size=20),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = c(0.25,0.7),
                   legend.key.width = unit(1.5,"cm"),
                   legend.background = element_rect(fill="transparent"))+
  labs(x="Date",color="Species",shape="Site",linetype="Site",
       y="LMA")

ggplot(meta(spectra_tagged),aes(x=date,y=LDMC,color=species_id,fill=site_id))+
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
       y="LDMC")

ggplot(meta(spectra_tagged),aes(x=date,y=EWT,color=species_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  coord_cartesian(ylim=c(0,0.009))+
  scale_color_manual(values=colors)+
  theme_bw()+theme(text=element_text(size=20),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = c(0.25,0.7),
                   legend.key.width = unit(1.5,"cm"),
                   legend.background = element_rect(fill="transparent"))+
  labs(x="Date",color="Species",shape="Site",linetype="Site",
       y="EWT")

ggplot(meta(spectra_tagged),aes(x=date,y=chl_area,color=species_id))+
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
       y="Chl per area")

ggplot(meta(spectra_tagged),aes(x=date,y=car_area,color=species_id))+
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
       y="Car per area")

ggplot(meta(spectra_tagged),aes(x=date,y=ARI2,color=species_id))+
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
       y="ARI")

ggplot(meta(spectra_tagged),aes(x=date,y=anth_est,color=species_id))+
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
       y="Estimated anthocyanins")

###########################################
## visualize tagged spectra

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