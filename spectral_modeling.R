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

meta(spectra_tagged)$LMA<-predict(LMA_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_LMA)[,,1]
meta(spectra_tagged)$LDMC<-predict(LDMC_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_LDMC)[,,1]
meta(spectra_tagged)$EWT<-predict(EWT_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_EWT)[,,1]
meta(spectra_tagged)$chl_area<-predict(chl_area_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_EWT)[,,1]^2
meta(spectra_tagged)$car_area<-predict(car_area_model,newdata=as.matrix(spectra_tagged),ncomp=ncomp_EWT)[,,1]^2

ggplot(meta(spectra_tagged),aes(x=date,y=LMA,color=species_id,fill=site_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  coord_cartesian(ylim=c(0,70))+theme_bw()

ggplot(meta(spectra_tagged),aes(x=date,y=LDMC,color=species_id,fill=site_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+theme_bw()

ggplot(meta(spectra_tagged),aes(x=date,y=EWT,color=species_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  theme_bw()+
  coord_cartesian(ylim=c(0,0.009))

ggplot(meta(spectra_tagged),aes(x=date,y=chl_area,color=species_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  theme_bw()

ggplot(meta(spectra_tagged),aes(x=date,y=car_area,color=species_id))+
  stat_summary(geom="point",fun="mean",size=4,aes(shape=site_id))+
  stat_summary(geom="line",fun="mean",size=2,aes(linetype=site_id))+
  stat_summary(geom="errorbar",size=2,fun.data="mean_se")+
  theme_bw()
