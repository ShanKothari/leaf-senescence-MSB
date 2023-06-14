library(spectrolab)
library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)

setwd("C:/Users/querc/Dropbox/PostdocProjects/LeafSenescence")

#####################################
## read spectra

spectra_dirs<-list.dirs(path="spectra")
spectra_list<-lapply(spectra_dirs[-1],read_spectra,format="sig")
spectra_all<-Reduce(spectrolab::combine,spectra_list)
names(spectra_all)<-gsub(pattern=".sig",replacement="",names(spectra_all))

####################################
## attach metadata

spectra_md<-read.csv("spectra/2021-LeafSenescence_spectra.csv")
spectra_md$file_number<-str_pad(spectra_md$file_number,4,pad="0")
spectra_md$file_name<-paste(spectra_md$base_file_name,spectra_md$file_number,sep="_")

spectra_match<-match(names(spectra_all),spectra_md$file_name)
meta(spectra_all)$plant_id<-spectra_md$plant_id[spectra_match]
meta(spectra_all)$bulk_id<-spectra_md$sample_id[spectra_match]
meta(spectra_all)$leaf_number<-spectra_md$leaf_number[spectra_match]
meta(spectra_all)$site_id<-spectra_md$site_id[spectra_match]
meta(spectra_all)$date<-as.Date(spectra_md$date_measured[spectra_match],"%m/%d/%Y")
meta(spectra_all)$quality<-spectra_md$quality[spectra_match]
meta(spectra_all)$sample_remarks<-spectra_md$sample_remarks[spectra_match]
meta(spectra_all)$week<-sapply(meta(spectra_all)$date,
                               function(x){
                                 if(is.na(x)){return(NA)}
                                 else if(x<"2021-08-22"){return("1")}
                                 else if(x<"2021-08-29"){return("2")}
                                 else if(x<"2021-09-12"){return("4")}
                                 else if(x<"2021-09-26"){return("6")}
                                 else if(x<"2021-10-03"){return("7")}
                                 else if(x<"2021-10-10"){return("8")}
                                 else if(x<"2021-10-17"){return("9")}
                                 else if(x<"2021-10-24"){return("10")}
                                 else if(x<"2021-10-31"){return("11")}
                                 else if(x<"2021-11-07"){return("12")}
                                 else if(x<"2021-11-14"){return("13")}
                                 else if(x<"2021-11-21"){return("14")}
                               })
meta(spectra_all)$week<-factor(meta(spectra_all)$week,levels=c(1,2,4,6,7,8,9,10,11,12,13,14))

plants<-read.csv("summary/plants.csv")
meta(spectra_all)$species_id<-plants$scientific_name[match(meta(spectra_all)$plant_id,plants$plant_id)]

#################################
## filter out bad/irrelevant spectra

spectra_all<-spectra_all[-which(toupper(meta(spectra_all)$quality)=="BAD"),]
## Spectralon scans
spectra_all<-spectra_all[-which(meta(spectra_all)$leaf_number=="white target"),]
spectra_all<-spectra_all[-which(meta(spectra_all)$leaf_number=="white reference"),]
## black references (perhaps correct for these later)
spectra_all<-spectra_all[-which(meta(spectra_all)$leaf_number=="black target"),]
## accidents and other 'irrelevant' spectra
spectra_all<-spectra_all[-which(meta(spectra_all)$leaf_number=="-")]
## missing a leaf number?
spectra_all<-spectra_all[-which(is.na(meta(spectra_all)$leaf_number))]
## senesced leaves scanned in the field
spectra_all<-spectra_all[-grep("senesced",meta(spectra_all)$sample_remarks)]

####################################
## matching sensors

sensor.ends<-NULL
## find indices of last wavelength for first two sensors
for(i in 1:(ncol(spectra_all)-1)){if(bands(spectra_all)[i]>bands(spectra_all)[i+1]) sensor.ends<-c(sensor.ends,i)}

spectra_s1_bands<-bands(spectra_all)[1:sensor.ends[1]]
spectra_s2_bands<-bands(spectra_all)[(sensor.ends[1]+1):sensor.ends[2]]
spectra_s3_bands<-bands(spectra_all)[(sensor.ends[2]+1):ncol(spectra_all)]

spectra_s1<-spectra_all[,spectra_s1_bands]
spectra_s2<-spectra_all[,spectra_s2_bands]
spectra_s3<-spectra_all[,spectra_s3_bands]
  
spectra_s1_resamp<-as.matrix(resample(spectra_s1,new_bands = 400:1010))
spectra_s2_resamp<-as.matrix(resample(spectra_s2,new_bands = 1000:1910))
spectra_s3_resamp<-as.matrix(resample(spectra_s3,new_bands = 1895:2400))
spectra_resamp<-do.call(cbind,args = list(spectra_s1_resamp,
                                          spectra_s2_resamp,
                                          spectra_s3_resamp))
# rownames(spectra_resamp)[1]<-"sample_id"
spectra_resamp_long<-melt(spectra_resamp,id.vars="sample_id")
colnames(spectra_resamp_long)<-c("sample_id","wavelength","reflectance")

spectra_resamp_long$wvl_id<-paste(spectra_resamp_long$sample_id,spectra_resamp_long$wavelength,sep="_")
dup_ids_ref<-spectra_resamp_long$wvl_id[duplicated(spectra_resamp_long$wvl_id)]
spectra_resamp_long_no_dups<-spectra_resamp_long[-which(spectra_resamp_long$wvl_id %in% dup_ids_ref),]
inter_wvls<-400:2400

## this function is for linear interpolation over the sensor overlap region
interpolate <- function(x) {
  wvls <- x$wavelength
  ref <- x$reflectance
  new_refs <- approx(wvls, ref, xout = inter_wvls)$y
  tmp <- data_frame(wavelength = inter_wvls, reflectance = new_refs)
  return(tmp)
}

## apply linear interpolation step over sensor overlap
spectra_resamp_long_cleaned <-spectra_resamp_long_no_dups%>%
  group_by(sample_id) %>%
  do(interpolate(.))

spectra_cleaned<-reshape(data.frame(spectra_resamp_long_cleaned),
                         idvar="sample_id",
                         timevar="wavelength",
                         direction="wide")

spectra_cleaned<-spectra(spectra_cleaned[,-1],
                         bands=400:2400,
                         names=spectra_cleaned[,1])
spectra_matched<-match_sensors(spectra_cleaned,splice_at = 1005,interpolate_wvl = 10)

## check that data are in the same order
sum(names(spectra_matched)!=names(spectra_all))
## and if so (returns 0), reattach metadata
meta(spectra_matched)<-meta(spectra_all)

###########################################
## split bulk samples from others

spectra_tagged<-spectra_matched[meta(spectra_matched)$leaf_number!="bulk"]
spectra_bulk<-spectra_matched[meta(spectra_matched)$leaf_number=="bulk"]
saveRDS(spectra_tagged,"processed_data/spectra_tagged.rds")
saveRDS(spectra_bulk,"processed_data/spectra_bulk.rds")

#######################################
## senesced, dried spectra

## these are spectra of senesced leaves measured in the lab
sen_spectra_raw<-read_spectra("spectra/2021-11-29_senesced",format="sig")

sensor.ends<-NULL
## find indices of last wavelength for first two sensors
for(i in 1:(ncol(sen_spectra_raw)-1)){if(bands(sen_spectra_raw)[i]>bands(sen_spectra_raw)[i+1]) sensor.ends<-c(sensor.ends,i)}

sen_s1_bands<-bands(sen_spectra_raw)[1:sensor.ends[1]]
sen_s2_bands<-bands(sen_spectra_raw)[(sensor.ends[1]+1):sensor.ends[2]]
sen_s3_bands<-bands(sen_spectra_raw)[(sensor.ends[2]+1):ncol(sen_spectra_raw)]

sen_s1<-sen_spectra_raw[,sen_s1_bands]
sen_s2<-sen_spectra_raw[,sen_s2_bands]
sen_s3<-sen_spectra_raw[,sen_s3_bands]

sen_s1_resamp<-as.matrix(resample(sen_s1,new_bands = 400:1010))
sen_s2_resamp<-as.matrix(resample(sen_s2,new_bands = 1000:1910))
sen_s3_resamp<-as.matrix(resample(sen_s3,new_bands = 1895:2400))
sen_resamp<-do.call(cbind,args = list(sen_s1_resamp,
                                          sen_s2_resamp,
                                          sen_s3_resamp))

sen_resamp_long<-melt(sen_resamp,id.vars="sample_id")
colnames(sen_resamp_long)<-c("sample_id","wavelength","reflectance")

sen_resamp_long$wvl_id<-paste(sen_resamp_long$sample_id,sen_resamp_long$wavelength,sep="_")
dup_ids_ref<-sen_resamp_long$wvl_id[duplicated(sen_resamp_long$wvl_id)]
sen_resamp_long_no_dups<-sen_resamp_long[-which(sen_resamp_long$wvl_id %in% dup_ids_ref),]
inter_wvls<-400:2400

## apply linear interpolation step over sensor overlap
sen_resamp_long_cleaned <-sen_resamp_long_no_dups%>%
  group_by(sample_id) %>%
  do(interpolate(.))

sen_cleaned<-reshape(data.frame(sen_resamp_long_cleaned),
                         idvar="sample_id",
                         timevar="wavelength",
                         direction="wide")

sen_cleaned<-spectra(sen_cleaned[,-1],
                         bands=400:2400,
                         names=sen_cleaned[,1])
sen_spectra<-match_sensors(sen_cleaned,splice_at = 1005,interpolate_wvl = 10)

######################################
## attach metadata to senesced leaves

## read metadata
sen_md<-read.csv("spectra/MSB_Senesced_11-29-2021.csv")
sen_md$Scan<-str_pad(sen_md$Scan,4,pad="0")
sen_md$file_name<-paste(sen_md$Basefile,sen_md$Scan,sep="_")

names(sen_spectra)<-gsub(pattern=".sig",replacement="",names(sen_spectra))
sen_match<-match(names(sen_spectra),sen_md$file_name)

meta(sen_spectra)$plant_id<-sen_md$Plant[sen_match]
meta(sen_spectra)$leaf_number<-sen_md$Leaf[sen_match]
meta(sen_spectra)$Remarks<-sen_md$Remarks[sen_match]
meta(sen_spectra)$leaf_id<-paste(meta(sen_spectra)$plant_id,meta(sen_spectra)$leaf_number,sep="_")

plants<-read.csv("summary/plants.csv")
meta(sen_spectra)$species_id<-plants$scientific_name[match(meta(sen_spectra)$plant_id,plants$plant_id)]
meta(sen_spectra)$site_id<-plants$site_id[match(meta(sen_spectra)$plant_id,plants$plant_id)]
## change to English for plots
meta(sen_spectra)$site_id[which(meta(sen_spectra)$site_id=="MSB_Tourbiere")]<-"MSB_Bog"

sen_spectra<-sen_spectra[-which(meta(sen_spectra)$leaf_number %in% c("WR","BR")),]
sen_spectra<-sen_spectra[-which(meta(sen_spectra)$leaf_id=="NA_NA"),]
sen_spectra<-sen_spectra[-grep("BAD",meta(sen_spectra)$Remarks),]

meta(sen_spectra)$ARI<-(1/(sen_spectra[,550])-1/sen_spectra[,700])/100
meta(sen_spectra)$PSRI<-(sen_spectra[,678]-sen_spectra[,500])/sen_spectra[,750]

colors<-c("Acer rubrum Linnaeus"="firebrick1",
          "Betula populifolia Marshall"="gold1")

ARI_senesced<-ggplot(meta(sen_spectra),
                     aes(x=site_id,y=ARI,color=species_id))+
  geom_violin(size=2)+
  theme_bw()+
  scale_color_manual(values=colors)+
  theme(text=element_text(size=20),
        legend.position = c(0.7,0.8),
        legend.background = element_rect(fill="transparent"))+
  labs(x="Site",color="Species",
       y="Anthocyanin Reflectance Index")+
  guides(color="none")

ggplot(meta(sen_spectra),aes(x=species_id,y=PSRI,fill=site_id))+
  geom_boxplot()+theme_bw()+theme(text=element_text(size=15))

pdf("manuscript/ARI_senesced.pdf",width=5,height=7)
ARI_senesced
dev.off()
