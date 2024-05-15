### This script tests for differences in subcortical volumes before and after
### a stay in Antarctica, using longitudinal freesurfer for site effects

### Kosha Ruparel

library('ggplot2') # v 3.3.2
library('ggseg') # v 1.1.5... 1.6.00
library('ggpubr') # v 0.4.0
library('dplyr') # v 1.0.2
library('lme4') # v 1.1-23... 1.1-26
library('pbkrtest') # v 0.4-8.6... 0.5-0.1
library('R.utils') # v 2.10.1
library('sjPlot') # v 2.8.4... 2.8.6
library('kableExtra') # v 1.3.1... 1.3.4
library('ggrepel') # v 0.8.2
library('tidyverse') # v 1.3.0
library('cowplot')
library('Rmisc')
#library('xtable')
#library('longCombat') # v 0.0.0.90000
#library('tableHTML')

set.seed(20)


################################### Data ###################################

# Load data

datacs <- read.csv('cross_sectional_output/concordia_crew_control_demo_freesurfer7.4_crosssectional.csv')

datalong <- read.csv('longitudinal_output/concordia_crew_control_demo_freesurfer7.4_longitudinal.csv')

# Subset volumes needed for analysis

subcortical <- c("Accumbens.area", "Amygdala",
  "Caudate", "Hippocampus", "Pallidum",
  "Putamen", "Thalamus")

cortical <- c("frontal_vol", "parietal_vol", "temporal_vol", "occipital_vol")

vol_cs <- datacs[, c("subject", "Time", "scanner", "group","prior_data_2019", "PatientAgeYears",subcortical, cortical)]

vol_long <- datalong[, c("subject", "Time", "scanner", "group","prior_data_2019", "PatientAgeYears",subcortical, cortical)]


# Recode time for phantoms to t0
#ind_data <- final_df[final_df$group %in% c("Crew", "Phantom"), ]
#final_df$Time <- recode(final_df$Time, "t1"="t0", "t2"="t0", "t3"="t0")

vol_cs$Time <- as.factor(vol_cs$Time)
vol_long$Time <- as.factor(vol_long$Time)

#long cambat the cs data
mod1 <- longCombat(idvar="subject", batchvar="scanner",
 features=c(subcortical, cortical), timevar="Time", formula="Time",
 ranef="(1|subject)", data=vol_cs[vol_cs$group=="Crew",])

vol_combat <- mod1$data_combat

#vol_combat <- vol_combat[, paste0(c(subcortical, cortical), ".combat")]

#all_data_crew <- cbind(vol_long[vol_cs$group=="Crew",], data_combat)

#row.names(all_data_crew) <- 1:nrow(all_data_crew)

###all_data matches BOX file Ellyn uploaded

#vol_combat
for (region in grep(".combat", names(vol_combat), value=TRUE)) {
  vol_combat[,paste0("perbase_", region)] <- NA
  for (sub in unique(vol_combat$subject)) {
    if ("t0" %in% vol_combat[vol_combat$subject == sub, "Time"]) {
      baseval <- vol_combat[vol_combat$subject == sub & vol_combat$Time == "t0", region]
      timepoints <- unique(vol_combat[vol_combat$subject == sub, "Time"])
      for (tp in timepoints) {
        vol_combat[vol_combat$subject == sub & vol_combat$Time == tp, paste0("perbase_", region)] <- vol_combat[vol_combat$subject == sub & vol_combat$Time == tp, region]/baseval
      }
    }
  }
}
#calculate percent change from baseline- crew
#vol_long
vol_long_crew<-vol_long[vol_long$group=="Crew",]
vol_long_control<-vol_long[vol_long$group=="Control",]

for (region in names(vol_long_crew)[7:17]) {
  vol_long_crew[,paste0("perbase_", region)] <- NA
  for (sub in unique(vol_long_crew$subject)) {
    if ("t0" %in% vol_long_crew[vol_long_crew$subject == sub, "Time"]) {
      baseval <- vol_long_crew[vol_long_crew$subject == sub & vol_long_crew$Time == "t0", region]
      timepoints <- unique(vol_long_crew[vol_long_crew$subject == sub, "Time"])
      for (tp in timepoints) {
        vol_long_crew[vol_long_crew$subject == sub & vol_long_crew$Time == tp, paste0("perbase_", region)] <- vol_long_crew[vol_long_crew$subject == sub & vol_long_crew$Time == tp, region]/baseval
      }
    }
  }
}

#calculate percent change from baseline- controls
for (region in names(vol_long_control)[7:17]) {
  vol_long_control[,paste0("perbase_", region)] <- NA
  for (sub in unique(vol_long_control$subject)) {
    if ("t0" %in% vol_long_control[vol_long_control$subject == sub, "Time"]) {
      baseval <- vol_long_control[vol_long_control$subject == sub & vol_long_control$Time == "t0", region]
      timepoints <- unique(vol_long_control[vol_long_control$subject == sub, "Time"])
      for (tp in timepoints) {
        vol_long_control[vol_long_control$subject == sub & vol_long_control$Time == tp, paste0("perbase_", region)] <- vol_long_control[vol_long_control$subject == sub & vol_long_control$Time == tp, region]/baseval
      }
    }
  }
}

#####
tmp_long<-rbind(vol_long_crew,vol_long_control)

vol_long2 <- reshape2::melt(tmp_long, c("subject", "Time", "scanner","group"),c(subcortical, cortical,paste0("perbase_", c(subcortical, cortical))))
names(vol_long2) <- c("Subject", "Time", "Scanner", "group","Region", "PercentBase")
vol_long2$PercentBase <- vol_long2$PercentBase*100

#combat - not used
vol_combat2 <- reshape2::melt(vol_combat, c("subject", "Time", "scanner"),c(paste0(c(subcortical, cortical),".combat"),paste0("perbase_", c(subcortical, cortical),".combat")))
names(vol_combat2) <- c("Subject", "Time", "Scanner", "Region", "PercentBase")
vol_combat2$PercentBase <- vol_combat2$PercentBase*100


#limit to percent base
vol_long2<-vol_long2[1409:2816,]

vol_long2$Region <- recode(vol_long2$Region,
    "perbase_Accumbens.area"="Accumbens",
    "perbase_Amygdala"="Amygdala",
    "perbase_Caudate"="Caudate",
    "perbase_Hippocampus"="Hippocampus",
    "perbase_Pallidum"="Pallidum",
    "perbase_Putamen"="Putamen",
    "perbase_Thalamus"="Thalamus",
    "perbase_frontal_vol"="Frontal",
    "perbase_parietal_vol"="Parietal",
    "perbase_temporal_vol"="Temporal",
    "perbase_occipital_vol"="Occipital")
vol_long2$Region<-factor(vol_long2$Region)
vol_long2$group<-factor(vol_long2$group)
vol_long2$Time<-factor(vol_long2$Time)
levels(vol_long2$Time)<-c("Pre","Post1","Post2")

vol_combat2$Time<-factor(vol_combat2$Time)
levels(vol_combat2$Time)<-c("Pre","Post1","Post2")


#smry
vol_long_smry <- summarySE(vol_long2, measurevar="PercentBase", groupvars=c("Region","Time","group"),na.rm=T)

vol_combat_smry <- summarySE(vol_combat2[738:1474,], measurevar="PercentBase", groupvars=c("Region","Time"),na.rm=T)

## bar plot ---Figure3
#sub-cortical and cortical- freesurfer long
vol_long_smry$cortical<-0
vol_long_smry$cortical[c(1:42)]<-1
vol_long_smry$cortical[vol_long_smry$cortical==1]<-"Subcortical"
vol_long_smry$cortical[vol_long_smry$cortical==0]<-"Cortical"

p<-ggplot(vol_long_smry[vol_long_smry$group=="Crew",], aes(x=Region, y=PercentBase, fill=Time)) +
    geom_bar(position=position_dodge(.9), colour="black", stat="identity",width=0.9) +
    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=PercentBase-ci, ymax=PercentBase+ci)) +
    coord_cartesian(ylim=c(60,110)) + labs(y= "%Baseline")+
    theme_bw() + scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")+
       theme(legend.position = "top")+  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))

p+facet_wrap( ~ cortical,ncol=1,scales = "free_x")


#controls
pc<-ggplot(vol_long_smry[vol_long_smry$group=="Control",], aes(x=Region, y=PercentBase, fill=Time)) +
    geom_bar(position=position_dodge(.9), colour="black", stat="identity",width=0.9) +
    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=PercentBase-ci, ymax=PercentBase+ci)) +
    coord_cartesian(ylim=c(60,110)) + labs(y= "%Baseline")+
    theme_bw() + scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")+
       theme(legend.position = "top")+  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))

pc+facet_wrap( ~ cortical,ncol=1,scales = "free_x")


## bar plot ---Figure3
#sub-cortical and cortical- combat values
vol_combat_smry$cortical<-0
vol_combat_smry$cortical[c(1:21)]<-1
vol_combat_smry$cortical[vol_combat_smry$cortical==1]<-"Subcortical"
vol_combat_smry$cortical[vol_combat_smry$cortical==0]<-"Cortical"

p<-ggplot(vol_combat_smry, aes(x=Region, y=PercentBase, fill=Time)) +
    geom_bar(position=position_dodge(.9), colour="black", stat="identity",width=0.9) +
    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=PercentBase-ci, ymax=PercentBase+ci)) +
    coord_cartesian(ylim=c(60,110)) + labs(y= "%Baseline")+
    theme_bw() + scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")+
       theme(legend.position = "top")+  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))

p+facet_wrap( ~ cortical,ncol=1,scales = "free_x")

#violin plot-crew
tmp<-vol_long2[vol_long2$group=="Crew",]
row.names(tmp) <- 1:nrow(tmp)
tmp$cortical<-0
tmp$cortical[c(1:469)]<-1

tmp3<-tmp[tmp$cortical==0,]

tmp4<-tmp[tmp$cortical==1,]

pdf(file = "Violin_plots_crew.pdf",width=12, height=8)

ggplot(aes(y = PercentBase,
       x = `Time`,
       fill = Time,color=Time),
       data = tmp3) +
geom_violin(position=position_dodge(),alpha=0.5) +
  geom_line(mapping = aes(group = Subject),
            position = position_dodge(0.1),
            alpha = 0.3) +
  geom_point(mapping = aes(fill = Time, group = Subject),
             size = 1.5, shape = 21,
             position = position_dodge(0.1)) +
labs(title="Crew: Cortical Regions",
     x="Region",
     y="Percent Base") +
     theme(axis.text.x = element_text(size=14,hjust = 1,angle=45),axis.text.y = element_text(size = 14,hjust = 1, vjust=.4),axis.title=element_text(size=15,face="bold"),strip.text = element_text(
       size = 15),legend.text=element_text(size=14)) +facet_grid(~`Region`)


ggplot(aes(y = PercentBase,
       x = `Time`,
       fill = Time,color=Time),
       data = tmp4) +
geom_violin(position=position_dodge(),alpha=0.5) +
  geom_line(mapping = aes(group = Subject),
            position = position_dodge(0.1),
            alpha = 0.3) +
  geom_point(mapping = aes(fill = Time, group = Subject),
             size = 1.5, shape = 21,
             position = position_dodge(0.1)) +
labs(title="Crew: Subcortical Regions",
     x="Region",
     y="Percent Base") +
     theme(axis.text.x = element_text(size=14,hjust = 1,angle=45),axis.text.y = element_text(size = 14,hjust = 1, vjust=.4),axis.title=element_text(size=15,face="bold"),strip.text = element_text(
       size = 15),legend.text=element_text(size=14)) +
facet_grid(~`Region`)
dev.off()


#violin plot-control
tmpc<-vol_long2[vol_long2$group=="Control",]
row.names(tmpc) <- 1:nrow(tmpc)
tmpc$cortical<-0
tmpc$cortical[c(1:427)]<-1

tmp1<-tmpc[tmpc$cortical==0,]

tmp2<-tmpc[tmpc$cortical==1,]

pdf(file = "Violin_plots_controls.pdf",width=12, height=8)

ggplot(aes(y = PercentBase,
       x = `Time`,
       fill = Time,color=Time),
       data = tmp1) +
geom_violin(position=position_dodge(),alpha=0.5) +
  geom_line(mapping = aes(group = Subject),
            position = position_dodge(0.1),
            alpha = 0.3) +
  geom_point(mapping = aes(fill = Time, group = Subject),
             size = 1.5, shape = 21,
             position = position_dodge(0.1)) +
labs(title="Controls: Cortical Regions",
     x="Region",
     y="Percent Base") +
     theme(axis.text.x = element_text(size=14,hjust = 1,angle=45),axis.text.y = element_text(size = 14,hjust = 1, vjust=.4),axis.title=element_text(size=15,face="bold"),strip.text = element_text(
       size = 15),legend.text=element_text(size=14)) +facet_grid(~`Region`)


ggplot(aes(y = PercentBase,
       x = `Time`,
       fill = Time,color=Time),
       data = tmp2) +
geom_violin(position=position_dodge(),alpha=0.5) +
  geom_line(mapping = aes(group = Subject),
            position = position_dodge(0.1),
            alpha = 0.3) +
  geom_point(mapping = aes(fill = Time, group = Subject),
             size = 1.5, shape = 21,
             position = position_dodge(0.1)) +
labs(title="Controls: Subcortical Regions",
     x="Region",
     y="Percent Base") +
     theme(axis.text.x = element_text(size=14,hjust = 1,angle=45),axis.text.y = element_text(size = 14,hjust = 1, vjust=.4),axis.title=element_text(size=15,face="bold"),strip.text = element_text(
       size = 15),legend.text=element_text(size=14))+
facet_grid(~`Region`)
dev.off()

##LMER models

library(lme4)
library(pbkrtest)
library(dplyr)
library(lmerTest)
library(effects)

#models on crew only
rois=c("Accumbens.area", "Amygdala", "Caudate", "Hippocampus",
  "Pallidum", "Putamen", "Thalamus", "frontal_vol", "parietal_vol", "temporal_vol", "occipital_vol")


data<-datalong[datalong$group=="Crew",]
data<-data[,c("subject", "Time", "scanner","group","prior_data_2019","PatientAgeYears","eTIV",rois)]
#data<-vol_long_crew[,1:17]
data$age<-floor(data$PatientAgeYears)
data$Time<-factor(data$Time)
data_long<-reshape2::melt(data, c("subject", "Time", "scanner","group","prior_data_2019","PatientAgeYears","eTIV"),c(subcortical, cortical))
names(data_long) <- c("subject", "time", "scanner", "group","prior_data_2019","patientageyears","eTIV","region", "vol")


results <- data.frame(Region=rois, Pval_agetimeint=rep(NA, length(rois)),
    Pval_etivtimeint=rep(NA, length(rois)),Pval_finalmodelTime=rep(NA, length(rois)))

results_table_t12_t0 <- data.frame(Region=rois, Estimate=rep(NA, length(rois)),StdErr=rep(NA, length(rois)),
                    Pval=rep(NA, length(rois)))

results_table_t18_t0 <- data.frame(Region=rois, Estimate=rep(NA, length(rois)),StdErr=rep(NA, length(rois)),
                Pval=rep(NA, length(rois)))

results_table_t12_t0_scaled <- data.frame(Region=rois, Estimate=rep(NA, length(rois)),StdErr=rep(NA, length(rois)),
                                    Pval=rep(NA, length(rois)))

results_table_t18_t0_scaled <- data.frame(Region=rois, Estimate=rep(NA, length(rois)),StdErr=rep(NA, length(rois)),
                                Pval=rep(NA, length(rois)))


                                ## running lmer without outlier
                                #sphaghetti plot
                                ggplot(data, aes(Time,Accumbens.area )) +
                                  geom_line(aes(group = factor(subject))) +
                                  geom_smooth()+
                                  facet_wrap(~ subject,ncol=3,scales="free")

                                ggplot(data_long, aes(time, vol)) +
                                  geom_line(aes(group = factor(subject))) +
                                  geom_smooth()+
                                  facet_grid(~ region)
  ## so removing 007 who has a 140% increase at t3
  data_scaled_no007<-data_scaled[-which(data_scaled$subject=="concordia_007"),]

  results_table_t12_t0_scaled_no007 <- data.frame(Region=rois, Estimate=rep(NA, length(rois)),StdErr=rep(NA, length(rois)),
                                      Pval=rep(NA, length(rois)))

  results_table_t18_t0_scaled_no007 <- data.frame(Region=rois, Estimate=rep(NA, length(rois)),StdErr=rep(NA, length(rois)),
                                  Pval=rep(NA, length(rois)))
pdf("Crew_time_effects.pdf")

#lmer throws a warning of covariates being on different scale
data_scaled<-data
data_scaled[c(7:18)] <- lapply(data_scaled[c(7:18)], function(x) c(scale(x)))

for (r in rois) {
  m1 <- lmer(formula(paste(r, "~ Time + age +eTIV+(1|subject) ")), data=data)
  m2 <- lmer(formula(paste(r, "~ Time * age +eTIV+(1|subject) ")), data=data)
  anova1<-anova(m1,m2)
  m3 <- lmer(formula(paste(r, "~ Time +eTIV+age+(1|subject) ")), data=data)
  m4 <- lmer(formula(paste(r, "~ Time*eTIV+age+(1|subject) ")), data=data)
  anova2<-anova(m3,m4)
  #above pvalues are not significant and twobarely significant will not survive fdr
  m5 <- lmer(formula(paste(r, "~ Time+eTIV+age+(1|subject) ")), data=data)
  m5_sc <- lmer(formula(paste(r, "~ Time+eTIV+age+(1|subject) ")), data=data_scaled)
  m5_sc_no007 <- lmer(formula(paste(r, "~ Time+eTIV+age+(1|subject) ")), data=data_scaled_no007)

  #estimate
  results_table_t12_t0[results$Region == r, "Estimate"]<-summary(m5)$coefficients[2,1]
  results_table_t18_t0[results$Region == r, "Estimate"]<-summary(m5)$coefficients[3,1]

  results_table_t12_t0_scaled[results$Region == r, "Estimate"]<-summary(m5_sc)$coefficients[2,1]
  results_table_t18_t0_scaled[results$Region == r, "Estimate"]<-summary(m5_sc)$coefficients[3,1]

  results_table_t12_t0_scaled_no007[results$Region == r, "Estimate"]<-summary(m5_sc_no007)$coefficients[2,1]
  results_table_t18_t0_scaled_no007[results$Region == r, "Estimate"]<-summary(m5_sc_no007)$coefficients[3,1]
  #std err
  results_table_t12_t0[results$Region == r, "StdErr"]<-summary(m5)$coefficients[2,2]
  results_table_t18_t0[results$Region == r, "StdErr"]<-summary(m5)$coefficients[3,2]
  results_table_t12_t0_scaled[results$Region == r, "StdErr"]<-summary(m5_sc)$coefficients[2,2]
  results_table_t18_t0_scaled[results$Region == r, "StdErr"]<-summary(m5_sc)$coefficients[3,2]
  results_table_t12_t0_scaled_no007[results$Region == r, "StdErr"]<-summary(m5_sc_no007)$coefficients[2,2]
  results_table_t18_t0_scaled_no007[results$Region == r, "StdErr"]<-summary(m5_sc_no007)$coefficients[3,2]
  #pval
  results_table_t12_t0[results$Region == r, "Pval"]<-format(summary(m5)$coefficients[2,5],scientific=FALSE)
  results_table_t18_t0[results$Region == r, "Pval"]<-format(summary(m5)$coefficients[3,5],scientific=FALSE)

  results_table_t12_t0_scaled[results$Region == r, "Pval"]<-format(summary(m5_sc)$coefficients[2,5],scientific=FALSE)
  results_table_t18_t0_scaled[results$Region == r, "Pval"]<-format(summary(m5_sc)$coefficients[3,5],scientific=FALSE)

  results_table_t12_t0_scaled_no007[results$Region == r, "Pval"]<-format(summary(m5_sc_no007)$coefficients[2,5],scientific=FALSE)
  results_table_t18_t0_scaled_no007[results$Region == r, "Pval"]<-format(summary(m5_sc_no007)$coefficients[3,5],scientific=FALSE)
  anova3<-anova(m5,ddf="Kenward-Roger")
  e <- allEffects(m5)
  plot(e)
  results[results$Region == r, "Pval_agetimeint"] <- anova1$P[2]
  results[results$Region == r, "Pval_etivtimeint"] <- anova2$P[2]
  results[results$Region == r, "Pval_finalmodelTime"] <- anova3$P[1]
}

dev.off()

results_table_t12_t0$Effect<-"Timet12"
results_table_t18_t0$Effect<-"Timet18"
results_table_out<-rbind(results_table_t12_t0,results_table_t18_t0)
results_table_out <- arrange(results_table_out,Region,Effect)

#reviewer has asked for raw mean and sd
results_table_raw_mean_sd<-as.data.frame(summarise(group_by(data_long[,c("region","time","vol")], region, time),
          mean=mean(vol), sd=sd(vol)))
library(dplyr)
raw_mean_sd_table<-data_long %>%
  group_by(region, time) %>%
  summarise(MEAN= mean(vol), SD = sd(vol), N = length(vol))


  #violin plot-crew minus 007
  tmp<-vol_long2[vol_long2$group=="Crew",]
  tmp<-tmp[-which(tmp$Subject=="concordia_007"),]
  row.names(tmp) <- 1:nrow(tmp)
  tmp$cortical<-0
  tmp$cortical[c(1:448)]<-1

  tmp3<-tmp[tmp$cortical==0,]
  tmp3$Region<-factor(tmp3$Region)
  tmp4<-tmp[tmp$cortical==1,]
  tmp4$Region<-factor(tmp4$Region)

  pdf(file = "Violin_plots_crew_no007.pdf",width=12, height=8)

  ggplot(aes(y = PercentBase,
         x = `Time`,
         fill = Time,color=Time),
         data = tmp3) +
  geom_violin(position=position_dodge(),alpha=0.5) +
    geom_line(mapping = aes(group = Subject),
              position = position_dodge(0.1),
              alpha = 0.3) +
    geom_point(mapping = aes(fill = Time, group = Subject),
               size = 1.5, shape = 21,
               position = position_dodge(0.1)) +
  labs(title="Crew: Cortical Regions",
       x="Region",
       y="Percent Base",size=15) +
  theme(axis.text.x = element_text(size=14,hjust = 1,angle=45),axis.text.y = element_text(size = 14,hjust = 1, vjust=.4),axis.title=element_text(size=15,face="bold"),strip.text = element_text(
    size = 15),legend.text=element_text(size=14)) +facet_grid(~`Region`)



  ggplot(aes(y = PercentBase,
         x = `Time`,
         fill = Time,color=Time),
         data = tmp4) +
  geom_violin(position=position_dodge(),alpha=0.5) +
    geom_line(mapping = aes(group = Subject),
              position = position_dodge(0.1),
              alpha = 0.3) +
    geom_point(mapping = aes(fill = Time, group = Subject),
               size = 1.5, shape = 21,
               position = position_dodge(0.1)) +
  labs(title="Crew: Subcortical Regions",
       x="Region",
       y="Percent Base") +
  theme(axis.text.x = element_text(size=14,hjust = 1,angle=45),axis.text.y = element_text(size = 14,hjust = 1, vjust=.4),axis.title=element_text(size=15,face="bold"),strip.text = element_text(
    size = 15),legend.text=element_text(size=14)) +
  facet_grid(~`Region`)
  dev.off()


#crew versus control-baseline only
vol_long_all <- reshape2::melt(tmp_long, c("subject", "Time", "scanner","group"),c(subcortical, cortical))
names(vol_long_all)[5:6]<-c("region","vol")
vol_long_all_t0<-vol_long_all[vol_long_all$Time=="t0",]

vol_long_all_smry <- summarySE(vol_long_all_t0, measurevar="vol", groupvars=c("region","group"),na.rm=T)
vol_long_all_smry$cortical<-1
vol_long_all_smry$cortical[1:14]<-0
levels(vol_long_all_smry$region)<-c("Accumbens","Amygdala","Caudate","Hippocampus","Pallidum","Putamen","Thalamus","Frontal","Parietal","Temporal","Occipital")

#bar plot

ggplot(vol_long_all_smry[vol_long_all_smry$cortical==0,], aes(x=region, y=vol, fill=group)) +
    geom_bar(position=position_dodge(.9), colour="black", stat="identity",width=0.9) +
    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=vol-ci, ymax=vol+ci)) + labs(y= "Volume")+
    theme_bw() + scale_fill_brewer(palette = "Dark2") +coord_cartesian(ylim=c(0,20000))+
  scale_color_brewer(palette = "Dark2")+
       theme(legend.position = "top")+  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))


        ggplot(vol_long_all_smry[vol_long_all_smry$cortical==1,], aes(x=region, y=vol, fill=group)) +
            geom_bar(position=position_dodge(.9), colour="black", stat="identity",width=0.9) +
            geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=vol-ci, ymax=vol+ci)) + labs(y= "Volume")+
            theme_bw() + scale_fill_brewer(palette = "Dark2") +coord_cartesian(ylim=c(0,200000))+
          scale_color_brewer(palette = "Dark2")+
               theme(legend.position = "top")+  theme(axis.text=element_text(size=12),
                axis.title=element_text(size=12,face="bold"))


### Figure 4


#data combined by Tyler and sent via secure share

cor_data<-read.csv("merged_concordia_neurocog_imaging_tymoore_may2024.csv")

#cor_data$SU_NotStressed_SLOPE<-cor_data$SU_Stressed_SLOPE*(-1)
#cor_data$SU_NotUnhappy_SLOPE<-cor_data$SU_Unhappy_SLOPE*(-1)
#cor_data$SU_NotPhysExhaust_SLOPE<-cor_data$SU_PhysExhaust_SLOPE*(-1)
#cor_data$SU_NotMentFatigued_SLOPE<-cor_data$SU_MentFatigued_SLOPE*(-1)
#cor_data$SU_NotLonely_SLOPE<-cor_data$SU_Lonely_SLOPE*(-1)
#cor_data$SU_NotSleepy_SLOPE<-cor_data$SU_Sleepy_SLOPE*(-1)
#cor_data$SU_GoodSleepQual_SLOPE<-cor_data$SU_PoorSleepQual_SLOPE*(-1)

#tmp<-read.csv("figure4_corr_info_2024.csv")
#keep_vars<-c(tmp$metric, colnames(tmp)[4:14])

#
library(reshape2)

# load library ggplot2
library(ggplot2)
library(ggpubr)

#cor_data2<-cor_data[,c("ID",keep_vars)]
#su<-cor_data2[,c(1,32:44)]
#vol<-cor_data2[,c(1,45:55)]

###sleep only
s_plot1<-ggplot( cor_data, aes( x=TIB_mean , y=CortexVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Time in Bed", y = "Cortical Volume PercentChange")
s_plot11<-ggplot( cor_data, aes( x=TIB_mean , y=SubCortGrayVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Time in Bed", y = "Subcortical Volume PercentChange")



s_plot2<-ggplot( cor_data, aes( x=TotalSleep_mean  , y=CortexVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Total Sleep Time", y = "Cortical Volume PercentChange")
s_plot22<-ggplot( cor_data, aes( x=TotalSleep_mean  , y=SubCortGrayVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Total Sleep Time", y = "Subcortical Volume PercentChange")


s_plot3<-ggplot( cor_data, aes( x=NightSleep_mean, y=CortexVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Night Sleep", y = "Cortical Volume PercentChange")
s_plot33<-ggplot( cor_data, aes( x=NightSleep_mean, y=SubCortGrayVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Night Sleep", y = "Subcortical Volume PercentChange")


s_plot4<-ggplot( cor_data, aes( x=DaySleep_mean , y=CortexVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Day Sleep", y = "Cortical Volume PercentChange")
s_plot44<-ggplot( cor_data, aes( x=DaySleep_mean , y=SubCortGrayVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Day Sleep", y = "Subcortical Volume PercentChange")

s_plot5<-ggplot( cor_data, aes( x=SleepEff_mean , y=CortexVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Sleep Efficiency", y = "Cortical Volume PercentChange")

s_plot55<-ggplot( cor_data, aes( x=SleepEff_mean , y=SubCortGrayVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="Sleep Efficiency", y = "Subcortical Volume PercentChange")


ggarrange(s_plot22, s_plot55 + rremove("x.text"), labels = c("A", "B", "C","D","E","F"),ncol = 2, nrow = 1)



c_plot22<-ggplot( cor_data, aes( x=MRT_pCorr_mean , y=SubCortGrayVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="MRT Accuracy", y = "Subcortical Volume PercentChange")

c_plot1<-ggplot( cor_data, aes( x=MP_Accuracy_mean , y=CortexVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="MP Accuracy", y = "Cortical Volume PercentChange")

c_plot11<-ggplot( cor_data, aes( x=MP_Accuracy_mean , y=SubCortGrayVol.y_percent_change ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson")+labs(x ="MP Accuracy", y = "Subcortical Volume PercentChange")



ggarrange(c_plot22, c_plot11, c_plot1, + rremove("x.text"), labels = c("A", "B", "C"),ncol = 2, nrow = 2)

###code below not used for 2024 resubmission

su_long<-melt(su, id.vars=c("ID"))
names(su_long)[2:3]<-c("Survey","survey_value")

vol_long<-melt(vol, id.vars=c("ID"))
names(vol_long)[2:3]<-c("Vol","vol_value")

#significant survey corr
su_vol<-merge(su_long,vol_long,by="ID")
su_vol_p1<-su_vol[which(su_vol$Survey=="SU_Healthy_SLOPE" & su_vol$Vol=="occipital_vol_Real_Change"),]
su_vol_p2<-su_vol[which(su_vol$Survey=="SU_Healthy_SLOPE" & su_vol$Vol=="frontal_vol_Real_Change"),]
su_vol_p3<-su_vol[which(su_vol$Survey=="SU_GoodSleepQual_SLOPE" & su_vol$Vol=="Pallidum_Real_Change"),]
su_vol_p4<-su_vol[which(su_vol$Survey=="SU_NotMentFatigued_SLOPE" & su_vol$Vol=="Pallidum_Real_Change"),]
su_vol_p5<-su_vol[which(su_vol$Survey=="SU_Fresh_SLOPE" & su_vol$Vol=="Pallidum_Real_Change"),]
su_vol_p6<-su_vol[which(su_vol$Survey=="SU_Fresh_SLOPE" & su_vol$Vol=="Caudate_Real_Change"),]

#plots: survey
su_plot1<-ggplot( su_vol_p1, aes( x=survey_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -1, label.y = 2)+labs(x ="Healthy", y = "Occipital Volume")

su_plot2<-ggplot( su_vol_p2, aes( x=survey_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -1, label.y = 2)+labs(x ="Healthy", y = "Frontal Volume")

su_plot3<-ggplot( su_vol_p3, aes( x=survey_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -1, label.y = 2)+labs(x ="Good Sleep Quality", y = "Pallidum Volume")

su_plot4<-ggplot( su_vol_p4, aes( x=survey_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -1, label.y = 2)+labs(x ="Not Mentally Fauigued", y = "Pallidum Volume")

su_plot5<-ggplot( su_vol_p5, aes( x=survey_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2)+labs(x ="Fresh", y = "Pallidum Volume")


su_plot6<-ggplot( su_vol_p6, aes( x=survey_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2)+labs(x ="Fresh", y = "Caudate Volume")


ggarrange(su_plot1, su_plot2, su_plot3,su_plot4,su_plot5,su_plot6 + rremove("x.text"), labels = c("A", "B", "C","D","E","F"),ncol = 3, nrow = 2)



###Activity only
ac<-cor_data2[,c(1,27:32)]
ac_long<-melt(ac, id.vars=c("ID"))
names(ac_long)[2:3]<-c("Activity","act_value")


#significant survey corr
ac_vol<-merge(ac_long,vol_long,by="ID")
ac_vol_p1<-ac_vol[which(ac_vol$Activity=="Light.Activity_SLOPE" & ac_vol$Vol=="parietal_vol_Real_Change"),]
ac_vol_p2<-ac_vol[which(ac_vol$Activity=="Light.Activity_SLOPE" & ac_vol$Vol=="temporal_vol_Real_Change"),]
ac_vol_p3<-ac_vol[which(ac_vol$Activity=="Light.Activity_SLOPE" & ac_vol$Vol=="Amygdala_Real_Change"),]
ac_vol_p4<-ac_vol[which(ac_vol$Activity=="Light.Activity_SLOPE" & ac_vol$Vol=="Caudate_Real_Change"),]
ac_vol_p5<-ac_vol[which(ac_vol$Activity=="Light.Activity_SLOPE" & ac_vol$Vol=="Thalamus_Real_Change"),]
ac_vol_p6<-ac_vol[which(ac_vol$Activity=="Moderate.Activity_SLOPE" & ac_vol$Vol=="occipital_vol_Real_Change"),]
ac_vol_p7<-ac_vol[which(ac_vol$Activity=="Moderate.Activity_SLOPE" & ac_vol$Vol=="Amygdala_Real_Change"),]
ac_vol_p8<-ac_vol[which(ac_vol$Activity=="Vigorous.Activity_SLOPE" & ac_vol$Vol=="temporal_vol_Real_Change"),]
ac_vol_p9<-ac_vol[which(ac_vol$Activity=="Vigorous.Activity_SLOPE" & ac_vol$Vol=="parietal_vol_Real_Change"),]
ac_vol_p10<-ac_vol[which(ac_vol$Activity=="Vigorous.Activity_SLOPE" & ac_vol$Vol=="frontal_vol_Real_Change"),]
ac_vol_p11<-ac_vol[which(ac_vol$Activity=="Vigorous.Activity_SLOPE" & ac_vol$Vol=="Amygdala_Real_Change"),]
ac_vol_p12<-ac_vol[which(ac_vol$Activity=="Vigorous.Activity_SLOPE" & ac_vol$Vol=="Thalamus_Real_Change"),]


#plots: activity
ac_plot1<-ggplot( ac_vol_p1, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Light Activity", y = "Parietal Volume")

ac_plot2<-ggplot( ac_vol_p1, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Light Activity", y = "Temporal Volume")

ac_plot3<-ggplot( ac_vol_p1, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Light Activity", y = "Amygdala Volume")
ac_plot4<-ggplot( ac_vol_p1, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Light Activity", y = "Caudate Volume")
ac_plot5<-ggplot( ac_vol_p1, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Light Activity", y = "Thalamus Volume")

ac_plot6<-ggplot( ac_vol_p1, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Moderate Activity", y = "Occiptital Volume")

ac_plot7<-ggplot( ac_vol_p1, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Moderate Activity", y = "Amygdala Volume")

ac_plot8<-ggplot( ac_vol_p2, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Vigorous Activity", y = "Temporal Volume")

ac_plot9<-ggplot( ac_vol_p2, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Vigorous Activity", y = "Parietal Volume")

ac_plot10<-ggplot( ac_vol_p2, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Vigorous Activity", y = "Frontal Volume")

ac_plot11<-ggplot( ac_vol_p2, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Vigorous Activity", y = "Amygdala Volume")


ac_plot12<-ggplot( ac_vol_p2, aes( x=act_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="Vigorous Activity", y = "Thalamus Volume")

ggarrange(ac_plot1, ac_plot2, ac_plot3,ac_plot4,ac_plot5,ac_plot6,ac_plot7,ac_plot8,ac_plot9,ac_plot10,ac_plot11,ac_plot12+ rremove("x.text"), labels = c("A", "B", "C","D","E","F","G","H","I","J","K","L"),ncol = 4, nrow = 3)


###Cognition only
cog<-cor_data2[,c(1:21)]
cog_long<-melt(cog, id.vars=c("ID"))
names(cog_long)[2:3]<-c("Cognition","cog_value")
#significant survey corr
cog_vol<-merge(cog_long,vol_long,by="ID")
cog_vol_p1<-cog_vol[which(cog_vol$Cognition=="BART_RiskScoreE_SLOPE" & cog_vol$Vol=="occipital_vol_Real_Change"),]
cog_vol_p2<-cog_vol[which(cog_vol$Cognition=="PVT_Slowness_SLOPE" & cog_vol$Vol=="Accumbens_area_Real_Change"),]
cog_vol_p3<-cog_vol[which(cog_vol$Cognition=="MP_AvRT_SLOPE" & cog_vol$Vol=="Accumbens_area_Real_Change"),]
cog_vol_p4<-cog_vol[which(cog_vol$Cognition=="DSST_pCorr_SLOPE" & cog_vol$Vol=="Hippocampus_Real_Change"),]
cog_vol_p5<-cog_vol[which(cog_vol$Cognition=="MP_Accuracy_SLOPE" & cog_vol$Vol=="Putamen_Real_Change"),]


#plots: cortical cognition
cog_plot1<-ggplot( cog_vol_p1, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="BART Accuracy", y = "Occipital Volume")

cog_plot2<-ggplot( cog_vol_p2, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="PVT Speed", y = "Accumbens Volume")

cog_plot3<-ggplot( cog_vol_p3, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 3.5)+labs(x ="MP Speed", y = "Accumbens Volume")

cog_plot4<-ggplot( cog_vol_p4, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2.5, label.y = 2.5)+labs(x ="DSST Accuracy", y = "Hippocampus Volume")

cog_plot5<-ggplot( cog_vol_p5, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2)+labs(x ="MP Accuracy", y = "Putamen Volume")


ggarrange(cog_plot1, cog_plot2, cog_plot3,cog_plot4,cog_plot5+ rremove("x.text"), labels = c("A", "B", "C","D","E"),ncol = 4, nrow = 2)


#code below not touched in 2024

cog_vol_p10<-cog_vol[which(cog_vol$Cognition=="BART_AvRT_SLOPE" & cog_vol$Vol=="combat_vol_miccai_ave_Accumbens_Area_REAL_CHANGE"),]
cog_vol_p11<-cog_vol[which(cog_vol$Cognition=="MRT_AvRT_SLOPE" & cog_vol$Vol=="combat_vol_miccai_ave_Accumbens_Area_REAL_CHANGE"),]
cog_vol_p12<-cog_vol[which(cog_vol$Cognition=="LOT_AvRT_SLOPE" & cog_vol$Vol=="combat_vol_miccai_ave_Accumbens_Area_REAL_CHANGE"),]
cog_vol_p13<-cog_vol[which(cog_vol$Cognition=="VOLT_AvRT_SLOPE" & cog_vol$Vol=="combat_vol_miccai_ave_Accumbens_Area_REAL_CHANGE"),]
cog_vol_p14<-cog_vol[which(cog_vol$Cognition=="ERT_AvRT_SLOPE" & cog_vol$Vol=="combat_vol_miccai_ave_Amygdala_REAL_CHANGE"),]
cog_vol_p15<-cog_vol[which(cog_vol$Cognition=="LOT_AvRT_SLOPE" & cog_vol$Vol=="combat_vol_miccai_ave_Amygdala_REAL_CHANGE"),]
cog_vol_p16<-cog_vol[which(cog_vol$Cognition=="AM_pCorr_SLOPE" & cog_vol$Vol=="combat_vol_miccai_ave_Amygdala_REAL_CHANGE"),]
cog_vol_p17<-cog_vol[which(cog_vol$Cognition=="BART_AvRT_SLOPE" & cog_vol$Vol=="combat_vol_miccai_ave_Putamen_REAL_CHANGE"),]
cog_vol_p18<-cog_vol[which(cog_vol$Cognition=="LOT_AvRT_SLOPE" & cog_vol$Vol=="combat_vol_miccai_ave_Putamen_REAL_CHANGE"),]
#plots: cortical cognition
cog_plot10<-ggplot( cog_vol_p10, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2.5, label.y = 2.5)+labs(x ="BART Speed", y = "Accumbens Volume")

cog_plot11<-ggplot( cog_vol_p11, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="MRT Speed", y = "Accumbens Volume")

cog_plot12<-ggplot( cog_vol_p12, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2.5, label.y = 2.5)+labs(x ="LOT Speed", y = "Accumbens Volume")

cog_plot13<-ggplot( cog_vol_p13, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 3.5)+labs(x ="VOLT Speed", y = "Accumbens Volume")

cog_plot14<-ggplot( cog_vol_p14, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -3.5, label.y = 3)+labs(x ="ERT pCorr", y = "Amygdala Volume")

cog_plot15<-ggplot( cog_vol_p15, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2, label.y = 2.5)+labs(x ="LOT Speed", y = "Amygdala Volume")

cog_plot16<-ggplot( cog_vol_p16, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2.5, label.y = 2.5)+labs(x ="AM pCorr", y = "Amygdala Volume")

cog_plot17<-ggplot( cog_vol_p17, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2.5, label.y = 2.5)+labs(x ="BART Speed", y = "Putamen Volume")

cog_plot18<-ggplot( cog_vol_p18, aes( x=cog_value, y=vol_value ))+ geom_smooth(method = "lm") +geom_point()+
stat_cor(method = "pearson", label.x = -2.5, label.y = 2.5)+labs(x ="LOT Speed", y = "Putamen Volume")

ggarrange(cog_plot10, cog_plot11, cog_plot12,cog_plot13,cog_plot14,cog_plot15,cog_plot16,cog_plot17,cog_plot18+ rremove("x.text"), labels = c("A", "B", "C","D","E","F","G","H","I"),ncol = 4, nrow = 3)
