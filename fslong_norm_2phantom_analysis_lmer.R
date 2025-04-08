#lmer
library(ggeffects)
library(lme4)
library(emmeans)
library(pbkrtest)
library(ggplot2)
library(visreg)
library(lmerTest)
library(reshape2)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(dplyr)
#getting pvalues out of lme4
#car::Anova and lmerTest::anova provide wrappers for Kenward-Roger-corrected tests using pbkrtest: lmerTest::anova also provides t tests via the Satterthwaite approximation (P,*)
#read data
data<-read.csv("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/Data/fslong_normed_2phantomsbyscanner_crew_control_12242024.csv")
euler<-read.csv("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/Data/concordia_crew_control_euler.csv")
whitematter<-read.csv("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/Data/concordia_normed_whitematter.csv")
data <- data %>%
  left_join(euler %>% select(scanid, euler_total), by = "scanid")
data <- data %>%
  left_join(whitematter %>% select(scanid, CerebralWhiteMatterVol), by = "scanid")


#crew data
data_crew<-data[data$group=="Crew",] #67
#uncomment below for controls
#data_controls<-data[data$group=="Control",]
#regions
rois=c("Accumbens_area", "Amygdala", "Caudate", "Hippocampus",
  "Pallidum", "Putamen", "Thalamus", "frontal_vol", "parietal_vol", "temporal_vol", "occipital_vol","CSF")

data_crew$Time<-as.factor(data_crew$Time)
#data_controls$Time<-as.factor(data_controls$Time)

#uncomment 2 lines below to remove non-CGN time 1 crew
#remove=c("concordia_007", "concordia_010", "concordia_107")
#data_crew=data_crew[!(data_crew$subject %in% remove),]
#above code removes subjects who did not start at CGN

#uncomment 2 lines below to remove controls who did not complete all sessions
#remove1=c("DLR_001", "DLR_002", "DLR_003","DLR_006","DLR_010","DLR_012","DLR_101","DLR_104","DLR_105","DLR_106","DLR_107","DLR_110")
#data_controls=data_controls[!(data_controls$subject %in% remove1),]

#listc <- vector(mode = "list")

#another univariate analyses needs to run with Time and PatientAge- are highly correlated but variability amongst timepoints exists in the crew
pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/lmer_output_crew_fslong_norm2phantom.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/lmer_output_controls_fslong_norm2phantom_COMPLETE_012425.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/lmer_output_crew_fslong_norm2phantom_removeNonCGNTime1.pdf")
#removal of non-CGN Time 1 subjects does not change results

out.lmer1 <- data.frame(Region=rois, F=rep(NA, length(rois)),P=rep(NA, length(rois)))
out.lmer1.posthoc<-as.list(1:length(rois))
names(out.lmer1.posthoc) <- rois

for (roi in rois) {
  f <- formula(paste(roi,"~Time + (1|subject)"))
  lmer.m1<- lmer(f, data=data_crew)
    out.lmer1[out.lmer1$Region == roi, "F"] <- anova(lmer.m1)[[5]]
    out.lmer1[out.lmer1$Region == roi, "P"] <- anova(lmer.m1)[[6]]
    posthoc.lmer1<- emmeans(lmer.m1, pairwise ~ Time, adjust = "bonferroni")
    posthoc.lmer1.contrasts<-as.data.frame(posthoc.lmer1$contrasts)
    posthoc.lmer1.contrasts$roi<-roi
    out.lmer1.posthoc[[roi]]<-posthoc.lmer1.contrasts
   
visreg(lmer.m1,"Time")
}

#redo above and generate effect sizes
library(lme4)
library(emmeans)
library(effectsize)
library(visreg)

# Initialize output data frames
out.lmer1 <- data.frame(Region = rois, F = rep(NA, length(rois)), P = rep(NA, length(rois)))
out.lmer1.posthoc <- as.list(1:length(rois))
names(out.lmer1.posthoc) <- rois
effect_sizes <- data.frame(Region = character(), Contrast = character(), Cohen_d = numeric(), stringsAsFactors = FALSE)

# Loop through each ROI
for (roi in rois) {
  # Define the model formula
  f <- formula(paste(roi, "~ Time + (1|subject)"))
  
  # Fit the linear mixed-effects model
  lmer.m1 <- lmer(f, data = data_crew)
  
  # Store F and P values
  out.lmer1[out.lmer1$Region == roi, "F"] <- anova(lmer.m1)[[5]]
  out.lmer1[out.lmer1$Region == roi, "P"] <- anova(lmer.m1)[[6]]
  
  # Perform post-hoc pairwise comparisons
  posthoc.lmer1 <- emmeans(lmer.m1, pairwise ~ Time, adjust = "bonferroni")
  posthoc.lmer1.contrasts <- as.data.frame(posthoc.lmer1$contrasts)
  posthoc.lmer1.contrasts$roi <- roi
  
  # Save post-hoc results
  out.lmer1.posthoc[[roi]] <- posthoc.lmer1.contrasts
  
  # Calculate Cohen's d for each contrast
  roi_effect_sizes <- eff_size(posthoc.lmer1$emmeans, sigma = sigma(lmer.m1), edf = df.residual(lmer.m1))
  
  # Add ROI and contrast information
  roi_effect_sizes <- as.data.frame(roi_effect_sizes)
  roi_effect_sizes$Region <- roi
  effect_sizes <- rbind(effect_sizes, roi_effect_sizes)
  
  # Visualize the model
  visreg(lmer.m1, "Time")
}

# Display results
print(out.lmer1)
print(effect_sizes)
###
##used https://effect-size-calculator.herokuapp.com/#paired-samples-t-test to produce effect size for the F tests. 


dev.off()

post<-out.lmer1.posthoc

out.lmer1$P_FDR <- round(p.adjust(out.lmer1$P, method='fdr'),4)
out.lmer1$P_Bonferroni <- round(p.adjust(out.lmer1$P, method='bonferroni'),4)

#added on 04/01/2025
#cerebral white matter
source("~/Dropbox/R/summary_se_function.R")
model.cWM <- lmer(CerebralWhiteMatterVol ~ Time + (1 | subject), data = data_crew)
anova(model.cWM)
posthoc.cWM<- emmeans(model.cWM, pairwise ~ Time, adjust = "bonferroni")
#use to calculate effect sizes of contrasts
library(effectsize)
contrast_effect_sizes <- eff_size(posthoc.cWM$emmeans, sigma = sigma(model.cWM), edf = df.residual(model.cWM))


cWM.vol <-summarySE(data=data_crew,measurevar="CerebralWhiteMatterVol",groupvars=c("Time"),na.rm=T) 


#csf lmer
source("~/Dropbox/R/summary_se_function.R")
model.csf <- lmer(CSF ~ Time + (1 | subject), data = data_crew)
anova(model.csf)
posthoc.csf<- emmeans(model.csf, pairwise ~ Time, adjust = "bonferroni")
#use to calculate effect sizes of contrasts
library(effectsize)
contrast_effect_sizes <- eff_size(posthoc.csf$emmeans, sigma = sigma(model.csf), edf = df.residual(model.csf))


CSF.vol <-summarySE(data=data_crew,measurevar="CSF",groupvars=c("Time"),na.rm=T) 


#lateral ventricle lmer
model.ventricle <- lmer(Lateral_Ventricle ~ Time + (1 | subject), data = data_crew)
anova(model.ventricle)
posthoc.ventricle<- emmeans(model.ventricle, pairwise ~ Time, adjust = "bonferroni")
#use to calculate effect sizes of contrasts
library(effectsize)
contrast_effect_sizes_vent <- eff_size(posthoc.ventricle$emmeans, sigma = sigma(model.csf), edf = df.residual(model.csf))


LV.vol <-summarySE(data=data_crew,measurevar="Lateral_Ventricle",groupvars=c("Time"),na.rm=T) 



#3rd ventricle lmer
model.3rdventricle <- lmer(X3rd_Ventricle ~ Time + (1 | subject), data = data_crew)
anova(model.3rdventricle)
posthoc.3rdventricle<- emmeans(model.3rdventricle, pairwise ~ Time, adjust = "bonferroni")
contrast_effect_sizes_3rdvent <- eff_size(posthoc.3rdventricle$emmeans, sigma = sigma(model.csf), edf = df.residual(model.csf))


Third_Ventrilce.vol <-summarySE(data=data_crew,measurevar="X3rd_Ventricle",groupvars=c("Time"),na.rm=T) 


#4th ventricle lmer
model.4thventricle <- lmer(X4th_Ventricle ~ Time + (1 | subject), data = data_crew)
anova(model.4thventricle)
posthoc.4thventricle<- emmeans(model.4thventricle, pairwise ~ Time, adjust = "bonferroni")
contrast_effect_sizes_4thvent <- eff_size(posthoc.4thventricle$emmeans, sigma = sigma(model.csf), edf = df.residual(model.csf))


Fourth_Ventrilce.vol <-summarySE(data=data_crew,measurevar="X4th_Ventricle",groupvars=c("Time"),na.rm=T) 

#p.adjust for CSF and ventricles for crew
p.orig<-c(2.21E-01, 3.78E-06, 5.07E-07, 1.78E-04)
adjust_pvals_csf_vent<-p.adjust(p.orig, method = "bonferroni")
adjust_pvals_csf_vent

#5th ventricle lmer (spinal column) not used
model.5thventricle <- lmer(X5th_Ventricle ~ Time + (1 | subject), data = data_crew)
anova(model.5thventricle)
posthoc.5thventricle<- emmeans(model.5thventricle, pairwise ~ Time, adjust = "bonferroni")
Fifth_Ventrilce.vol <-summarySE(data=data_crew,measurevar="X5th_Ventricle",groupvars=c("Time"),na.rm=T) 



#plot spaghetti
#accumbnens
plot.spaghetti.AA <- ggplot(data_crew, aes(x = Time, y = Accumbens_area, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Accumbens",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")

plot.spaghetti.AA

#mean line only

meanline.acc<-ggplot(data_crew, aes(x = Time, y = Accumbens_area, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Accumbens",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()

#amygdala

plot.spaghetti.AY <- ggplot(data_crew, aes(x = Time, y = Amygdala, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Amygdala",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")


plot.spaghetti.AY

#mean line only
meanline.amy<-ggplot(data_crew, aes(x = Time, y = Amygdala, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Amygdala",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()

#caudate
plot.spaghetti.C <- ggplot(data_crew, aes(x = Time, y = Caudate, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Caudate",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")

plot.spaghetti.C

#mean line only
meanline.cau<-ggplot(data_crew, aes(x = Time, y = Caudate, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Caudate",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()

#hippocampus
plot.spaghetti.H <- ggplot(data_crew, aes(x = Time, y = Hippocampus, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Hippocampus",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")

plot.spaghetti.H

#mean line only
meanline.hipp<-ggplot(data_crew, aes(x = Time, y = Hippocampus, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Hippocampus",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()

#pallidum
plot.spaghetti.PL <- ggplot(data_crew, aes(x = Time, y = Pallidum, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Pallidum",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")


plot.spaghetti.PL

#mean line only
meanline.pall<-ggplot(data_crew, aes(x = Time, y = Pallidum, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Pallidum",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()


#putamen
plot.spaghetti.PM <- ggplot(data_crew, aes(x = Time, y = Putamen, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Putamen",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")

plot.spaghetti.PM

#mean line only
meanline.put<-ggplot(data_crew, aes(x = Time, y = Putamen, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Putamen",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()

#thalamus
plot.spaghetti.T <- ggplot(data_crew, aes(x = Time, y = Thalamus, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Thalamus",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")


#mean line only

meanline.thal<-ggplot(data_crew, aes(x = Time, y = Thalamus, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Thalamus",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()


plot.spaghetti.T

#parietal lobe
plot.spaghetti.PV <- ggplot(data_crew, aes(x = Time, y = parietal_vol, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Parietal",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")

plot.spaghetti.PV

#mean line only
meanline.pv<-ggplot(data_crew, aes(x = Time, y = parietal_vol, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Parietal Lobe",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()

#frontal
plot.spaghetti.FV <- ggplot(data_crew, aes(x = Time, y = frontal_vol, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Frontal",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")


plot.spaghetti.FV

#mean line only
meanline.fv<-ggplot(data_crew, aes(x = Time, y = frontal_vol, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Frontal Lobe",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()


#occipital
plot.spaghetti.OV <- ggplot(data_crew, aes(x = Time, y = occipital_vol, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Occipital",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")


plot.spaghetti.OV

#mean line only
meanline.ov<-ggplot(data_crew, aes(x = Time, y = occipital_vol, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Occipital Lobe",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()

#temporal
plot.spaghetti.TV <- ggplot(data_crew, aes(x = Time, y = temporal_vol, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Temporal",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")

plot.spaghetti.TV

#mean line only
meanline.tv<-ggplot(data_crew, aes(x = Time, y = temporal_vol, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Temporal Lobe",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()


#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalROIs_byTime-CREW.pdf")
pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalROIs_byTime-CONTROLS-ALL_012425.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalROIs_byTime-CONTROLS-completers_012425.pdf")

ggarrange(plot.spaghetti.AA, plot.spaghetti.AY, plot.spaghetti.C, plot.spaghetti.H, plot.spaghetti.PL,plot.spaghetti.PM, plot.spaghetti.T,
          labels = c("A", "B", "C", "D", "E", "F", "G"),
          ncol = 3, nrow = 3)
dev.off()

#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/MEAN-SubCorticalROIs_byTime-CREW.pdf")
pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/MEAN-SubCorticalROIs_byTime-CONTROLS-ALL_012425.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/MEAN-SubCorticalROIs_byTime-CONTROLS-completers_012425.pdf")

ggarrange(meanline.acc, meanline.amy, meanline.cau, meanline.hipp, meanline.pall,meanline.put, meanline.thal,
          #labels = c("E", "F", "G", "H", "I", "J", "K"),
          ncol = 3, nrow = 3)
dev.off()

#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalROIs_byTimeONLY_CREW.pdf")
pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalROIs_byTimeONLY_CONTROLS-ALL_012425.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalROIs_byTimeONLY_CONTROLS-completers_012425.pdf")

ggarrange(plot.spaghetti.FV, plot.spaghetti.TV, plot.spaghetti.PV, plot.spaghetti.OV,
          #labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()


#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/MEAN-CorticalROIs_byTimeONLY_CREW.pdf")
pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/MEAN-CorticalROIs_byTimeONLY_CONTROLS-ALL_012425.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/MEAN-CorticalROIs_byTimeONLY_CONTROLS-completers_012425.pdf")


ggarrange(meanline.fv, meanline.tv, meanline.pv, meanline.ov,
          #labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

dev.off()

#recomment line 332
#CSF

plot.spaghetti.CSF <- ggplot(data_crew, aes(x = Time, y = CSF, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "CSF",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")
#dev.off()


#mean line only
meanline.csf<-ggplot(data_crew, aes(x = Time, y = CSF, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "CSF",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()


#raw ICV
icv<-read.csv("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/Data/concordia_crew_control_demo_freesurfer7.4_longitudinal.csv")
icv<-icv[icv$group=="Crew",]
icv$Time<-as.factor(icv$Time)

model.icv <- lm(EstimatedTotalIntraCranialVol ~ scanner , data = icv)
anova(model.icv)
posthoc.icv<- emmeans(model.csf, pairwise ~ Time, adjust = "bonferroni")



#lateral ventricle

#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/Ventricle_byTimeONLY_CREW.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CSF_byTimeONLY_CONTROLS.pdf")

plot.spaghetti.LV <- ggplot(data_crew, aes(x = Time, y = Lateral_Ventricle, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Lateral Ventricle",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")
#dev.off()

#mean line only
meanline.LV<-ggplot(data_crew, aes(x = Time, y = Lateral_Ventricle, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Lateral Ventricle",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()


#3rd Ventricle
plot.spaghetti.3V <- ggplot(data_crew, aes(x = Time, y = X3rd_Ventricle, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "3rd Ventricle",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")

#mean line only
meanline.3rdV<-ggplot(data_crew, aes(x = Time, y = X3rd_Ventricle, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "3rd Ventricle",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()

#4th Ventricle

plot.spaghetti.4V <- ggplot(data_crew, aes(x = Time, y = X4th_Ventricle, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "4th Ventricle",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")

#mean line only
meanline.4thV<-ggplot(data_crew, aes(x = Time, y = X4th_Ventricle, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "4th Ventricle",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()

###added 04/01/2025
#cereberal white matter

plot.spaghetti.cWM <- ggplot(data_crew, aes(x = Time, y = CerebralWhiteMatterVol, group = subject, color = subject)) +
  geom_line() +
  geom_point() +
  ylab("X") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, aes(group = 1), color = "black") + # Confidence intervals
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, aes(group = 1), color = "black") + # Mean squares
  stat_summary(fun = mean, geom = "line", lwd = 1, aes(group = 1), color = "black") + # Mean line
  labs(
    title = "Cerebral White Matter",
    x = "Time",
    y = "Normalized Volume"
  ) +
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2")) +
  scale_color_grey(start = 0.8, end = 0.2) +  # Grayscale color mappingtheme_minimal() for controls
  theme_minimal() +
  theme(legend.position = "none")
#dev.off()


#mean line only
meanline.cWM<-ggplot(data_crew, aes(x = Time, y = CerebralWhiteMatterVol, group = 1)) + # Ensure the line connects the points
  stat_summary(fun = mean, geom = "line", color = "black", lwd = 1) + # Mean line
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) + # Mean points
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") + # Confidence intervals
  scale_x_discrete(labels = c("t0" = "Pre", "t12" = "Post1", "t18" = "Post2"))+
  labs(title = "Cerebral White Matter",
       x = "Time",
       y = "Normalized Volume") +
  theme_minimal()




pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CSF_and_LV_3V_4V_CREW.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CSF_and_LV_3V_4V_CREW_CONTROLS-ALL_012425.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CSF_and_LV_3V_4V_CREW_CONTROLS-completers_040125.pdf")

ggarrange(plot.spaghetti.CSF, plot.spaghetti.LV, plot.spaghetti.3V, plot.spaghetti.4V,plot.spaghetti.cWM,
          labels = c("L", "M", "N", "O", "P"),
          ncol = 2, nrow = 3)
dev.off()

pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CSF_and_LV_3V_4V_CREW_Mean.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CSF_and_LV_3V_4V_CREW_Mean_CONTROLS-ALL_040125.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CSF_and_LV_3V_4V_CREW_Mean_CONTROLS-completers_040125.pdf")

ggarrange(meanline.csf, meanline.LV, meanline.3rdV, meanline.4thV, meanline.cWM,
          #labels = c("L", "M", "N", "O", "P"),
          ncol = 2, nrow = 3)
dev.off()





#prepare for other sensitivity analyses

data_crew$TOD<-as.factor(data_crew$TOD)
data_crew$TOD<-as.factor(data_crew$TOD)
data_crew$Time<-as.factor(data_crew$Time)
data_crew$scanner<-as.factor(data_crew$scanner)
data_crew$PatientSex<-as.factor(data_crew$PatientSex)
data_crew$PatientSex<-ordered(data_crew$PatientSex, levels = c("M", "F"))

####time of day#####
tod_lm_temporal<-lm(temporal_vol~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_temporal)
tod_lm_frontal<-lm(frontal_vol~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_frontal)
tod_lm_occipital<-lm(occipital_vol~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_occipital)
tod_lm_parietal<-lm(parietal_vol~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_parietal)



# Load necessary library
library(ggplot2)
library(reshape2)
library(dplyr)


# Melt the data to make it long format for plotting
melted_data <- melt(data_crew, id.vars = "TOD", 
                    measure.vars = c("temporal_vol", "frontal_vol", "parietal_vol", "occipital_vol"),
                    variable.name = "ROI", value.name = "Volume")

melted_data <- melted_data[!(is.na(melted_data$TOD) | melted_data$TOD == ""), ]

# Relabel Volume_Type
melted_data <- melted_data %>%
  mutate(ROI = case_when(
    ROI == "temporal_vol" ~ "Temporal",
    ROI == "frontal_vol" ~ "Frontal",
    ROI == "parietal_vol" ~ "Parietal",
    ROI == "occipital_vol" ~ "Occipital",
    TRUE ~ ROI
  ))

summary_stats <- melted_data %>%
  group_by(TOD, ROI) %>%
  summarise(
    Mean = mean(Volume, na.rm = TRUE),
    SD = sd(Volume, na.rm = TRUE)
  )

# Define dodge position
dodge <- position_dodge(width = 0.5)


b<-ggplot(summary_stats, aes(x = TOD, y = Mean, color = ROI, group = ROI)) +
  geom_point(size = 3, position = dodge) +
  geom_line(size = 1, position = dodge) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2,  position = dodge) +
  labs(title = "Cortical Volumes by TOD", 
       x = "TOD", 
       y = "Mean Volume (±SD)", 
       color = "Volume Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalTOD.pdf")
b
dev.off()

tod_lm_thal<-lm(Thalamus~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_thal)
tod_lm_caud<-lm(Caudate~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_caud)
tod_lm_pallidum<-lm(Pallidum~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_pallidum)
tod_lm_putamen<-lm(Putamen~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_putamen)
tod_lm_hippocampus<-lm(Hippocampus~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_hippocampus)
tod_lm_amygdala<-lm(Amygdala~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_amygdala)
tod_lm_accumbens<-lm(Accumbens_area~TOD, data=data_crew, na.action=na.omit)
summary(tod_lm_accumbens)




# Melt the data to make it long format for plotting
melted_data.sub <- melt(data_crew, id.vars = "TOD", 
                    measure.vars = c("Accumbens_area", "Amygdala", "Caudate", "Hippocampus", "Pallidum", "Putamen", "Thalamus" ),
                    variable.name = "ROI", value.name = "Volume")

melted_data.sub <- melted_data.sub[!(is.na(melted_data.sub$TOD) | melted_data.sub$TOD == ""), ]

# Relabel Volume_Type
melted_data.sub <- melted_data.sub %>%
  mutate(ROI = case_when(
    ROI == "Accumbens_area" ~ "Accumbens",
    TRUE ~ ROI
  ))

summary_stats.sub <- melted_data.sub %>%
  group_by(TOD, ROI) %>%
  summarise(
    Mean = mean(Volume, na.rm = TRUE),
    SD = sd(Volume, na.rm = TRUE)
  )

# Define dodge position
dodge <- position_dodge(width = 0.5)


c<-ggplot(summary_stats.sub, aes(x = TOD, y = Mean, color = ROI, group = ROI)) +
  geom_point(size = 3, position = dodge) +
  geom_line(size = 1, position = dodge) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2,  position = dodge) +
  labs(title = "Subcortical Volumes by TOD", 
       x = "TOD", 
       y = "Mean Volume (±SD)", 
       color = "Volume Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalTOD.pdf")
c
dev.off()


#age
# Create a scatter plot with correlation lines
#frontal Accumbens_area
age.fv<-ggplot(data_crew, aes(x = PatientAgeYears, y = frontal_vol, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Frontal Lobe",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())

#Temproal
age.tv<-ggplot(data_crew, aes(x = PatientAgeYears, y = temporal_vol, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Temporal Lobe",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())

#Parietal
age.pv<-ggplot(data_crew, aes(x = PatientAgeYears, y = parietal_vol, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Parietal Lobe",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())

#Occipital
age.ov<-ggplot(data_crew, aes(x = PatientAgeYears, y = occipital_vol, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Occipital Lobe",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())


pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalROIs_byAge-CREW.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalROIs_byTimeONLY_CONTROLS.pdf")

ggarrange(age.fv, age.tv, age.pv, age.ov,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()


#accumbens Accumbens_area
age.acc<-ggplot(data_crew, aes(x = PatientAgeYears, y = Accumbens_area, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Accumbens",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())
#amygdala
age.amy<-ggplot(data_crew, aes(x = PatientAgeYears, y = Amygdala, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Amygdala",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())

#caudate
age.cau<-ggplot(data_crew, aes(x = PatientAgeYears, y = Caudate, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Caudate",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())

#hippocampus
age.hipp<-ggplot(data_crew, aes(x = PatientAgeYears, y = Hippocampus, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Hippocampus",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())

#pallidum
age.pal<-ggplot(data_crew, aes(x = PatientAgeYears, y = Pallidum, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Pallidum",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())

#putamen
age.put<-ggplot(data_crew, aes(x = PatientAgeYears, y = Putamen, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Putamen",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())

#thalamus
age.thal<-ggplot(data_crew, aes(x = PatientAgeYears, y = Thalamus, color = scanner)) +
  geom_point() +        # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add correlation (linear regression) lines
  labs(title = "Thalamus",
       x = "Participant Age (Years)",
       y = "Normalized Volume") +
  theme_minimal() +
  theme(legend.title = element_blank())


pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalROIs_byAge-CREW.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalROIs_byTime-CONTROLS.pdf")

ggarrange(age.acc, age.amy, age.cau, age.hipp, age.pal,age.put, age.thal,
          labels = c("E", "F", "G", "H", "I", "J", "K"),
          ncol = 3, nrow = 3)
dev.off()


#density plots
# Create density plots

#frontal
sex.density.fv<-ggplot(data_crew, aes(x = frontal_vol, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Frontal Lobe",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#temporal
sex.density.tv<-ggplot(data_crew, aes(x = temporal_vol, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Temporal Lobe",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#parietal
sex.density.pv<-ggplot(data_crew, aes(x = parietal_vol, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Parietal Lobe",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#occipital
sex.density.ov<-ggplot(data_crew, aes(x = occipital_vol, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Occipital Lobe",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())

pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalROIs_DensitybySex-CREW.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalROIs_byTimeONLY_CONTROLS.pdf")

ggarrange(sex.density.fv, sex.density.tv, sex.density.pv, sex.density.ov,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()

#accumbens
sex.density.acc<-ggplot(data_crew, aes(x = Accumbens_area, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Accumbens",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#amygdala
sex.density.amy<-ggplot(data_crew, aes(x = Amygdala, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Amygdala",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#caudate
sex.density.cau<-ggplot(data_crew, aes(x = Caudate, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Caudate",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#hippocampus
sex.density.hipp<-ggplot(data_crew, aes(x = Hippocampus, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Hippocampus",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#pallidum
sex.density.pall<-ggplot(data_crew, aes(x = Pallidum, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Pallidum",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#putamen
sex.density.put<-ggplot(data_crew, aes(x = Putamen, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Putamen",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#thalamus
sex.density.thal<-ggplot(data_crew, aes(x = Thalamus, fill = PatientSex)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Thalamus",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalROIs_DensitybySex-CREW.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalROIs_byTime-CONTROLS.pdf")

ggarrange(sex.density.acc, sex.density.amy, sex.density.cau, sex.density.hipp, sex.density.pall,sex.density.put, sex.density.thal,
          labels = c("E", "F", "G", "H", "I", "J", "K"),
          ncol = 3, nrow = 3)
dev.off()



#Density by Site
#frontal
site.density.fv<-ggplot(data_crew, aes(x = frontal_vol, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Frontal Lobe",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#temporal
site.density.tv<-ggplot(data_crew, aes(x = temporal_vol, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Temporal Lobe",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#parietal
site.density.pv<-ggplot(data_crew, aes(x = parietal_vol, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Parietal Lobe",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#occipital
site.density.ov<-ggplot(data_crew, aes(x = occipital_vol, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Occipital Lobe",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())

pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalROIs_DensitybySite-CREW.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/CorticalROIs_byTimeONLY_CONTROLS.pdf")

ggarrange(site.density.fv, site.density.tv, site.density.pv, site.density.ov,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()


#accumbens
site.density.acc<-ggplot(data_crew, aes(x = Accumbens_area, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Accumbens",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#amygdala
site.density.amy<-ggplot(data_crew, aes(x = Amygdala, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Amygdala",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#caudate
site.density.cau<-ggplot(data_crew, aes(x = Caudate, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Caudate",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#hippocampus
site.density.hipp<-ggplot(data_crew, aes(x = Hippocampus, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Hippocampus",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#pallidum
site.density.pall<-ggplot(data_crew, aes(x = Pallidum, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Pallidum",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#putamen
site.density.put<-ggplot(data_crew, aes(x = Putamen, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Putamen",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
#thalamus
site.density.thal<-ggplot(data_crew, aes(x = Thalamus, fill = scanner)) +
  geom_density(alpha = 0.6) +  # Add density plots with transparency
  labs(title = "Thalamus",
       x = "Normalized Volume",
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())

pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalROIs_DensitybySite-CREW.pdf")
#pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/SubCorticalROIs_byTime-CONTROLS.pdf")

ggarrange(site.density.acc, site.density.amy, site.density.cau, site.density.hipp, site.density.pall,site.density.put, site.density.thal,
          labels = c("E", "F", "G", "H", "I", "J", "K"),
          ncol = 3, nrow = 3)
dev.off()
#Euler plots
# Recode and order the Time variable
data <- data %>%
  mutate(Time = factor(recode(Time, 
                              "t0" = "Pre", 
                              "t12" = "Post1", 
                              "t18" = "Post2"), 
                       levels = c("Pre", "Post1", "Post2")))

# Create the box plot
pdf("~/Library/CloudStorage/Box-Box/Concordia_Manuscript/WorkingDraft/npjMicrogravity/Revisions/2ndRevision/EulerBoxPlots.pdf")
ggplot(data, aes(x = interaction(scanner, Time, sep = " - "), y = euler_total, fill = Time)) +
  geom_boxplot(outlier.shape =NA) +  # Add boxplot with styled outliers
  labs(title = "Box Plot of Data Quality by Site and Time",
       x = "Scanner and Time",
       y = "Euler Total") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity
dev.off()



################

#below includes nested full model. 
#convert wide to long
data_crew_wide<-data_crew[,c("subject","Time","PatientSex","PatientAgeYears","frontal_vol", "parietal_vol", "temporal_vol", "occipital_vol","Accumbens_area", "Amygdala", "Caudate", "Hippocampus",
  "Pallidum", "Putamen", "Thalamus")]
data_crew_long<-melt(data_crew_wide, id.vars=c("subject", "Time","PatientSex","PatientAgeYears"))
names(data_crew_long)[5:6]<-c("roi","vol")
data_crew_long$roi<-factor(data_crew_long$roi)
data_crew_long$Time<-factor(data_crew_long$Time)

#not significant
lmer_longcortical1<- lmer(formula(paste("vol~ Time*roi +(1|subject) ")), data=data_crew_long)
emmeans(lmer_longcortical1,pirwise ~ Time|roi)

#when specifying nested, it is significant, however the q is whether roi can be treated as hierachcial or nested
lmer_longcortical2<- lmer(formula(paste("vol~ Time*roi +(1|subject/roi) ")), data=data_crew_long)

anova(lmer_longcortical1,lmer_longcortical2)

emmeans(lmer_longcortical2,pairwise ~ Time|roi)


#adding Patient Age as covariate- results stay the same as expected, sex is same over time but age does change

# Time and PatientAge- are highly correlated but variability amongst timepoints exists in the crew

out.lmer2 <- data.frame(Region=rois, F=rep(NA, length(rois)),P=rep(NA, length(rois)))
out.lmer2.posthoc<-as.list(1:length(rois))
names(out.lmer2.posthoc) <- rois

for (roi in rois) {
  f <- formula(paste(roi,"~Time + PatientAgeYears+PatientSex+(1|subject)"))
  lmer.m2<- lmer(f, data=data_crew)
    out.lmer2[out.lmer2$Region == roi, "F"] <- anova(lmer.m2)[[5]][1]
    out.lmer2[out.lmer2$Region == roi, "P"] <- anova(lmer.m2)[[6]][1]
    posthoc.lmer2<- emmeans(lmer.m2, pairwise ~ Time, adjust = "bonferroni")
    posthoc.lmer2.contrasts<-as.data.frame(posthoc.lmer2$contrasts)
    posthoc.lmer2.contrasts$roi<-roi
    out.lmer2.posthoc[[roi]]<-posthoc.lmer2.contrasts
}

out.lmer2$P_FDR <- round(p.adjust(out.lmer2$P, method='fdr'),4)
out.lmer2$P_Bonferroni <- round(p.adjust(out.lmer2$P, method='bonferroni'),4)

