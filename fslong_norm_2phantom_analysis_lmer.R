#lmer
library(ggeffects)
library(lme4)
library(emmeans)
library(pbkrtest)
library(ggplot2)
library(visreg)
library(lmerTest)
#getting pvalues out of lme4
#car::Anova and lmerTest::anova provide wrappers for Kenward-Roger-corrected tests using pbkrtest: lmerTest::anova also provides t tests via the Satterthwaite approximation (P,*)
#read data
data<-read.csv("fslong_normed_2phantomsbyscanner_crew_control.csv")

#crew data
data_crew<-data[data$group=="Crew",] #67

#regions
rois=c("Accumbens_area", "Amygdala", "Caudate", "Hippocampus",
  "Pallidum", "Putamen", "Thalamus", "frontal_vol", "parietal_vol", "temporal_vol", "occipital_vol","CSF")

data_crew$Time<-as.factor(data_crew$Time)

#listc <- vector(mode = "list")

#another univariate analyses needs to run with Time and PatientAge- are highly correlated but variability amongst timepoints exists in the crew
pdf("lmer_output_crew_fslong_norm2phantom.pdf")

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
dev.off()

out.lmer1$P_FDR <- round(p.adjust(out.lmer1$P, method='fdr'),4)
out.lmer1$P_Bonferroni <- round(p.adjust(out.lmer1$P, method='bonferroni'),4)

out.lmer1.posthoc
#plot spaghetti
plot.spaghetti.H<-ggplot(data_crew, aes (x = Time, y = Hippocampus, group = subject, color = subject))
plot.spaghetti.H <-plot.spaghetti.H + geom_line() + geom_point() + ylab("X") + stat_summary(fun.y=mean,geom="line",lwd=2,aes(group=1))+
  theme_minimal() +theme(legend.position = "none")+  labs(title = "ROI Volume by Time",x = "Time",y = "Volume")

plot.spaghetti.H

plot.spaghetti.P<-ggplot(data_crew, aes (x = Time, y = Pallidum, group = subject, color = subject))
plot.spaghetti.P <-plot.spaghetti.P + geom_line() + geom_point() + ylab("X") + stat_summary(fun.y=mean,geom="line",lwd=2,aes(group=1))+
  theme_minimal() +theme(legend.position = "none")+  labs(title = "ROI Volume by Time",x = "Time",y = "Volume")

plot.spaghetti.P

plot.spaghetti.T<-ggplot(data_crew, aes (x = Time, y = Thalamus, group = subject, color = subject))
plot.spaghetti.T <-plot.spaghetti.T + geom_line() + geom_point() + ylab("X") + stat_summary(fun.y=mean,geom="line",lwd=2,aes(group=1))+
  theme_minimal() +theme(legend.position = "none")+  labs(title = "ROI Volume by Time",x = "Time",y = "Volume")

plot.spaghetti.T


plot.spaghetti.PV<-ggplot(data_crew, aes (x = Time, y = parietal_vol, group = subject, color = subject))
plot.spaghetti.PV<-plot.spaghetti.PV + geom_line() + geom_point() + ylab("X") + stat_summary(fun.y=mean,geom="line",lwd=2,aes(group=1))+
  theme_minimal() +theme(legend.position = "none")+  labs(title = "ROI Volume by Time",x = "Time",y = "Volume")

plot.spaghetti.PV

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


#bar plots

