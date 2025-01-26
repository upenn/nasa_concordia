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
  "Pallidum", "Putamen", "Thalamus", "frontal_vol", "parietal_vol", "temporal_vol", "occipital_vol","CSF","Lateral_Ventricle","X3rd_Ventricle","X4th_Ventricle","X5th_Ventricle")

data_crew$Time<-as.factor(data_crew$Time)


out.lmer1 <- data.frame(Region=rois, F=rep(NA, length(rois)),P=rep(NA, length(rois)))
out.lmer1.posthoc<-as.list(1:length(rois))
names(out.lmer1.posthoc) <- rois
#Effect sizes for metric data can be calculated with r = √(t²/(t^2+df)) (Rosenthal, 1991, p. 19)
#lower 0.1: no effect
#0.1-0.29: small effect
##0.3-0.49: medium effect
#0.5-1: large effect
#Rosenthal, R. 1991. Meta-analytic Procedures for Social Research. 2nd ed., Sage Publications, Newbury Park. Cohen J 1992 A power primer. Psychological Bulletin 112: 155-159.

library(effectsize) #https://cran.r-project.org/web/packages/effectsize/vignettes/from_test_statistics.html
# issue raised by author of emmeans - calculating effect size
options(es.use_symbols = TRUE)
for (roi in rois) {
  f <- formula(paste(roi,"~Time + (1|subject)"))
  lmer.m1<- lmer(f, data=data_crew)
    out.lmer1[out.lmer1$Region == roi, "F"] <- anova(lmer.m1)[[5]]
    out.lmer1[out.lmer1$Region == roi, "P"] <- anova(lmer.m1)[[6]]
    posthoc.lmer1<- emmeans(lmer.m1, pairwise ~ Time, adjust = "bonferroni")
    posthoc.lmer1.contrasts<-as.data.frame(posthoc.lmer1$contrasts)
    posthoc.lmer1.contrasts$roi<-roi
    out.lmer1.posthoc[[roi]]<-posthoc.lmer1.contrasts
    #effect size table for each pairwise contrast
     posthoc.t.effectsizeest<-as.data.frame(t_to_eta2(as.data.frame(posthoc.lmer1$contrasts)$t.ratio,as.data.frame(posthoc.lmer1$contrasts[1:3])$df))
     posthoc.f.effectsizeest<-as.data.frame(F_to_eta2(anova(lmer.m1)[1,5],anova(lmer.m1)[1,3],anova(lmer.m1)[1,4]))
    print(roi)
    print("------------Tvalue Effectsize----------------------------")
    print(posthoc.t.effectsizeest)
    print("------------Fvalue Effectsize----------------------------")
    print(posthoc.f.effectsizeest)
    #emmeans to compute effect size
    #eff_size(posthoc.lmer1$emmeans, sigma=sigma(lmer.m1), edf = df.residual(lmer.m1)) #not sure which sigma - residual sigma is incorrect here
    #https://stats.stackexchange.com/questions/509283/what-df-do-you-use-for-effect-size-calculation-using-emmeans-with-a-crossed-mixe
  
    out.em1 <-emmeans(lmer.m1, "Time")

#visreg(lmer.m1,"Time")
}
dev.off()

