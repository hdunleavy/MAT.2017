setwd("~/Documents/Advisor.meetings/Summer2017Fieldwork/MAT")


##### Root tip enzymes by treatment #########
rt.enz <- read.csv('MAT.enzymes.csv')
library(lme4)
library(nlme)
library(MuMIn)
library(emmeans)
#  Root enzyme activity models 
#LAP, sig interaction
mixed <- lme(fixed = log(LAP.mm2.min+1)~Warmed*Fertilized,
    data = rt.enz,
    random = ~1|Plot/Core,
    na.action= na.omit)
mixed.no.int <- lme(fixed = log(LAP.mm2.min+1)~Warmed+Fertilized,
             data = rt.enz,
             random = ~1|Plot/Core,
             na.action= na.omit)
reduced <- lme(fixed = log(LAP.mm2.min+1)~1,
                    data = rt.enz,
                    random = ~1|Plot/Core,
                    na.action= na.omit)
summary(mixed)
plot(mixed)
r.squaredGLMM(mixed)
leastsquare.LAP.mm <- lsmeans(mixed, 
                              pairwise ~ Warmed*Fertilized, 
                              adjust = "tukey")
CLD.LAP.mm <- cld(leastsquare.LAP.mm$lsmeans,
             alpha=0.05,
             Letters=letters,
             adjust="tukey")

#NAG, sig interaction
mixed <- lme(fixed = log(NAG.mm2.min+1)~Warmed*Fertilized,
             data = rt.enz,
             random = ~1|Plot/Core,
             na.action= na.omit)
mixed.no.int <- lme(fixed = log(NAG.mm2.min+1)~Warmed+Fertilized,
                    data = rt.enz,
                    random = ~1|Plot/Core,
                    na.action= na.omit)
reduced <- lme(fixed = log(NAG.mm2.min+1)~1,
               data = rt.enz,
               random = ~1|Plot/Core,
               na.action= na.omit)
summary(mixed)
plot(mixed)
r.squaredGLMM(mixed)
leastsquare.NAG.mm <- lsmeans(mixed, 
                              pairwise ~ Warmed*Fertilized, 
                              adjust = "tukey")
CLD.NAG.mm <- cld(leastsquare.NAG.mm$lsmeans,
                  alpha=0.05,
                  Letters=letters,
                  adjust="tukey")

#PHOS, sig interaction
mixed <- lme(fixed = log(PHOS.mm2.min+1)~Warmed*Fertilized,
             data = rt.enz,
             random = ~1|Plot/Core,
             na.action= na.omit)
mixed.no.int <- lme(fixed = log(PHOS.mm2.min+1)~Warmed+Fertilized,
                    data = rt.enz,
                    random = ~1|Plot/Core,
                    na.action= na.omit)
reduced <- lme(fixed = log(PHOS.mm2.min+1)~1,
               data = rt.enz,
               random = ~1|Plot/Core,
               na.action= na.omit)
summary(mixed)
plot(mixed)
r.squaredGLMM(mixed)
leastsquare.PHOS.mm <- lsmeans(mixed, 
                              pairwise ~ Warmed*Fertilized, 
                              adjust = "tukey")
CLD.PHOS.mm <- cld(leastsquare.PHOS.mm$lsmeans,
                  alpha=0.05,
                  Letters=letters,
                  adjust="tukey")

#ABTS, sig interaction
mixed <- lme(fixed = log(ABTS.mm2.min+1)~Warmed*Fertilized,
             data = rt.enz,
             random = ~1|Plot/Core,
             na.action= na.omit)
mixed.no.int <- lme(fixed = log(ABTS.mm2.min+1)~Warmed+Fertilized,
                    data = rt.enz,
                    random = ~1|Plot/Core,
                    na.action= na.omit)
reduced <- lme(fixed = log(ABTS.mm2.min+1)~1,
               data = rt.enz,
               random = ~1|Plot/Core,
               na.action= na.omit)
summary(mixed)
plot(mixed)
r.squaredGLMM(mixed)
leastsquare.ABTS.mm <- lsmeans(mixed, 
                               pairwise ~ Warmed*Fertilized, 
                               adjust = "tukey")
CLD.ABTS.mm <- cld(leastsquare.ABTS.mm$lsmeans,
                   alpha=0.05,
                   Letters=letters,
                   adjust="tukey")

#TMB, sig 
mixed <- lme(fixed = log(TMB.mm2.min+1)~Warmed*Fertilized,
             data = rt.enz,
             random = ~1|Plot/Core,
             na.action= na.omit)
mixed.no.int <- lme(fixed = log(TMB.mm2.min+1)~Warmed+Fertilized,
                    data = rt.enz,
                    random = ~1|Plot/Core,
                    na.action= na.omit)
reduced <- lme(fixed = log(TMB.mm2.min+1)~1,
               data = rt.enz,
               random = ~1|Plot/Core,
               na.action= na.omit)
summary(mixed.no.int)
plot(mixed.no.int)
r.squaredGLMM(mixed.no.int)
leastsquare.TMB.mm <- lsmeans(mixed.no.int, 
                               pairwise ~ Warmed*Fertilized, 
                               adjust = "tukey")
CLD.TMB.mm <- cld(leastsquare.TMB.mm$lsmeans,
                   alpha=0.05,
                   Letters=letters,
                   adjust="tukey")



##### Community enzymes by treatment #########
mat.summary <- read.csv('MAT.summarytable.csv')
library(lme4)
library(nlme)
library(MuMIn)
library(lsmeans)
#  Community enzyme activity models 
#LAP, sig interaction
mixed <- lme(fixed = LAP.cm3.min~Warmed*Fertilized,
             data = mat.summary,
             random = ~1|Plot,
             na.action= na.omit)
mixed.no.int <- lme(fixed = LAP.cm3.min~Warmed+Fertilized,
                    data = mat.summary,
                    random = ~1|Plot,
                    na.action= na.omit)
reduced <- lme(fixed = LAP.cm3.min~1,
               data = mat.summary,
               random = ~1|Plot,
               na.action= na.omit)
summary(mixed)
plot(mixed)
r.squaredGLMM(mixed)
leastsquare.LAP.c <- lsmeans(mixed, 
                              pairwise ~ Warmed*Fertilized, 
                              adjust = "tukey")
CLD.LAP.c <- cld(leastsquare.LAP.c,
                  alpha=0.05,
                  Letters=letters,
                  adjust="tukey")

#NAG, sig
mixed <- lme(fixed = log(NAG.cm3.min+1)~Warmed*Fertilized,
             data = mat.summary,
             random = ~1|Plot,
             na.action= na.omit)
mixed.no.int <- lme(fixed = log(NAG.cm3.min+1)~Warmed+Fertilized,
                    data = mat.summary,
                    random = ~1|Plot,
                    na.action= na.omit)
reduced <- lme(fixed = log(NAG.cm3.min+1)~1,
               data = mat.summary,
               random = ~1|Plot,
               na.action= na.omit)
summary(mixed.no.int)
plot(mixed.no.int)
r.squaredGLMM(mixed.no.int)
leastsquare.NAG.c <- lsmeans(mixed.no.int, 
                             pairwise ~ Warmed*Fertilized, 
                             adjust = "tukey")
CLD.NAG.c <- cld(leastsquare.NAG.c,
                 alpha=0.05,
                 Letters=letters,
                 adjust="tukey")

#PHOS, sig
mixed <- lme(fixed = log(PHOS.cm3.min+1)~Warmed*Fertilized,
             data = mat.summary,
             random = ~1|Plot,
             na.action= na.omit)
mixed.no.int <- lme(fixed = log(PHOS.cm3.min+1)~Warmed+Fertilized,
                    data = mat.summary,
                    random = ~1|Plot,
                    na.action= na.omit)
reduced <- lme(fixed = log(PHOS.cm3.min+1)~1,
               data = mat.summary,
               random = ~1|Plot,
               na.action= na.omit)
summary(mixed.no.int)
plot(mixed.no.int)
r.squaredGLMM(mixed.no.int)
leastsquare.PHOS.c <- lsmeans(mixed.no.int, 
                             pairwise ~ Warmed*Fertilized, 
                             adjust = "tukey")
CLD.PHOS.c <- cld(leastsquare.PHOS.c,
                 alpha=0.05,
                 Letters=letters,
                 adjust="tukey")

#ABTS, sig interaction
mixed <- lme(fixed = log(ABTS.cm3.min+1)~Warmed*Fertilized,
             data = mat.summary,
             random = ~1|Plot,
             na.action= na.omit)
mixed.no.int <- lme(fixed = log(ABTS.cm3.min+1)~Warmed+Fertilized,
                    data = mat.summary,
                    random = ~1|Plot,
                    na.action= na.omit)
reduced <- lme(fixed = log(ABTS.cm3.min+1)~1,
               data = mat.summary,
               random = ~1|Plot,
               na.action= na.omit)
summary(mixed)
plot(mixed)
r.squaredGLMM(mixed)
leastsquare.ABTS.c <- lsmeans(mixed, 
                              pairwise ~ Warmed*Fertilized, 
                              adjust = "tukey")
CLD.ABTS.c <- cld(leastsquare.ABTS.c,
                  alpha=0.05,
                  Letters=letters,
                  adjust="tukey")

#TMB, sig 
mixed <- lme(fixed = log(TMB.cm3.min+1)~Warmed*Fertilized,
             data = mat.summary,
             random = ~1|Plot,
             na.action= na.omit)
mixed.no.int <- lme(fixed = log(TMB.cm3.min+1)~Warmed+Fertilized,
                    data = mat.summary,
                    random = ~1|Plot,
                    na.action= na.omit)
reduced <- lme(fixed = log(TMB.cm3.min+1)~1,
               data = mat.summary,
               random = ~1|Plot,
               na.action= na.omit)
summary(mixed.no.int)
plot(mixed.no.int)
r.squaredGLMM(mixed.no.int)
leastsquare.TMB.c <- lsmeans(mixed.no.int, 
                              pairwise ~ Warmed*Fertilized, 
                              adjust = "tukey")
CLD.TMB.c <- cld(leastsquare.TMB.c,
                  alpha=0.05,
                  Letters=letters,
                  adjust="tukey")



##### Root tip enzyme plots #####
library(ggplot2)
library(ggpubr)
library(dplyr)
   ###### DF making #####
LAP.mm.df <- data.frame(c(1:4))
LAP.mm.df$means <- CLD.LAP.mm$lsmean
LAP.mm.df$SE <- CLD.LAP.mm$SE
LAP.mm.df$warmed <- CLD.LAP.mm$Warmed
LAP.mm.df$fert <- CLD.LAP.mm$Fertilized
LAP.mm.df <- LAP.mm.df %>% arrange(warmed, fert)

NAG.mm.df <- data.frame(c(1:4))
NAG.mm.df$means <- CLD.NAG.mm$lsmean
NAG.mm.df$SE <- CLD.NAG.mm$SE
NAG.mm.df$warmed <- CLD.NAG.mm$Warmed
NAG.mm.df$fert <- CLD.NAG.mm$Fertilized
NAG.mm.df <- NAG.mm.df %>% arrange(warmed, fert)

PHOS.mm.df <- data.frame(c(1:4))
PHOS.mm.df$means <- CLD.PHOS.mm$lsmean
PHOS.mm.df$SE <- CLD.PHOS.mm$SE
PHOS.mm.df$warmed <- CLD.PHOS.mm$Warmed
PHOS.mm.df$fert <- CLD.PHOS.mm$Fertilized
PHOS.mm.df <- PHOS.mm.df %>% arrange(warmed, fert)

ABTS.mm.df <- data.frame(c(1:4))
ABTS.mm.df$means <- CLD.ABTS.mm$lsmean
ABTS.mm.df$SE <- CLD.ABTS.mm$SE
ABTS.mm.df$warmed <- CLD.ABTS.mm$Warmed
ABTS.mm.df$fert <- CLD.ABTS.mm$Fertilized
ABTS.mm.df <- ABTS.mm.df %>% arrange(warmed, fert)

TMB.mm.df <- data.frame(c(1:4))
TMB.mm.df$means <- CLD.TMB.mm$lsmean
TMB.mm.df$SE <- CLD.TMB.mm$SE
TMB.mm.df$warmed <- CLD.TMB.mm$Warmed
TMB.mm.df$fert <- CLD.TMB.mm$Fertilized
TMB.mm.df <- TMB.mm.df %>% arrange(warmed, fert)

MAT.hydro <- data.frame(c(1:12))
MAT.hydro$means[1:4] <- LAP.mm.df$means
MAT.hydro$means[5:8] <- NAG.mm.df$means
MAT.hydro$means[9:12] <- PHOS.mm.df$means
MAT.hydro$se[1:4] <- LAP.mm.df$SE
MAT.hydro$se[5:8] <- NAG.mm.df$SE
MAT.hydro$se[9:12] <- PHOS.mm.df$SE
MAT.hydro$Treatment <- rep(c("Control", "Fertilized", "Warmed", "Warmed & Fertilized"), 3)
MAT.hydro$Enzyme[1:4] <- rep("Leucine aminopeptidase", 4)
MAT.hydro$Enzyme[5:8] <- rep("Chitinase", 4)
MAT.hydro$Enzyme[9:12] <- rep("Phosphatase", 4)
MAT.hydro$Enzyme <- factor(MAT.hydro$Enzyme, levels = c("Leucine aminopeptidase", "Chitinase", "Phosphatase"))

MAT.oxid <- data.frame(c(1:12))
MAT.oxid$means[1:4] <- ABTS.mm.df$means
MAT.oxid$means[5:8] <- TMB.mm.df$means
MAT.oxid$means[9:12] <- rep(NA,4)
MAT.oxid$se[1:4] <- ABTS.mm.df$SE
MAT.oxid$se[5:8] <- TMB.mm.df$SE
MAT.oxid$se[9:12] <- rep(NA,4)
MAT.oxid$Treatment <- rep(c("Control", "Fertilized", "Warmed", "Warmed & Fertilized"), 3)
MAT.oxid$Enzyme[1:4] <- rep("Phenol oxidase", 4)
MAT.oxid$Enzyme[5:8] <- rep("Peroxidase", 4)
MAT.oxid$Enzyme[9:12] <- rep("x",4)
MAT.oxid$Enzyme <- factor(MAT.oxid$Enzyme, levels = c("Phenol oxidase", "Peroxidase", "x"))


   ###### plot making #####
p.hydro <- ggplot(MAT.hydro, aes(x=Enzyme, y=means, 
                             ymin=means-se, ymax=means+se, fill = Treatment))+
  geom_bar(stat="identity", position = position_dodge(.95), color = "black") + 
  geom_errorbar(position = position_dodge(0.95), width = 0.4, color="black")+
  labs(title = "Root tip enzymes", x = "Treatment", y = "pmol/min/mm2") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values =  c("Control" = "#2b83ba", "Fertilized" = "#abdda4", "Warmed" = "#d73027", "Warmed & Fertilized" = "#ffffbf"),
                    guide = F)+
  ylim(0,7)
p.oxid <- ggplot(MAT.oxid, aes(x=Enzyme, y=means, 
                                 ymin=(means-se), ymax=(means+se), fill = Treatment))+
  geom_bar(stat="identity", position = position_dodge(.95), color = "black") + 
  geom_errorbar( color="black",  position = position_dodge(.95), width = 0.4)+
  labs(title = "Root tip enzymes", x = "Treatment", y = "pmol/min/mm2") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values =  c("Control" = "#2b83ba", "Fertilized" = "#abdda4", "Warmed" = "#d73027", "Warmed & Fertilized" = "#ffffbf"), guide = F)+
  ylim(0,7)
ggarrange(p.hydro,p.oxid)


##### Community enzyme plots #####
  ##### DF making ######
LAP.c.df <- data.frame(c(1:4))
LAP.c.df$means <- CLD.LAP.c$lsmean
LAP.c.df$SE <- CLD.LAP.c$SE
LAP.c.df$warmed <- CLD.LAP.c$Warmed
LAP.c.df$fert <- CLD.LAP.c$Fertilized
LAP.c.df <- LAP.c.df %>% arrange(warmed, fert)

NAG.c.df <- data.frame(c(1:4))
NAG.c.df$means <- CLD.NAG.c$lsmean
NAG.c.df$SE <- CLD.NAG.c$SE
NAG.c.df$warmed <- CLD.NAG.c$Warmed
NAG.c.df$fert <- CLD.NAG.c$Fertilized
NAG.c.df <- NAG.c.df %>% arrange(warmed, fert)

PHOS.c.df <- data.frame(c(1:4))
PHOS.c.df$means <- CLD.PHOS.c$lsmean
PHOS.c.df$SE <- CLD.PHOS.c$SE
PHOS.c.df$warmed <- CLD.PHOS.c$Warmed
PHOS.c.df$fert <- CLD.PHOS.c$Fertilized
PHOS.c.df <- PHOS.c.df %>% arrange(warmed, fert)

ABTS.c.df <- data.frame(c(1:4))
ABTS.c.df$means <- CLD.ABTS.c$lsmean
ABTS.c.df$SE <- CLD.ABTS.c$SE
ABTS.c.df$warmed <- CLD.ABTS.c$Warmed
ABTS.c.df$fert <- CLD.ABTS.c$Fertilized
ABTS.c.df <- ABTS.c.df %>% arrange(warmed, fert)

TMB.c.df <- data.frame(c(1:4))
TMB.c.df$means <- CLD.TMB.c$lsmean
TMB.c.df$SE <- CLD.TMB.c$SE
TMB.c.df$warmed <- CLD.TMB.c$Warmed
TMB.c.df$fert <- CLD.TMB.c$Fertilized
TMB.c.df <- TMB.c.df %>% arrange(warmed, fert)

MAT.hydro.c <- data.frame(c(1:12))
MAT.hydro.c$means[1:4] <- LAP.c.df$means
MAT.hydro.c$means[5:8] <- NAG.c.df$means
MAT.hydro.c$means[9:12] <- PHOS.c.df$means
MAT.hydro.c$se[1:4] <- LAP.c.df$SE
MAT.hydro.c$se[5:8] <- NAG.c.df$SE
MAT.hydro.c$se[9:12] <- PHOS.c.df$SE
MAT.hydro.c$Treatment <- rep(c("Control", "Fertilized", "Warmed", "Warmed & Fertilized"), 3)
MAT.hydro.c$Enzyme[1:4] <- rep("Leucine aminopeptidase", 4)
MAT.hydro.c$Enzyme[5:8] <- rep("Chitinase", 4)
MAT.hydro.c$Enzyme[9:12] <- rep("Phosphatase", 4)
MAT.hydro.c$Enzyme <- factor(MAT.hydro.c$Enzyme, levels = c("Leucine aminopeptidase", "Chitinase", "Phosphatase"))

MAT.oxid.c <- data.frame(c(1:12))
MAT.oxid.c$means[1:4] <- ABTS.c.df$means
MAT.oxid.c$means[5:8] <- TMB.c.df$means
MAT.oxid.c$means[9:12] <- rep(NA,4)
MAT.oxid.c$se[1:4] <- ABTS.c.df$SE
MAT.oxid.c$se[5:8] <- TMB.c.df$SE
MAT.oxid.c$se[9:12] <- rep(NA,4)
MAT.oxid.c$Treatment <- rep(c("Control", "Fertilized", "Warmed", "Warmed & Fertilized"), 3)
MAT.oxid.c$Enzyme[1:4] <- rep("Phenol oxidase", 4)
MAT.oxid.c$Enzyme[5:8] <- rep("Peroxidase", 4)
MAT.oxid.c$Enzyme[9:12] <- rep("x",4)
MAT.oxid.c$Enzyme <- factor(MAT.oxid.c$Enzyme, levels = c("Phenol oxidase", "Peroxidase", "x"))
  ##### plot making ######

p.hydro.c <- ggplot(MAT.hydro.c, aes(x=Enzyme, y=means, 
                                 ymin=means-se, ymax=means+se, fill = Treatment))+
  geom_bar(stat="identity", position = position_dodge(.95), color = "black") + 
  geom_errorbar( color="black",  position = position_dodge(.95), width = 0.4)+
  labs(title = "Community enzymes", x = "Treatment", y = "pmol/min/cm3") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values =  c("Control" = "#2b83ba", "Fertilized" = "#abdda4", "Warmed" = "#d73027", "Warmed & Fertilized" = "#ffffbf"),
                    guide = F)+
  ylim(0,7.5)
p.oxid.c <- ggplot(MAT.oxid.c, aes(x=Enzyme, y=means, 
                               ymin=(means-se), ymax=(means+se), fill = Treatment))+
  geom_bar(stat="identity", position = position_dodge(.95), color = "black") + 
  geom_errorbar(color="black",  position = position_dodge(.95), width = 0.4)+
  labs(title = "Community enzymes", x = "Treatment", y = "pmol/min/cm3") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values =  c("Control" = "#2b83ba", "Fertilized" = "#abdda4", "Warmed" = "#d73027", "Warmed & Fertilized" = "#ffffbf"), guide = F)+
  ylim(0,7.5)

ggarrange(p.hydro.c,p.oxid.c)

###### one big plot #######

MAT.c <- data.frame(c(1:20))
MAT.c$means <- c(1:20)
MAT.c$se <- c(1:20)
MAT.c$Enzyme <- c(1:20)
MAT.c$means[1:12] <- MAT.hydro.c$means
MAT.c$means[13:20] <- MAT.oxid.c$means[1:8]
MAT.c$se[1:12] <- MAT.hydro.c$se
MAT.c$se[13:20] <- MAT.oxid.c$se[1:8]
MAT.c$Enzyme <- rep(c("Leucine aminopeptidase", "Chitinase", "Phosphatase","Phenol oxidase", "Peroxidase"), each = 4)
MAT.c$Enzyme <- factor(MAT.c$Enzyme, levels = c("Leucine aminopeptidase", "Chitinase", "Phosphatase","Phenol oxidase", "Peroxidase"))
MAT.c$Treatment <- rep(c("Control", "Fertilized", "Warmed", "Warmed & Fertilized"), 5)

p.c <- ggplot(MAT.c, aes(x=Enzyme, y=means, 
                                     ymin=means-se, ymax=means+se, fill = Treatment))+
  geom_bar(stat="identity", position = position_dodge(.9), color = "black") + 
  geom_errorbar( color="black",  position = position_dodge(.9), width = 0.4)+
  labs(title = "Community enzymes", x = "Enzyme", y = expression(paste("log(Enzyme activity [",pmol/min/cm^{3},"] )"))) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values =  c("Control" = "#2b83ba", "Fertilized" = "#abdda4", "Warmed" = "#d73027", "Warmed & Fertilized" = "#ffffbf"),
                    guide = F)+
  ylim(0,7.5)

p.c
MAT.mm <- data.frame(c(1:20))
MAT.mm$means <- c(1:20)
MAT.mm$se <- c(1:20)
MAT.mm$Enzyme <- c(1:20)
MAT.mm$means[1:12] <- MAT.hydro$means
MAT.mm$means[13:20] <- MAT.oxid$means[1:8]
MAT.mm$se[1:12] <- MAT.hydro$se
MAT.mm$se[13:20] <- MAT.oxid$se[1:8]
MAT.mm$Enzyme <- rep(c("Leucine aminopeptidase", "Chitinase", "Phosphatase","Phenol oxidase", "Peroxidase"), each = 4)
MAT.mm$Enzyme <- factor(MAT.mm$Enzyme, levels = c("Leucine aminopeptidase", "Chitinase", "Phosphatase","Phenol oxidase", "Peroxidase"))
MAT.mm$Treatment <- rep(c("Control", "Fertilized", "Warmed", "Warmed & Fertilized"), 5)

p.mm <- ggplot(MAT.mm, aes(x=Enzyme, y=means, 
                         ymin=means-se, ymax=means+se, fill = Treatment))+
  geom_bar(stat="identity", position = position_dodge(.9), color = "black") + 
  geom_errorbar( color="black",  position = position_dodge(.9), width = 0.4)+
  labs(title = "Root tip enzymes", x = "Enzyme", y = expression(paste("log(Enzyme activity [",pmol/min/mm^{2},"] )"))) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values =  c("Control" = "#2b83ba", "Fertilized" = "#abdda4", "Warmed" = "#d73027", "Warmed & Fertilized" = "#ffffbf"))+
  ylim(0,7.5)+
  theme(legend.position = c(0.15,0.9),
        legend.key.size = unit(15, "point"), 
        legend.title = element_text(NA))

p.mm

### just fertilized
MAT.mm.fert <- MAT.mm[MAT.mm$Treatment == "Control" |MAT.mm$Treatment == "Fertilized",]
MAT.mm.fert$letters <- c("a","b","a","b","a","b","a","b","a","b")
p.fert.mm <- ggplot(MAT.mm.fert, aes(x=Enzyme, y=means, 
                           ymin=means-se, ymax=means+se, fill = Treatment))+
  geom_bar(stat="identity", position = position_dodge(.9), color = "black") + 
  geom_errorbar( color="black",  position = position_dodge(.9), width = 0.4)+
  geom_text(position = position_dodge(0.9), 
            aes(x = Enzyme, y = means+se+0.2, label = letters))+
  labs(title = "28-year experiment", x = "Enzyme", y = expression(paste("log(Enzyme activity [",pmol/min/mm^{2},"] )"))) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values =  c("Control" = "#2b83ba", "Fertilized" = "#abdda4"))+
  ylim(0,7.5)+
  theme(legend.position = "bottom",
        legend.key.size = unit(15, "point"), 
        legend.title = element_text(NA))+
  ylim(0,7.5)
p.fert.mm
##### Ordination #######
library(vegan)
library(ggplot2)
   ##### Enzyme fingerprint #####
mat.ord <- read.csv('MAT.enz.ord.csv')
("Control" = "#2b83ba", "Fertilized" = "#abdda4", "Warmed" = "#d73027", "Warmed & Fertilized" = "#ffffbf")
mat.ord$col <- rep(c("#2b83ba","#d73027","#ffffbf", "#abdda4"), each = 4)
mat.ord$sum <- 4*(mat.summary$LAP.mm2.min+mat.summary$NAG.mm2.min+mat.summary$PHOS.mm2.min+mat.summary$ABTS.mm2.min+mat.summary$TMB.mm2.min)/
  max((mat.summary$LAP.mm2.min+mat.summary$NAG.mm2.min+mat.summary$PHOS.mm2.min+mat.summary$ABTS.mm2.min+mat.summary$TMB.mm2.min))
mat.ord$Ammonium

ord.enz <- metaMDS(mat.ord[,20:24], autotransform = F,trymax = 200, k =2, distance = "bray")
ord <- metaMDS(mat.ord[,20:24], autotransform = F,trymax = 200, k =3, distance = "bray")
stressplot(ord.enz)
save(ord.enz, file = "MAT.ord.enz.RData" )
load("MAT.ord.enz.RData")
str(ord.enz)
mat.ord$NMDS1 <- ord.enz$points[,1]
mat.ord$NMDS2 <- ord.enz$points[,2]
env.scores.mat <- as.data.frame(scores(ord.enz.fit, display = "vectors"))
env.scores.mat$env.variables <- c("Ammonium(ug N/g soil)", "Nitrate (ug N/g soil)", "C:N")
ord.enz$species

original.dist<-vegdist(mat.ord[,20:24], method = "bray") 
r2<-numeric(2) 
stress<-numeric(2) 
nmds.scores<-scores(ord.enz) 
nmds.dist<-dist(nmds.scores[,2],method = "euclidean") 
r2[2]<-summary(lm(original.dist~nmds.dist))[[8]] 
nmds.dist<-dist(nmds.scores[,1],method = "euclidean") 
r2[1]<-summary(lm(original.dist~nmds.dist))[[8]]
r2

ggplot(data = mat.ord)+
  geom_point(aes(x = NMDS1, y = NMDS2, shape = 1, fill = Treatment, size = sum), colour = "black", pch = 21)+
  scale_fill_manual(values = c("Control" = "#2b83ba",  "Warmed" = "#d73027", "Warmed.Fertilized" = "#ffffbf","Fertilized" = "#abdda4"))+
  geom_segment(data = env.scores.mat,
               aes(x = 0, xend = 1*NMDS1, y = 0, yend = 1*NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + 
  geom_text(data = env.scores.mat, 
            aes(x = 1*NMDS1, y = 1*NMDS2, label=env.variables),
            size = 4, 
            hjust = 0.8,
            vjust = 1.5)+
  ylim(-0.7, 1)+
  xlim(-0.7,1)+
  scale_size(range = 2.5*c(0.8198643, 4))+
  guides(fill = guide_legend(override.aes = list(size = 5)), size = FALSE)+
  labs(x = expression(paste("NMDS1, ",r^{2}, "= 0.78")), y = expression(paste("NMDS2, ",r^{2}, "= 0.08")))+
  theme(panel.border= element_rect(fill = NA, colour = "black"),
        panel.background = element_blank(),
        legend.position = c(0.84,0.85),
        legend.background = element_rect(fill = NA, colour = "black", size = 0.2),
        legend.key = element_rect(fill = NA))
  
plot(ord.enz, type = "n",,
     xlim = c(-0.7, 0.7), ylim = c(-0.7,0.7))
points(ord.enz, display = "sites", col = 1, pch = 21, bg = mat.ord$col, cex = 5*mat.ord$sum)
text(ord.enz, display = "species")
ord.enz.fit <- envfit(ord.enz~
                        aver.NH4.ugNgsoil+aver.NO3.ugNgsoil+aver.CN, na.rm=T, data = mat.ord, perm = 999)
plot(ord.enz.fit, col = "black", labels = c("Ammonium (ug N/g soil)", "Nitrate (ug N/g soil)", "C:N"))

perman.ord <- adonis(vegdist(mat.ord[,20:24], method = "bray") ~ mat.ord$Warmed*mat.ord$Fertilized#+
                       #aver.NH4.ugNgsoil+aver.NO3.ugNgsoil+aver.CN
                     ,
                     data = mat.ord)
perman.ord <- adonis(vegdist(mat.ord[,20:24], method = "bray") ~ , data = mat.ord)

beta.disp.ord <- betadisper(vegdist(mat.ord[,20:24]),group = mat.ord$Warmed:mat.ord$Fertilized)
TukeyHSD(beta.disp.ord)
save(perman.ord, beta.disp.ord, file = "MAT.enz.perMANOVA.RData")



library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
pairwise.adonis(original.dist,list(mat.ord$Warmed,mat.ord$Fertilized))
summary(simper(mat.ord[,20:24],mat.ord$Warmed:mat.ord$Fertilized))

