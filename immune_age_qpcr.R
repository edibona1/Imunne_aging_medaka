##Stats for RTqPCR data (as delta ct values) for immune aging study
library (lme4)
library(emmeans)
library(dplyr)
library(ggplot2)
library(moments)
library(dunn.test)
#### Assumptions checked####

#both sexes together 
shapiro.test(stats_immuno_aging$n_dct_c3)
#C3 not normal 
shapiro.test(stats_immuno_aging$n_dct_c8)
#C8 not normal
shapiro.test(stats_immuno_aging$n_dct_creact)
#CRP not normal
shapiro.test(stats_immuno_aging$n_dct_tlr5s)
#TLR5s is normal 
shapiro.test(stats_immuno_aging$n_dct_tlr5m)
#TLR5m is normal 
shapiro.test(stats_immuno_aging$n_dct_lyzc)
#LYZ not normal
shapiro.test(stats_immuno_aging$n_dct_myd88)
#MYD88 not normal
shapiro.test(stats_immuno_aging$n_dct_nfkb)
#NFKb not normal 
shapiro.test(stats_immuno_aging$n_dct_tcrb)
#TCRb not normal 
shapiro.test(stats_immuno_aging$n_dct_il1b)
#IL1b not normal
shapiro.test(stats_immuno_aging$n_dct_mhcii)
#MHCII not normal

####subset data #####
#subset female and male 
male_data <- subset(stats_immuno_aging, sex == "M")
female_data <- subset(stats_immuno_aging, sex == "F")
mo3_data <- subset(stats_immuno_aging, age== "3")
mo4_data <- subset(stats_immuno_aging, age== "4")
mo5_data <- subset(stats_immuno_aging, age== "5")
mo6_data <- subset(stats_immuno_aging, age== "6")
mo7_data <- subset(stats_immuno_aging, age== "7")
mo11_data <- subset(stats_immuno_aging, age== "11")
mo13_data <- subset(stats_immuno_aging, age== "13")
mo20_data <- subset(stats_immuno_aging, age== c("20", "23"))

#############differences based on sex within each age#############

##3mo normality check and ANOVA/KW ------ TLR5M difference####
shapiro.test(mo3_data$n_dct_c3)
#C3 is normal 
mo3_c3 <- aov(n_dct_c3 ~ sex, data = mo3_data)
summary(mo3_c3)
TukeyHSD(mo3_c3) #0.7045978
kruskal.test(n_dct_c3 ~ sex, data = mo3_data) #0.5127

shapiro.test(mo3_data$n_dct_c8)
#C8 is normal
mo3_c8 <- aov(n_dct_c8 ~ sex, data = mo3_data)
summary(mo3_c8)
TukeyHSD(mo3_c8) #0.6326547

shapiro.test(mo3_data$n_dct_creact)
#CRP is normal
mo3_crp <- aov(n_dct_creact ~ sex, data = mo3_data)
summary(mo3_crp)
TukeyHSD(mo3_crp) #0.4052164

shapiro.test(mo3_data$n_dct_tlr5s)
#TLR5s is normal 
mo3_tlr5s <- aov(n_dct_tlr5s ~ sex, data = mo3_data)
summary(mo3_tlr5s)
TukeyHSD(mo3_tlr5s) #0.5839716

shapiro.test(mo3_data$n_dct_tlr5m)
#TLR5m is normal 
mo3_tlr5m <- aov(n_dct_tlr5m ~ sex, data = mo3_data)
summary(mo3_tlr5m)
TukeyHSD(mo3_tlr5m) #0.0075006

shapiro.test(mo3_data$n_dct_lyzc)
#LYZ is normal
mo3_lyzc <- aov(n_dct_lyzc ~ sex, data = mo3_data)
summary(mo3_lyzc)
TukeyHSD(mo3_lyzc) #0.2066933

shapiro.test(mo3_data$n_dct_myd88)
#MYD88 is normal
mo3_myd88 <- aov(n_dct_myd88 ~ sex, data = mo3_data)
summary(mo3_myd88)
TukeyHSD(mo3_myd88) #0.9760482

shapiro.test(mo3_data$n_dct_nfkb)
#NFKb is normal 
mo3_nfkb <- aov(n_dct_nfkb ~ sex, data = mo3_data)
summary(mo3_nfkb)
TukeyHSD(mo3_nfkb) #0.2077675

shapiro.test(mo3_data$n_dct_tcrb)
#TCRb is normal 
mo3_tcrb <- aov(n_dct_tcrb ~ sex, data = mo3_data)
summary(mo3_tcrb)
TukeyHSD(mo3_tcrb) #0.1906069

shapiro.test(mo3_data$n_dct_il1b)
#IL1b is normal
mo3_il1b <- aov(n_dct_il1b ~ sex, data = mo3_data)
summary(mo3_il1b)
TukeyHSD(mo3_il1b) #0.8769968

shapiro.test(mo3_data$n_dct_mhcii)
#MHCII is normal
mo3_mhcii <- aov(n_dct_mhcii ~ sex, data = mo3_data)
summary(mo3_mhcii)
TukeyHSD(mo3_mhcii) #0.9070596

##4mo normality check and ANOVA/KW ------ TLR5S difference ####
shapiro.test(mo4_data$n_dct_c3)
#C3 is normal 
mo4_c3 <- aov(n_dct_c3 ~ sex, data = mo4_data)
summary(mo4_c3)
TukeyHSD(mo4_c3) #0.205282

shapiro.test(mo4_data$n_dct_c8)
#C8 is normal
mo4_c8 <- aov(n_dct_c8 ~ sex, data = mo4_data)
summary(mo4_c8)
TukeyHSD(mo4_c8) #0.3843548

shapiro.test(mo4_data$n_dct_creact)
#CRP is not normal
kruskal.test(n_dct_creact ~ sex, data = mo4_data) #0.7237

shapiro.test(mo4_data$n_dct_tlr5s)
#TLR5s is normal 
mo4_tlr5s <- aov(n_dct_tlr5s ~ sex, data = mo4_data)
summary(mo4_tlr5s)
TukeyHSD(mo4_tlr5s) #0.0416405

shapiro.test(mo4_data$n_dct_tlr5m)
#TLR5m is normal 
mo4_tlr5m <- aov(n_dct_tlr5m ~ sex, data = mo4_data)
summary(mo4_tlr5m)
TukeyHSD(mo4_tlr5m) #0.1479417

shapiro.test(mo4_data$n_dct_lyzc)
#LYZ is not normal
kruskal.test(n_dct_lyzc ~ sex, data = mo4_data) #1

shapiro.test(mo4_data$n_dct_myd88)
#MYD88 is normal
mo4_myd88 <- aov(n_dct_myd88 ~ sex, data = mo4_data)
summary(mo4_myd88)
TukeyHSD(mo4_myd88) #0.4125937

shapiro.test(mo4_data$n_dct_nfkb)
#NFKb is not normal 
kruskal.test(n_dct_nfkb ~ sex, data = mo4_data) #0.0771

shapiro.test(mo4_data$n_dct_tcrb)
#TCRb is normal 
mo4_tcrb <- aov(n_dct_tcrb ~ sex, data = mo4_data)
summary(mo4_tcrb)
TukeyHSD(mo4_tcrb) #0.9373318

shapiro.test(mo4_data$n_dct_il1b)
#IL1b is normal
mo4_il1b <- aov(n_dct_il1b ~ sex, data = mo4_data)
summary(mo4_il1b)
TukeyHSD(mo4_il1b) #0.432765

shapiro.test(mo4_data$n_dct_mhcii)
#MHCII is normal
mo4_mhcii <- aov(n_dct_mhcii ~ sex, data = mo4_data)
summary(mo4_mhcii)
TukeyHSD(mo4_mhcii) #0.4563179

##5mo normality check and ANOVA/KW ------ TLR5M LYZC difference####
shapiro.test(mo5_data$n_dct_c3)
#C3 is normal 
mo5_c3 <- aov(n_dct_c3 ~ sex, data = mo5_data)
summary(mo5_c3)
TukeyHSD(mo5_c3) #0.9453199

shapiro.test(mo5_data$n_dct_c8)
#C8 is normal
mo5_c8 <- aov(n_dct_c8 ~ sex, data = mo5_data)
summary(mo5_c8)
TukeyHSD(mo5_c8) #0.655671

shapiro.test(mo5_data$n_dct_creact)
#CRP is normal
mo5_creact <- aov(n_dct_creact ~ sex, data = mo5_data)
summary(mo5_creact)
TukeyHSD(mo5_creact) #0.7389883

shapiro.test(mo5_data$n_dct_tlr5s)
#TLR5s is normal 
mo5_tlr5s <- aov(n_dct_tlr5s ~ sex, data = mo5_data)
summary(mo5_tlr5s)
TukeyHSD(mo5_tlr5s) #0.1717852

shapiro.test(mo5_data$n_dct_tlr5m)
#TLR5m is normal 
mo5_tlr5m <- aov(n_dct_tlr5m ~ sex, data = mo5_data)
summary(mo5_tlr5m)
TukeyHSD(mo5_tlr5m) #0.0284941

shapiro.test(mo5_data$n_dct_lyzc)
#LYZ is normal
mo5_lyzc <- aov(n_dct_lyzc ~ sex, data = mo5_data)
summary(mo5_lyzc)
TukeyHSD(mo5_lyzc) #0.0226444

shapiro.test(mo5_data$n_dct_myd88)
#MYD88 is not normal
kruskal.test(n_dct_myd88 ~ sex, data = mo5_data) #1

shapiro.test(mo5_data$n_dct_nfkb)
#NFKb is normal 
mo5_nfkb <- aov(n_dct_nfkb ~ sex, data = mo5_data)
summary(mo5_nfkb)
TukeyHSD(mo5_nfkb) #0.6536827

shapiro.test(mo5_data$n_dct_tcrb)
#TCRb is normal 
mo5_tcrb <- aov(n_dct_tcrb ~ sex, data = mo5_data)
summary(mo5_tcrb)
TukeyHSD(mo5_tcrb) #0.1265109

shapiro.test(mo5_data$n_dct_il1b)
#IL1b is normal
mo5_il1b <- aov(n_dct_il1b ~ sex, data = mo5_data)
summary(mo5_il1b)
TukeyHSD(mo5_il1b) # 0.5930349

shapiro.test(mo5_data$n_dct_mhcii)
#MHCII is not normal
kruskal.test(n_dct_mhcii ~ sex, data = mo5_data) #0.7237

##6mo normality check and ANOVA/KW ------ NO differences####
shapiro.test(mo6_data$n_dct_c3)
#C3 is normal 
mo6_c3 <- aov(n_dct_c3 ~ sex, data = mo6_data)
summary(mo6_c3)
TukeyHSD(mo6_c3) #0.3240474

shapiro.test(mo6_data$n_dct_c8)
#C8 is normal
mo6_c8 <- aov(n_dct_c8 ~ sex, data = mo6_data)
summary(mo6_c8)
TukeyHSD(mo6_c8) #0.3754836

shapiro.test(mo6_data$n_dct_creact)
#CRP is normal
mo6_creact <- aov(n_dct_creact ~ sex, data = mo6_data)
summary(mo6_creact)
TukeyHSD(mo6_creact) #0.7106695

shapiro.test(mo6_data$n_dct_tlr5s)
#TLR5s is normal 
mo6_tlr5s <- aov(n_dct_tlr5s ~ sex, data = mo6_data)
summary(mo6_tlr5s)
TukeyHSD(mo6_tlr5s) #0.9310296

shapiro.test(mo6_data$n_dct_tlr5m)
#TLR5m is normal 
mo6_tlr5m <- aov(n_dct_tlr5m ~ sex, data = mo6_data)
summary(mo6_tlr5m)
TukeyHSD(mo6_tlr5m) #0.1472619

shapiro.test(mo6_data$n_dct_lyzc)
#LYZ is normal
mo6_lyzc <- aov(n_dct_lyzc ~ sex, data = mo6_data)
summary(mo6_lyzc)
TukeyHSD(mo6_lyzc) #0.4510335

shapiro.test(mo6_data$n_dct_myd88)
#MYD88 is normal
mo6_myd88 <- aov(n_dct_myd88 ~ sex, data = mo6_data)
summary(mo6_myd88)
TukeyHSD(mo6_myd88) #0.2766851

shapiro.test(mo6_data$n_dct_nfkb)
#NFKb is normal 
mo6_nfkb <- aov(n_dct_nfkb ~ sex, data = mo6_data)
summary(mo6_nfkb)
TukeyHSD(mo6_nfkb) #0.2076655

shapiro.test(mo6_data$n_dct_tcrb)
#TCRb is normal 
mo6_tcrb <- aov(n_dct_tcrb ~ sex, data = mo6_data)
summary(mo6_tcrb)
TukeyHSD(mo6_tcrb) #0.6804078

shapiro.test(mo6_data$n_dct_il1b)
#IL1b is normal
mo6_il1b <- aov(n_dct_il1b ~ sex, data = mo6_data)
summary(mo6_il1b)
TukeyHSD(mo6_il1b) #0.592892

shapiro.test(mo6_data$n_dct_mhcii)
#MHCII is normal
mo6_mhcii <- aov(n_dct_mhcii ~ sex, data = mo6_data)
summary(mo6_mhcii)
TukeyHSD(mo6_mhcii) #0.5310951

##7mo normality check and ANOVA/KW ------ TLR5S MYD88 difference####
shapiro.test(mo7_data$n_dct_c3)
#C3 is normal 
mo7_c3 <- aov(n_dct_c3 ~ sex, data = mo7_data)
summary(mo7_c3)
TukeyHSD(mo7_c3) #0.5562437

shapiro.test(mo7_data$n_dct_c8)
#C8 is normal
mo7_c8 <- aov(n_dct_c8 ~ sex, data = mo7_data)
summary(mo7_c8)
TukeyHSD(mo7_c8) #0.220753

shapiro.test(mo7_data$n_dct_creact)
#CRP is normal
mo7_creact <- aov(n_dct_creact ~ sex, data = mo7_data)
summary(mo7_creact)
TukeyHSD(mo7_creact) #0.801229

shapiro.test(mo7_data$n_dct_tlr5s)
#TLR5s is normal 
mo7_tlr5s <- aov(n_dct_tlr5s ~ sex, data = mo7_data)
summary(mo7_tlr5s)
TukeyHSD(mo7_tlr5s) #0.0039265

shapiro.test(mo7_data$n_dct_tlr5m)
#TLR5m is normal 
mo7_tlr5m <- aov(n_dct_tlr5m ~ sex, data = mo7_data)
summary(mo7_tlr5m)
TukeyHSD(mo7_tlr5m) #0.7542527

shapiro.test(mo7_data$n_dct_lyzc)
#LYZ is normal
mo7_lyzc <- aov(n_dct_lyzc ~ sex, data = mo7_data)
summary(mo7_lyzc)
TukeyHSD(mo7_lyzc) #0.7204956

shapiro.test(mo7_data$n_dct_myd88)
#MYD88 is normal
mo7_myd88 <- aov(n_dct_myd88 ~ sex, data = mo7_data)
summary(mo7_myd88)
TukeyHSD(mo7_myd88) #0.0287877

shapiro.test(mo7_data$n_dct_nfkb)
#NFKb is normal 
mo7_nfkb <- aov(n_dct_nfkb ~ sex, data = mo7_data)
summary(mo7_nfkb)
TukeyHSD(mo7_nfkb) #0.4866381

shapiro.test(mo7_data$n_dct_tcrb)
#TCRb is not normal 
kruskal.test(n_dct_tcrb ~ sex, data = mo7_data) #0.9168

shapiro.test(mo7_data$n_dct_il1b)
#IL1b is not normal
kruskal.test(n_dct_il1b ~ sex, data = mo7_data) #0.6015

shapiro.test(mo7_data$n_dct_mhcii)
#MHCII is normal
mo7_mhcii <- aov(n_dct_mhcii ~ sex, data = mo7_data)
summary(mo7_mhcii)
TukeyHSD(mo7_mhcii) #0.1967557

##11mo normality check and ANOVA/KW ------ C3 C8 difference####
shapiro.test(mo11_data$n_dct_c3)
#C3 is normal 
mo11_c3 <- aov(n_dct_c3 ~ sex, data = mo11_data)
summary(mo11_c3)
TukeyHSD(mo11_c3) #0.0306487

shapiro.test(mo11_data$n_dct_c8)
#C8 is normal
mo11_c8 <- aov(n_dct_c8 ~ sex, data = mo11_data)
summary(mo11_c8)
TukeyHSD(mo11_c8) #0.0097883

shapiro.test(mo11_data$n_dct_creact)
#CRP is  normal
mo11_creact <- aov(n_dct_creact ~ sex, data = mo11_data)
summary(mo11_creact)
TukeyHSD(mo11_creact) #0.0568824

shapiro.test(mo11_data$n_dct_tlr5s)
#TLR5s is normal 
mo11_tlr5s <- aov(n_dct_tlr5s ~ sex, data = mo11_data)
summary(mo11_tlr5s)
TukeyHSD(mo11_tlr5s) #0.3096492

shapiro.test(mo11_data$n_dct_tlr5m)
#TLR5m is normal 
mo11_tlr5m <- aov(n_dct_tlr5m ~ sex, data = mo11_data)
summary(mo11_tlr5m)
TukeyHSD(mo11_tlr5m) #0.2156222

shapiro.test(mo11_data$n_dct_lyzc)
#LYZ is not normal
kruskal.test(n_dct_lyzc ~ sex, data = mo11_data) #0.2087

shapiro.test(mo11_data$n_dct_myd88)
#MYD88 is normal
mo11_myd88 <- aov(n_dct_myd88 ~ sex, data = mo11_data)
summary(mo11_myd88)
TukeyHSD(mo11_myd88) # 0.212561

shapiro.test(mo11_data$n_dct_nfkb)
#NFKb is normal 
mo11_nfkb <- aov(n_dct_nfkb ~ sex, data = mo11_data)
summary(mo11_nfkb)
TukeyHSD(mo11_nfkb) #0.0535811

shapiro.test(mo11_data$n_dct_tcrb)
#TCRb is not normal 
kruskal.test(n_dct_tcrb ~ sex, data = mo11_data) #0.4647

shapiro.test(mo11_data$n_dct_il1b)
#IL1b is normal
mo11_il1b <- aov(n_dct_il1b ~ sex, data = mo11_data)
summary(mo11_il1b)
TukeyHSD(mo11_il1b) #0.0626344

shapiro.test(mo11_data$n_dct_mhcii)
#MHCII is normal
mo11_mhcii <- aov(n_dct_mhcii ~ sex, data = mo11_data)
summary(mo11_mhcii)
TukeyHSD(mo11_mhcii) #0.2248476

##13mo normality check and ANOVA/KW ------ C3 MYD88 TCRb MHCII difference####
shapiro.test(mo13_data$n_dct_c3)
#C3 is normal 
mo13_c3 <- aov(n_dct_c3 ~ sex, data = mo13_data)
summary(mo13_c3)
TukeyHSD(mo13_c3) #0.0461172

shapiro.test(mo13_data$n_dct_c8)
#C8 is not normal
kruskal.test(n_dct_c8 ~ sex, data = mo13_data) # 0.1172

shapiro.test(mo13_data$n_dct_creact)
#CRP is normal
mo13_creact <- aov(n_dct_creact ~ sex, data = mo13_data)
summary(mo13_creact)
TukeyHSD(mo13_creact) #0.081681

shapiro.test(mo13_data$n_dct_tlr5s)
#TLR5s is normal 
mo13_tlr5s <- aov(n_dct_tlr5s ~ sex, data = mo13_data)
summary(mo13_tlr5s)
TukeyHSD(mo13_tlr5s) #0.0651129

shapiro.test(mo13_data$n_dct_tlr5m)
#TLR5m is normal 
mo13_tlr5m <- aov(n_dct_tlr5m ~ sex, data = mo13_data)
summary(mo13_tlr5m)
TukeyHSD(mo13_tlr5m) # 0.4096812

shapiro.test(mo13_data$n_dct_lyzc)
#LYZ is normal
mo13_lyzc <- aov(n_dct_lyzc ~ sex, data = mo13_data)
summary(mo13_lyzc)
TukeyHSD(mo13_lyzc) #0.068046

shapiro.test(mo13_data$n_dct_myd88)
#MYD88 is normal
mo13_myd88 <- aov(n_dct_myd88 ~ sex, data = mo13_data)
summary(mo13_myd88)
TukeyHSD(mo13_myd88) #0.0018817

shapiro.test(mo13_data$n_dct_nfkb)
#NFKb is normal 
mo13_nfkb <- aov(n_dct_nfkb ~ sex, data = mo13_data)
summary(mo13_nfkb)
TukeyHSD(mo13_nfkb) #0.0562601

shapiro.test(mo13_data$n_dct_tcrb)
#TCRb is normal 
mo13_tcrb <- aov(n_dct_tcrb ~ sex, data = mo13_data)
summary(mo13_tcrb)
TukeyHSD(mo13_tcrb) #0.0414957

shapiro.test(mo13_data$n_dct_il1b)
#IL1b is not normal
kruskal.test(n_dct_il1b ~ sex, data = mo13_data) # 0.0758

shapiro.test(mo13_data$n_dct_mhcii)
#MHCII is normal
mo13_mhcii <- aov(n_dct_mhcii ~ sex, data = mo13_data)
summary(mo13_mhcii)
TukeyHSD(mo13_mhcii) #0.0136381

##20/23mo normality check and ANOVA/KW ------ MYD88 difference####
shapiro.test(mo20_data$n_dct_c3)
#C3 is normal 
mo20_c3 <- aov(n_dct_c3 ~ sex, data = mo20_data)
summary(mo20_c3)
TukeyHSD(mo20_c3) #0.3403747

shapiro.test(mo20_data$n_dct_c8)
#C8 is normal
mo20_c8 <- aov(n_dct_c8 ~ sex, data = mo20_data)
summary(mo20_c8)
TukeyHSD(mo20_c8) #0.2647842

shapiro.test(mo20_data$n_dct_creact)
#CRP is not normal
kruskal.test(n_dct_creact ~ sex, data = mo20_data) #0.4386

shapiro.test(mo20_data$n_dct_tlr5s)
#TLR5s is not normal 
kruskal.test(n_dct_tlr5s ~ sex, data = mo20_data) #0.1213

shapiro.test(mo20_data$n_dct_tlr5m)
#TLR5m is normal 
mo20_tlr5m <- aov(n_dct_tlr5m ~ sex, data = mo20_data)
summary(mo20_tlr5m)
TukeyHSD(mo20_tlr5m) #0.1269645

shapiro.test(mo20_data$n_dct_lyzc)
#LYZ is normal
mo20_lyzc <- aov(n_dct_lyzc ~ sex, data = mo20_data)
summary(mo20_lyzc)
TukeyHSD(mo20_lyzc) #0.677358

shapiro.test(mo20_data$n_dct_myd88)
#MYD88 is normal
mo20_myd88 <- aov(n_dct_myd88 ~ sex, data = mo20_data)
summary(mo20_myd88)
TukeyHSD(mo20_myd88) #0.025438

shapiro.test(mo20_data$n_dct_nfkb)
#NFKb is not normal 
kruskal.test(n_dct_nfkb ~ sex, data = mo20_data) #0.1213

shapiro.test(mo20_data$n_dct_tcrb)
#TCRb is normal 
mo20_tcrb <- aov(n_dct_tcrb ~ sex, data = mo20_data)
summary(mo20_tcrb)
TukeyHSD(mo20_tcrb) #0.1172404

shapiro.test(mo20_data$n_dct_il1b)
#IL1b is normal
mo20_il1b <- aov(n_dct_il1b ~ sex, data = mo20_data)
summary(mo20_il1b)
TukeyHSD(mo20_il1b) #0.3722477

shapiro.test(mo20_data$n_dct_mhcii)
#MHCII is normal
mo20_mhcii <- aov(n_dct_mhcii ~ sex, data = mo20_data)
summary(mo20_mhcii)
TukeyHSD(mo20_mhcii) #0.2557472


###GRAPHS for sex difference ####
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggplot2)
library(dplyr)
#3mo graph####
# Assuming mo3_data is your dataframe for 3-month-old samples
summary_stats <- mo3_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups
gene_groups <- data.frame(
  gene = c("n_dct_creact", "n_dct_tlr5s", "n_dct_tlr5m", "n_dct_tcrb", "n_dct_mhcii", 
           "n_dct_myd88", "n_dct_nfkb", "n_dct_c3", "n_dct_il1b", "n_dct_c8", "n_dct_lyzc"),
  group = c("immune initiators", "immune initiators", "immune initiators", "immune initiators", "immune initiators",
            "mediators", "mediators", "mediators", "mediators", "effectors", "effectors")
)

# Add the group information to long_data
long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = gene_name_changes))

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune Genes", y = "-∆ Ct", title = "3 mph") +
  theme_bw() +
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("blue", "pink")) + # Customize colors if needed
  facet_wrap(~ group, scales = "free_x", ncol = 1) # Facet by gene group

####new 3mo graph with facets ####
# Assuming mo3_data is your dataframe for 3-month-old samples
summary_stats <- mo3_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))
View(summary_stats)

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups using updated gene names
gene_groups <- data.frame(
  gene = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
           "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ"),
  group = c("Initiators", "Initiators", "Initiators", "Initiators", "Initiators",
            "Mediators", "Mediators", "Mediators", "Mediators", "Effectors", "Effectors")
)

long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
                                        "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ")),
         group = factor(group, levels = c("Initiators", "Mediators", "Effectors")))  # Set the order of facets

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune genes", y = "-∆Ct", title = "3 mph") +
  theme_bw(base_size = 20) +  # Set base font size for theme_bw
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 20)  # Facet label text size
  ) +
  scale_fill_manual(values = c("blue", "pink")) +  # Customize colors if needed
  scale_x_discrete(labels = function(x) {
    # Display gene names as they are, without italicizing
    x
  }) +
  facet_grid(~ group, scales = "free_x")  # Facet by gene group

#4mo graph####
# Assuming mo4_data is your dataframe for 3-month-old samples
summary_stats <- mo4_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups
gene_groups <- data.frame(
  gene = c("n_dct_creact", "n_dct_tlr5s", "n_dct_tlr5m", "n_dct_tcrb", "n_dct_mhcii", 
           "n_dct_myd88", "n_dct_nfkb", "n_dct_c3", "n_dct_il1b", "n_dct_c8", "n_dct_lyzc"),
  group = c("immune initiators", "immune initiators", "immune initiators", "immune initiators", "immune initiators",
            "mediators", "mediators", "mediators", "mediators", "effectors", "effectors")
)

# Add the group information to long_data
long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = gene_name_changes))

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune Genes", y = "-∆ Ct", title = "4 mph") +
  theme_bw() +
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("blue", "pink")) + # Customize colors if needed
  facet_wrap(~group, ncol = 3) # Facet by gene group


####new 4mo graph with facets ####
# Assuming mo4_data is your dataframe for 4-month-old samples
summary_stats <- mo4_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))
View(summary_stats)

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups using updated gene names
gene_groups <- data.frame(
  gene = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
           "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ"),
  group = c("Initiators", "Initiators", "Initiators", "Initiators", "Initiators",
            "Mediators", "Mediators", "Mediators", "Mediators", "Effectors", "Effectors")
)

long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
                                        "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ")),
         group = factor(group, levels = c("Initiators", "Mediators", "Effectors")))  # Set the order of facets

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune genes", y = "-∆Ct", title = "4 mph") +
  theme_bw(base_size = 20) +  # Set base font size for theme_bw
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 20)  # Facet label text size
  ) +
  scale_fill_manual(values = c("blue", "pink")) +  # Customize colors if needed
  scale_x_discrete(labels = function(x) {
    # Display gene names as they are, without italicizing
    x
  }) +
  facet_grid(~ group, scales = "free_x")  # Facet by gene group
#5mo graph####
# Assuming mo5_data is your dataframe for 3-month-old samples
summary_stats <- mo5_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups
gene_groups <- data.frame(
  gene = c("n_dct_creact", "n_dct_tlr5s", "n_dct_tlr5m", "n_dct_tcrb", "n_dct_mhcii", 
           "n_dct_myd88", "n_dct_nfkb", "n_dct_c3", "n_dct_il1b", "n_dct_c8", "n_dct_lyzc"),
  group = c("immune initiators", "immune initiators", "immune initiators", "immune initiators", "immune initiators",
            "mediators", "mediators", "mediators", "mediators", "effectors", "effectors")
)

# Add the group information to long_data
long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = gene_name_changes))

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune Genes", y = "-∆ Ct", title = "5 mph") +
  theme_bw() +
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("blue", "pink")) + # Customize colors if needed
  facet_wrap(~group, ncol = 3) # Facet by gene group
####new 5mo graph with facets ####
# Assuming mo5_data is your dataframe for 5-month-old samples
summary_stats <- mo5_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))
View(summary_stats)

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups using updated gene names
gene_groups <- data.frame(
  gene = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
           "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ"),
  group = c("Initiators", "Initiators", "Initiators", "Initiators", "Initiators",
            "Mediators", "Mediators", "Mediators", "Mediators", "Effectors", "Effectors")
)

long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
                                        "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ")),
         group = factor(group, levels = c("Initiators", "Mediators", "Effectors")))  # Set the order of facets

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune genes", y = "-∆Ct", title = "5 mph") +
  theme_bw(base_size = 20) +  # Set base font size for theme_bw
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 20)  # Facet label text size
  ) +
  scale_fill_manual(values = c("blue", "pink")) +  # Customize colors if needed
  scale_x_discrete(labels = function(x) {
    # Display gene names as they are, without italicizing
    x
  }) +
  facet_grid(~ group, scales = "free_x")  # Facet by gene group
#6mo graph####
# Assuming mo6_data is your dataframe for 3-month-old samples
summary_stats <- mo6_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups
gene_groups <- data.frame(
  gene = c("n_dct_creact", "n_dct_tlr5s", "n_dct_tlr5m", "n_dct_tcrb", "n_dct_mhcii", 
           "n_dct_myd88", "n_dct_nfkb", "n_dct_c3", "n_dct_il1b", "n_dct_c8", "n_dct_lyzc"),
  group = c("immune initiators", "immune initiators", "immune initiators", "immune initiators", "immune initiators",
            "mediators", "mediators", "mediators", "mediators", "effectors", "effectors")
)

# Add the group information to long_data
long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = gene_name_changes))

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune Genes", y = "-∆ Ct", title = "6 mph") +
  theme_bw() +
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("blue", "pink")) + # Customize colors if needed
  facet_wrap(~group, ncol = 3) # Facet by gene group

####new 6mo graph with facets ####
# Assuming mo6_data is your dataframe for 6-month-old samples
summary_stats <- mo6_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups using updated gene names
gene_groups <- data.frame(
  gene = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
           "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ"),
  group = c("Initiators", "Initiators", "Initiators", "Initiators", "Initiators",
            "Mediators", "Mediators", "Mediators", "Mediators", "Effectors", "Effectors")
)

long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
                                        "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ")),
         group = factor(group, levels = c("Initiators", "Mediators", "Effectors")))  # Set the order of facets

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune genes", y = "-∆Ct", title = "6 mph") +
  theme_bw(base_size = 20) +  # Set base font size for theme_bw
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 20)  # Facet label text size
  ) +
  scale_fill_manual(values = c("blue", "pink")) +  # Customize colors if needed
  scale_x_discrete(labels = function(x) {
    # Display gene names as they are, without italicizing
    x
  }) +
  #scale_y_continuous(limits = c(-10, 20)) + #set y axis scale
  facet_grid(~ group, scales = "free_x")  # Facet by gene group
#7mo graph####
# Assuming mo7_data is your dataframe for 3-month-old samples
summary_stats <- mo7_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups
gene_groups <- data.frame(
  gene = c("n_dct_creact", "n_dct_tlr5s", "n_dct_tlr5m", "n_dct_tcrb", "n_dct_mhcii", 
           "n_dct_myd88", "n_dct_nfkb", "n_dct_c3", "n_dct_il1b", "n_dct_c8", "n_dct_lyzc"),
  group = c("immune initiators", "immune initiators", "immune initiators", "immune initiators", "immune initiators",
            "mediators", "mediators", "mediators", "mediators", "effectors", "effectors")
)

# Add the group information to long_data
long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = gene_name_changes))

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune Genes", y = "-∆ Ct", title = "7 mph") +
  theme_bw() +
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("blue", "pink")) + # Customize colors if needed
  facet_wrap(~group, ncol = 3) # Facet by gene group

####new 7mo graph with facets ####
# Assuming mo7_data is your dataframe for 7-month-old samples
summary_stats <- mo7_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups using updated gene names
gene_groups <- data.frame(
  gene = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
           "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ"),
  group = c("Initiators", "Initiators", "Initiators", "Initiators", "Initiators",
            "Mediators", "Mediators", "Mediators", "Mediators", "Effectors", "Effectors")
)

long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
                                        "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ")),
         group = factor(group, levels = c("Initiators", "Mediators", "Effectors")))  # Set the order of facets

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune genes", y = "-∆Ct", title = "7 mph") +
  theme_bw(base_size = 20) +  # Set base font size for theme_bw
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 20)  # Facet label text size
  ) +
  scale_fill_manual(values = c("blue", "pink")) +  # Customize colors if needed
  scale_x_discrete(labels = function(x) {
    # Display gene names as they are, without italicizing
    x
  }) +
  #scale_y_continuous(limits = c(-10, 20)) + #set y axis scale
  facet_grid(~ group, scales = "free_x")  # Facet by gene group
#11mo graph####
# Assuming mo11_data is your dataframe for 3-month-old samples
summary_stats <- mo11_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups
gene_groups <- data.frame(
  gene = c("n_dct_creact", "n_dct_tlr5s", "n_dct_tlr5m", "n_dct_tcrb", "n_dct_mhcii", 
           "n_dct_myd88", "n_dct_nfkb", "n_dct_c3", "n_dct_il1b", "n_dct_c8", "n_dct_lyzc"),
  group = c("immune initiators", "immune initiators", "immune initiators", "immune initiators", "immune initiators",
            "mediators", "mediators", "mediators", "mediators", "effectors", "effectors")
)

# Add the group information to long_data
long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = gene_name_changes))

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune Genes", y = "-∆ Ct", title = "11 mph") +
  theme_bw() +
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("blue", "pink")) + # Customize colors if needed
  facet_wrap(~group, ncol = 3) # Facet by gene group
####new 11mo graph with facets ####
# Assuming mo11_data is your dataframe for 11-month-old samples
summary_stats <- mo11_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups using updated gene names
gene_groups <- data.frame(
  gene = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
           "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ"),
  group = c("Initiators", "Initiators", "Initiators", "Initiators", "Initiators",
            "Mediators", "Mediators", "Mediators", "Mediators", "Effectors", "Effectors")
)

long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
                                        "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ")),
         group = factor(group, levels = c("Initiators", "Mediators", "Effectors")))  # Set the order of facets

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune genes", y = "-∆Ct", title = "11 mph") +
  theme_bw(base_size = 20) +  # Set base font size for theme_bw
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 20)  # Facet label text size
  ) +
  scale_fill_manual(values = c("blue", "pink")) +  # Customize colors if needed
  scale_x_discrete(labels = function(x) {
    # Display gene names as they are, without italicizing
    x
  }) +
  #scale_y_continuous(limits = c(-10, 20)) + #set y axis scale
  facet_grid(~ group, scales = "free_x")  # Facet by gene group
#13mo graph####
####new 13mo graph with facets ####
# Assuming mo13_data is your dataframe for 13-month-old samples
summary_stats <- mo13_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups using updated gene names
gene_groups <- data.frame(
  gene = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
           "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ"),
  group = c("Initiators", "Initiators", "Initiators", "Initiators", "Initiators",
            "Mediators", "Mediators", "Mediators", "Mediators", "Effectors", "Effectors")
)

long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
                                        "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ")),
         group = factor(group, levels = c("Initiators", "Mediators", "Effectors")))  # Set the order of facets

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune genes", y = "-∆Ct", title = "13 mph") +
  theme_bw(base_size = 20) +  # Set base font size for theme_bw
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 20)  # Facet label text size
  ) +
  scale_fill_manual(values = c("blue", "pink")) +  # Customize colors if needed
  scale_x_discrete(labels = function(x) {
    # Display gene names as they are, without italicizing
    x
  }) +
  #scale_y_continuous(limits = c(-10, 20)) + #set y axis scale
  facet_grid(~ group, scales = "free_x")  # Facet by gene group
#20-23mo graph####
# Assuming mo20_data is your dataframe for 3-month-old samples
summary_stats <- mo20_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups
gene_groups <- data.frame(
  gene = c("n_dct_creact", "n_dct_tlr5s", "n_dct_tlr5m", "n_dct_tcrb", "n_dct_mhcii", 
           "n_dct_myd88", "n_dct_nfkb", "n_dct_c3", "n_dct_il1b", "n_dct_c8", "n_dct_lyzc"),
  group = c("immune initiators", "immune initiators", "immune initiators", "immune initiators", "immune initiators",
            "mediators", "mediators", "mediators", "mediators", "effectors", "effectors")
)

# Add the group information to long_data
long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = gene_name_changes))

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune Genes", y = "-∆ Ct", title = "20-23 mph") +
  theme_bw() +
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("blue", "pink")) + # Customize colors if needed
  facet_wrap(~group, ncol = 3) # Facet by gene group

####new 20-23mo graph with facets ####
# Assuming mo20_data is your dataframe for 20-23-month-old samples
summary_stats <- mo20_data %>%
  group_by(sex) %>%
  summarise(across(starts_with("n_dct_"), list(mean = mean, sd = sd)))

gene_columns <- grep("_mean$", colnames(summary_stats), value = TRUE)

# Create a long format data frame
long_data <- data.frame()
for (gene_col in gene_columns) {
  # Extract the gene name
  gene <- sub("_mean$", "", gene_col)
  # Create temporary data frame for this gene
  temp_df <- data.frame(
    sex = summary_stats$sex,
    gene = gene,
    mean = summary_stats[[gene_col]],
    sd = summary_stats[[paste0(gene, "_sd")]]
  )
  # Append to the long format data frame
  long_data <- rbind(long_data, temp_df)
}

# Define gene name changes
gene_name_changes <- c(
  "n_dct_creact" = "CRP",
  "n_dct_tlr5s" = "TLR5-s",
  "n_dct_tlr5m" = "TLR5-m",
  "n_dct_tcrb" = "TCRb",
  "n_dct_mhcii" = "MHCII",
  "n_dct_myd88" = "MYD88",
  "n_dct_nfkb" = "Nf-kß",
  "n_dct_c3" = "C3",
  "n_dct_il1b" = "IL1b",
  "n_dct_lyzc" = "LYZ",
  "n_dct_c8" = "C8"
)

# Replace gene names in long_data
long_data$gene <- factor(long_data$gene, levels = names(gene_name_changes))
levels(long_data$gene) <- gene_name_changes[levels(long_data$gene)]

# Define gene groups using updated gene names
gene_groups <- data.frame(
  gene = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
           "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ"),
  group = c("Initiators", "Initiators", "Initiators", "Initiators", "Initiators",
            "Mediators", "Mediators", "Mediators", "Mediators", "Effectors", "Effectors")
)

long_data <- long_data %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(gene = factor(gene, levels = c("CRP", "TLR5-s", "TLR5-m", "TCRb", "MHCII", 
                                        "MYD88", "Nf-kß", "C3", "IL1b", "C8", "LYZ")),
         group = factor(group, levels = c("Initiators", "Mediators", "Effectors")))  # Set the order of facets

# Create the plot with faceting
ggplot(long_data, aes(x = gene, y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7), width = 0.2) +
  labs(x = "Immune genes", y = "-∆Ct", title = "20-23 mph") +
  theme_bw(base_size = 20) +  # Set base font size for theme_bw
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 20)  # Facet label text size
  ) +
  scale_fill_manual(values = c("blue", "pink")) +  # Customize colors if needed
  scale_x_discrete(labels = function(x) {
    # Display gene names as they are, without italicizing
    x
  }) +
  #scale_y_continuous(limits = c(-10, 20)) + #set y axis scale
  facet_grid(~ group, scales = "free_x")  # Facet by gene group
#############differences based on age within sex#############
##Female normality and differences -------- ####
#female normality check
shapiro.test(female_data$n_dct_c3)
#C3 not normal 
kruskal.test(n_dct_c3 ~ age, data = female_data) #0.0014
c3_dunn<-dunn.test(female_data,n_dct_c3 ~ age, p.adjust.method = "BH", detailed = FALSE)
View(c3_dunn)
shapiro.test(female_data$n_dct_c8)
#C8 not normal
shapiro.test(female_data$n_dct_creact)
#CRP is normal
shapiro.test(female_data$n_dct_tlr5s)
#TLR5s is normal 
shapiro.test(female_data$n_dct_tlr5m)
#TLR5m is normal 
shapiro.test(female_data$n_dct_lyzc)
#LYZ not normal
shapiro.test(female_data$n_dct_myd88)
#MYD88 not normal
shapiro.test(female_data$n_dct_nfkb)
#NFKb is normal 
shapiro.test(female_data$n_dct_tcrb)
#TCRb is normal 
shapiro.test(female_data$n_dct_il1b)
#IL1b not normal
shapiro.test(female_data$n_dct_mhcii)
#MHCII is normal

##Male normality and differences -------- ####
#male normality check
shapiro.test(male_data$n_dct_c3)
#C3 is normal 
shapiro.test(male_data$n_dct_c8)
#C8 is normal
shapiro.test(male_data$n_dct_creact)
#CRP not normal
shapiro.test(male_data$n_dct_tlr5s)
#TLR5s is normal 
shapiro.test(male_data$n_dct_tlr5m)
#TLR5m is normal 
shapiro.test(male_data$n_dct_lyzc)
#LYZ not normal
shapiro.test(male_data$n_dct_myd88)
#MYD88 not normal
shapiro.test(male_data$n_dct_nfkb)
#NFKb not normal 
shapiro.test(male_data$n_dct_tcrb)
#TCRb not normal 
shapiro.test(male_data$n_dct_il1b)
#IL1b not normal
shapiro.test(male_data$n_dct_mhcii)
#MHCII not normal



####ANOVAS/KW####
male_data$age <- factor(male_data$age)
female_data$age <-factor(female_data$age)
male_data$sex <- factor(male_data$sex)
female_data$sex <-factor(female_data$sex)
stats_immuno_aging$age <-factor(stats_immuno_aging$age)
stats_immuno_aging$sex <- factor(stats_immuno_aging$sex)

kruskal.test(n_dct_c3 ~ sex, data = stats_immuno_aging)
#p-value = 0.6259
sex_anova<- aov(n_dct_c3 ~ age * sex, data=stats_immuno_aging)
summary(sex_anova)
TukeyHSD(sex_anova)
#no differences in expression within each age based on sex 


kruskal.test(n_dct_c3 ~ age, data = stats_immuno_aging)
#p-value = 1.882e-07
kruskal.test(n_dct_c3 ~ age, data = female_data)
#p-value = 0.0014
c3_dunn<-dunn_test(female_data,n_dct_c3 ~ age, p.adjust.method = "BH", detailed = FALSE)
View(c3_dunn)
male_data$age <- factor(male_data$age)
c3_m_aov<-aov(n_dct_c3 ~ age, data = male_data)
summary(c3_m_aov)
c3_m_tuk<-TukeyHSD(c3_m_aov)

kruskal.test(n_dct_tlr5m ~age, data=female_data)
tlr5m_dunn<-dunn_test(female_data,n_dct_tlr5m ~ age, p.adjust.method = "BH", detailed = FALSE)
View(tlr5m_dunn)
tlr5m_f_aov<-aov(n_dct_tlr5m ~ age, female_data)
summary(tlr5m_f_aov)
TukeyHSD(tlr5m_f_aov)


kruskal.test(n_dct_c8 ~ sex, data = stats_immuno_aging)
#p-value = 0.6259
kruskal.test(n_dct_c8 ~ age, data = stats_immuno_aging)
#p-value = 1.882e-07
kruskal.test(n_dct_c8 ~ age, data = female_data)
#p-value = 0.0014
c8_dunn<-dunn_test(female_data,n_dct_c8 ~ age, p.adjust.method = "BH", detailed = FALSE)
View(c8_dunn)









##effect size Cohen d test ####
library(dplyr)
library(effsize)
library(writexl)

# Define the list of genes
genes <- c("n_dct_creact", "n_dct_tlr5m", "n_dct_tlr5s", "n_dct_tcrb", 
           "n_dct_mhcii", "n_dct_myd88", "n_dct_nfkb", "n_dct_c3", 
           "n_dct_il1b", "n_dct_c8", "n_dct_lyzc")

# Initialize an empty data frame to store the results
cohen_d_results <- data.frame()

# Function to calculate Cohen's d between two age groups
calculate_cohen_d_pairwise <- function(data, gene, age1, age2) {
  group1 <- data %>% filter(age == age1) %>% pull(!!sym(gene))
  group2 <- data %>% filter(age == age2) %>% pull(!!sym(gene))
  
  if(length(group1) > 0 && length(group2) > 0) {
    # Calculate Cohen's d for the two groups
    cohen_d <- effsize::cohen.d(group1, group2)$estimate
    return(data.frame(gene = gene, age_comparison = paste(age1, "vs", age2), cohen_d = cohen_d))
  } else {
    return(data.frame(gene = gene, age_comparison = paste(age1, "vs", age2), cohen_d = NA))
  }
}

# Loop over each gene and calculate Cohen's d for each pair of ages
for (gene in genes) {
  # Get unique age groups
  age_groups <- unique(female_data$age)
  
  # Compare each pair of ages
  for (i in 1:(length(age_groups) - 1)) {
    for (j in (i + 1):length(age_groups)) {
      age1 <- age_groups[i]
      age2 <- age_groups[j]
      
      # Calculate Cohen's d for females
      female_cohen_d <- calculate_cohen_d_pairwise(female_data, gene, age1, age2)
      female_cohen_d$sex <- "Female"
      
      # Calculate Cohen's d for males
      male_cohen_d <- calculate_cohen_d_pairwise(male_data, gene, age1, age2)
      male_cohen_d$sex <- "Male"
      
      # Combine results
      cohen_d_results <- rbind(cohen_d_results, female_cohen_d, male_cohen_d)
    }
  }
}

# Print the results
print(cohen_d_results)
View(cohen_d_results)
write.csv(cohen_d_results, file = "cohend.csv", row.names = TRUE)
