library (lme4)
library(emmeans)
library(e1071)
library(dplyr)
library(stats)


####GRANULOCYTES STATS ####
#outlier test and removal and skewness
########
# Calculate the quartiles and IQR
Q1 <- quantile(FCM_immuneaging$freq_granulocyte, 0.25)
Q3 <- quantile(FCM_immuneaging$freq_granulocyte, 0.75)
IQR <- Q3 - Q1

# Define the lower and upper bounds for outliers
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Identify outliers
outliers <- FCM_immuneaging$freq_granulocyte < lower_bound | FCM_immuneaging$freq_granulocyte > upper_bound

# Count the number of outliers
num_outliers <- sum(outliers)
print(paste("Number of outliers:", num_outliers))
#one outlier 

#remove the outlier 
# Remove outliers from the data
clean_data_FCM<- FCM_immuneaging[!outliers, ]

skew1 <- skewness(clean_data_FCM$freq_granulocyte)
print(skew1) #0.06327293 skewed positive
###########
#Run an ANOVA 
##########
clean_data_FCM$age <- factor(clean_data_FCM$age)
clean_data_FCM$sex <- factor(clean_data_FCM$sex)
# Fit ANOVA model for frequency of granulocytes
anova_gran <- aov(freq_granulocyte ~ age + sex + age:sex, data = clean_data_FCM)

# Print summary of the ANOVA
summary(anova_gran)

# If ANOVA shows significant differences among age groups, conduct post-hoc tests
# Example with Tukey's HSD post-hoc test
TukeyHSD(anova_gran)


####LYMPHOCYTES STATS ####
#outlier test and removal and skewness
########
# Calculate the quartiles and IQR
Q1 <- quantile(FCM_immuneaging$freq_lymp, 0.25)
Q3 <- quantile(FCM_immuneaging$freq_lymp, 0.75)
IQR <- Q3 - Q1

# Define the lower and upper bounds for outliers
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Identify outliers
outliers <- FCM_immuneaging$freq_lymp < lower_bound | FCM_immuneaging$freq_lymp > upper_bound

# Count the number of outliers
num_outliers <- sum(outliers)
print(paste("Number of outliers:", num_outliers))
#one outlier 

#remove the outlier 
# Remove outliers from the data
clean_data_FCM<- FCM_immuneaging[!outliers, ]

skew2 <- skewness(clean_data_FCM$freq_lymp)
print(skew2) #0.6594122 skewed positive (-0.5 to 0.5 considered relativley normal)

#########
#Run ANOVA
#################
clean_data_FCM$age <- factor(clean_data_FCM$age)
clean_data_FCM$sex <- factor(clean_data_FCM$sex)
# Fit ANOVA model for frequency of granulocytes
anova_lymp <- aov(freq_lymp ~ age + sex + age:sex, data = clean_data_FCM)

# Print summary of the ANOVA
summary(anova_lymp)

# If ANOVA shows significant differences among age groups, conduct post-hoc tests
# Example with Tukey's HSD post-hoc test
TukeyHSD(anova_lymp)
#no differences 