#in this script I will use the tidyverse to run LMM on my mass loss and mineralization data
#the code will be reproducible
#it will use the tidy version of the combined data
install.packages("lsmeans")
install.packages("multcomp")


library(tidyverse)
library(ggplot2)
library(readxl)
library(lme4)
library(lmerTest)
library(effects)
library(arms)
library(lsmeans)
library(multcomp)


tidy_data <- read_xlsx("C:\\Users\\kelle\\OneDrive\\Thesis\\CNAnalysis\\tidy_mass_cn_icp_pctchange_data.xlsx")
View(tidy_data)

glimpse(tidy_data)

tidy_data <- tidy_data %>% 
  mutate(stand_age = case_when(
    stand %in% c("C9", "HBM") ~"mid",
    stand %in% c("C8", "HBO", "JBO") ~ "old",
    TRUE ~ NA_character_
  ))



tidy_data <- tidy_data %>% 
  mutate(destination = case_when(
    tolower(destination) %in% c("control") ~ "C",
    tolower(destination) %in% c("n+p") ~ "NP",
    tolower(destination) %in% c("0.87") ~ "C",
    TRUE ~ destination  # Keep other values unchanged
  )) %>% 
  mutate(fert_type = case_when(
    (source == "N" | destination == "N") ~ "N",
    (source == "P" | destination == "P") ~ "P",
    (source == "NP" | destination == "NP") ~ "NP",
    (source == "Ca" | destination == "Ca") ~ "Ca",
    (source == "C" & destination == "C") ~ "Control",
    TRUE ~ NA  # New condition for none of the specified conditions
  )) %>% 
  mutate(fert_location = case_when(
    destination %in% c("N", "P", "NP", "Ca") ~ "Soil",
    source %in% c("N", "P", "NP") ~ "Leaf",
    source == "C" & destination == "C" ~ "Control",
    TRUE ~ NA  # New condition for none of the specified conditions
  )) %>% 
  filter(!is.na(fert_type)) %>% 
  mutate(pct_element = (mgelement_gleaf / 10), #1000 mg in 1 gram so mg/1000, then multiply by 100 to make a percentage
         pct_c_pct_element_ratio_decomposed = pct_c / pct_element,
         initial_pct_c_pct_element_ratio = average_pctc / (initial_conc_icp/10),
         pct_change_c_element_ratio = ((pct_c_pct_element_ratio_decomposed - initial_pct_c_pct_element_ratio) / initial_pct_c_pct_element_ratio)) %>% 
  filter(id != "146")


unique(tidy_data$stand_age)

treated_source <- tidy_data %>% 
  filter(!(destination %in% c("N", "P", "NP", "Ca"))) %>% 
  filter(source != "?")

treated_destination <- tidy_data %>% 
  filter(!(source %in% c("N", "P", "NP"))) %>% 
  filter(destination != "Ca") %>% 
  filter(source != "?")


#confirm n and p columns are factors
class(tidy_data$p)

one_element <- treated_source %>% 
  filter(element == "K") 

Mg_negative_pct_change <- one_element %>% 
  filter(percent_change_icp < 0 & year == "2023")

View(Mg_negative_pct_change)
View(one_element)

################### treated source mineralization model #######################3333
source_model <- lmer(percent_change_icp ~ source_n*source_p + stand_age + (1|site/stand/destination) + (1|days) + (1|pct_water), data = one_element)

summary(source_model)

anova_result <- anova(model, type = 3) #main effects, check interaction 
print(anova_result)

shapiro.test(resid(model))

#log transform if necessary
source_model <- lmer(log(percent_change_icp) ~ source_n*source_p + stand_age + (1|site/stand/destination) + (1|days) + (1|pct_water), data = one_element)


simple_effects <- lsmeans(source_model, pairwise ~ source_p | source_n)
summary(simple_effects)
simple_effects <- lsmeans(source_model, pairwise ~ source_n | source_p)

### Interaction Plot ###
# Create an interaction plot
interaction_plot <- allEffects(model, xlevels = list(source_n = unique(one_element$source_n), source_p = unique(one_element$source_p)))
print(interaction_plot)
# Plot the interaction
plot(interaction_plot, multiline = TRUE, col = c("red", "blue"), lwd = 2, xlab = "source_n", ylab = "Percent Change ICP")

View(one_element)

################################################## treated destination model ###################################

one_element <- treated_destination %>% 
  filter( element == "K")

dest_model <- lmer(percent_change_icp ~ destination_n*destination_p + stand_age + (1|site/stand/destination) + (1|days) + (1|pct_water), data = one_element)

#log transformation if needed
dest_model <- lmer(log(percent_change_icp) ~ destination_n*destination_p + stand_age + (1|site/stand/destination) + (1|days) + (1|pct_water), data = one_element)


summary(dest_model)

anova_result <- anova(dest_model, type = 3) #main effects, check interaction 
print(anova_result)

shapiro.test(resid(dest_model))

######### simple effects
simple_effects <- lsmeans(dest_model, pairwise ~ destination_p | destination_n, ddf = "satterthwaite")
summary(simple_effects)
simple_effects <- lsmeans(dest_model, pairwise ~ destination_n | destination_p, ddf = "satterthwaite")
simple_effects <- lsmeans(dest_model, pairwise ~ stand_age, ddf = "satterthwaite")


########################################### comparing source and dest models
install.packages("MuMIn")
library(MuMIn)


r.squaredGLMM(source_model)
#R2m = proportion of total variance explained by the fixed effects - 5%
#R2c = proportion of total variance explained by the fixed and random effects together - 5%

r.squaredGLMM(dest_model)




##################################### Calculation Effect Size for each component of the anova #######################3333
########################## Omega Squared Value #######################################
###### doesn't work yet
install.packages("rstatix")

library(rstatix)
library(lme4)
library(purrr)

################getting means and variances for the Cohen's D calculation 

one_element <- treated_destination %>% 
  filter(element == "Mn")

one_element <- treated_destination %>% 
  distinct(id, .keep_all = TRUE) %>% 
  filter(!is.na(percent_change_cn)

View(one_element)

control <- one_element %>% 
  filter( destination == "C" & source == "C")

mean(control$percent_change_icp)
sqrt(var(control$percent_change_icp))
nrow(control$percent_change_icp)

View(control)

treatment <- one_element %>% 
  filter(destination == "P" | destination == "NP")

treatment <- one_element %>% 
  filter(destination == "N" | destination == "NP")

treatment <- one_element %>% 
  filter(destination == "P")

mean(treatment$percent_change_icp)
sqrt(var(treatment$percent_change_icp))
nrow(treatment$percent_change_icp)

View(treatment)
##################### ##############################treated source
one_element <- treated_source %>% 
  filter(element == "Mg")

one_element <- treated_source %>% 
  filter(!is.na(percent_change_icp))

View(one_element)

unique(one_element$id)

control <- one_element %>% 
  filter( source == "C")

control_log <- control %>% 
  mutate(log_pct_change_icp = log(percent_change_icp)) %>% 
  filter(!is.na(log_pct_change_icp))

mean(control_log$log_pct_change_icp)


(mean(control$percent_change_icp))
sqrt(var(control$percent_change_icp))
nrow(control$percent_change_icp)
View(control)


treatment <- one_element %>% 
  filter(source == "P" | source == "NP")

treatment <- one_element %>% 
  filter(source == "N" | source == "NP")

treatment <- one_element %>% 
  filter(source == "N")



mean(treatment$percent_change_icp)
sqrt(var(treatment$percent_change_icp))
nrow(treatment$percent_change_icp)

View(treatment)
View(treatment)


##################################if you had to log transform the model
control_log <- control %>% 
  mutate(log_pct_change_icp = log(percent_change_icp)) %>% 
  filter(!is.na(log_pct_change_icp))

mean(control_log$log_pct_change_icp)
sqrt(var(control_log$log_pct_change_icp))
nrow(control_log$log_pct_change_icp)

View(control_log)

log_one_element <- treated_source %>% 
  filter(element == "Mn") %>% 
  mutate(log_pct_change_icp = log(percent_change_icp)) %>% 
  filter(!is.na(log_pct_change_icp))

treatment <- log_one_element %>% 
  filter(source == "P" | source == "NP")

treatment <- log_one_element %>% 
  filter(source == "N" | source == "NP")

treatment <- log_one_element %>% 
  filter(source == "N")


mean(treatment$log_pct_change_icp)
sqrt(var(treatment$log_pct_change_icp))
nrow(treatment$log_pct_change_icp)

View(control_log)




################################mass loss model ##################################
##################################################################################

# treated source mass loss model

mass_loss_data <- treated_source %>% 
  filter(pct_mass_remaining != "#N/A") %>% 
  mutate(pct_mass_remaining = as.numeric(pct_mass_remaining)) %>% 
  filter(!is.na(pct_mass_remaining)) %>% 
  filter(element == "K") #filter to one element to remove duplicates

View(mass_loss_data)
class(mass_loss_data$pct_mass_remaining)

model <- lmer(pct_mass_remaining ~ source_n*source_p + (1|site/stand) + (1|days) + (1|pct_water), data = mass_loss_data)

anova(model, type = 3)
shapiro.test(resid(model))


### treated destination mass loss model ####
mass_loss_data <- treated_destination %>% 
  filter(pct_mass_remaining != "#N/A") %>% 
  mutate(pct_mass_remaining = as.numeric(pct_mass_remaining)) %>% 
  filter(!is.na(pct_mass_remaining)) %>% 
  filter(element == "K")

model1 <- lmer(pct_mass_remaining ~ destination_n*destination_p + (1|site/stand) + (1|days) + (1|pct_water), data = mass_loss_data)

anova(model1)
shapiro.test(resid(model))

View(tidy_data)






########################## change in CN model
treated_source_cn <- tidy_data %>% 
  filter(!(destination %in% c("N", "P", "NP", "Ca"))) %>% 
  filter(source != "?") %>% 
  distinct(id, .keep_all = TRUE)


cn_source_model <- lmer(percent_change_cn ~ source_n*source_p + stand_age + (1|site/stand/destination) + (1|days) + (1|pct_water), data = treated_source_cn)
shapiro.test(resid(cn_source_model))                      
anova(cn_source_model, type = 3) #n is significant, p and interaction are not

summary(cn_source_model)

emmeans(cn_source_model, specs = "source_n", levels = "0", df = "satterthwaite")

emmeans(cn_source_model, ~ stand_age)

#####################333 treated destination model

treated_dest_cn <- tidy_data %>% 
  filter(!(source %in% c("N", "P", "NP"))) %>% 
  filter(!(destination %in% "Ca")) %>% 
  distinct(id, .keep_all = TRUE)

cn_dest_model <- lmer(percent_change_cn ~ destination_n*destination_p + stand_age + (1|site/stand/destination) + (1|days) + (1|pct_water), data = treated_dest_cn)

shapiro.test(resid(cn_dest_model)) #significantly different but W is 0.94 meaning it is close to normal
plot(resid(cn_dest_model))
hist(resid(cn_dest_model))
cn_dest_anova <- anova(cn_dest_model, type = 3)
print(cn_dest_anova)

summary(cn_dest_model)

library(lmerTest)

r.squaredGLMM(cn_source_model)
#R2m = proportion of total variance explained by the fixed effects - 5%
#R2c = proportion of total variance explained by the fixed and random effects together - 5%

r.squaredGLMM(cn_dest_model)


######################333 POST HOC CODE!!!!!!!!!!!!1 ########################################
library(emmeans)

#main effects
emmeans(cn_dest_model, list(pairwise ~ destination_n), adjust = "tukey", lmer.df = "satterthwaite") 
emmeans(cn_dest_model, list(pairwise ~ destination_p), adjust = "tukey", lmer.df = "satterthwaite")

emmeans(cn_source_model, list(pairwise ~ source_n), adjust = "tukey", lmer.df = "satterthwaite") 
emmeans(cn_source_model, list(pairwise ~ source_p), adjust = "tukey", lmer.df = "satterthwaite")



######################################using emmeans for post hoc on simple effects
library(emmeans)
one_element <- treated_destination %>% 
  filter(element == "P")

model <- lmer(percent_change_icp ~ destination_n*destination_p + stand_age + (1|site/stand/destination) + (1|days) + (1|pct_water), data = one_element)


print(means <- emmeans(cn_dest_model, ~destination_n*destination_p))

print(contrasts <- contrast(means, method = "pairwise", by = "destination_n", interaction = TRUE, lmer.df ="satterthwaite"))

print(contrasts <- contrast(means, method = "pairwise", by = "destination_p", interaction = TRUE, lmer.df ="satterthwaite"))

################################## ####################### post hoc MAIN EFFECTS 

summary(model)

anova_result <- anova(model, type = 3) #main effects, check interaction 
print(anova_result)


print(means <- emmeans(model, ~destination_n*destination_p))

print(contrasts <- contrast(means, method = "pairwise", by = "destination_n", interaction = FALSE, lmer.df ="satterthwaite"))

print(contrasts <- contrast(means, method = "pairwise", by = "destination_p", interaction = FALSE, lmer.df ="satterthwaite"))





#####################################################treated source model
one_element <- treated_source %>% 
  filter(element == "P")

model <- lmer(percent_change_icp ~ source_n*source_p + stand_age + (1|site/stand/destination) + (1|days) + (1|pct_water), data = one_element)


emmeans(model, list(pairwise ~ source_n), adjust = "tukey", lmer.df = "satterthwaite") 
emmeans(model, list(pairwise ~ source_p), adjust = "tukey", lmer.df = "satterthwaite")
emmeans(model, list(pairwise ~ stand_age), adjust = "tukey", lmer.df = "satterthwaite")

####################################### post hoc test on the lmm SIMPLE EFFECTS 
print(means <- emmeans(model, ~source_n*source_p))

print(contrasts <- contrast(means, method = "pairwise", by = "source_n", lmer.df ="satterthwaite" ))

print(contrasts <- contrast(means, method = "pairwise", by = "source_p", lmer.df = "satterthwaite"))






#################################################### Calcium results###############################################33


calcium_data <- filter(tidy_data, destination == "Ca" | destination == "C")
calcium_data <- filter(calcium_data, source == "C")
calcium_data <- mutate(calcium_data, ca_trmt = case_when(
  destination == "Ca" ~ "1",
  destination == "C" ~ "0"
))

ca_one_element <- calcium_data %>% 
  filter(element == "")

ca_model <- lmer(percent_change_icp ~ ca_trmt + (1|site/stand/destination) + (1|days), data = ca_one_element) #log transformed is slightly better but still significantly different and W less than 0.9

#log transformed
ca_model <- lmer(log(percent_change_icp) ~ ca_trmt + (1|site/stand/destination) + (1|days), data = ca_one_element) #log transformed is slightly better but still significantly different and W less than 0.9


shapiro.test(resid(ca_model)) 

summary(ca_model)

anova_result <- anova(ca_model, type = 3) #main effects, check interaction 
print(anova_result)

View(ca_one_element)


get_dupes(calcium_data$id)

unique(calcium_data$element)

######################################################### mass loss model

mass_loss_data <- tidy_data %>% 
  distinct(id, .keep_all = TRUE) %>% 
  filter(pct_mass_remaining != "#N/A" & pct_mass_remaining != "#VALUE!" & pct_mass_remaining != "NA") %>% 
  mutate(pct_mass_remaining = as.numeric(pct_mass_remaining))

View(mass_loss_data)

mass_loss_treated_source <- mass_loss_data %>% 
  filter(destination == "C" & destination != "Ca")

mass_loss_treated_destination <-  mass_loss_data %>% 
  filter(source == "C" & destination != "Ca")

mass_loss_source_model <- lmer(log(pct_mass_remaining) ~ source_n*source_p + (1|site/stand/destination) + (1|days), data = mass_loss_treated_source) 


shapiro.test(resid(mass_loss_source_model)) #log transformation is better
plot(resid(mass_loss_source_model))
hist(resid(mass_loss_source_model))
mass_loss_source_anova <- anova(mass_loss_source_model, type = 3) #no sig interaction 
print(mass_loss_source_anova)
print(cn_dest_anova)

mass_loss_dest_model <- lmer(log(pct_mass_remaining) ~ destination_n*destination_p + (1|site/stand/destination) + (1|days), data = mass_loss_treated_destination) 


shapiro.test(resid(mass_loss_dest_model)) #log transformation is better
plot(resid(mass_loss_dest_model))
hist(resid(mass_loss_dest_model))
mass_loss_dest_anova <- anova(mass_loss_dest_model, type = 3) #no sig interaction 
print(mass_loss_dest_anova)
print(cn_dest_anova)

##################################################################################################
###################################################################################################
##################################### Ca model######################################################

ca_data <- tidy_data %>% 
  filter(fert_type == "Ca" | fert_type == "Control") %>% 
  mutate(fert_type = as.factor(fert_type))
 

P_ca_data <- ca_data %>% 
  filter(element == "P")


#### P conc ca model
P_conc_ca_model <- lmer(percent_change_icp ~ fert_type + (1|site/stand) + (1|days), data = P_ca_data) 
shapiro.test(resid(P_conc_ca_model)) #log transformation is better
anova(P_conc_ca_model, type = 3) #F test is p  = 0.063 


summary(P_conc_ca_model)
####### Mg conc Ca model

Mg_ca_data <- ca_data %>% 
  filter(element == "P")


Mg_conc_ca_model <- lmer(percent_change_icp ~ fert_type + (1|site/stand) + (1|days), data = Mg_ca_data) 
shapiro.test(resid(P_conc_ca_model)) #log transformation is better
anova(Mg_conc_ca_model, type = 3) #F test is p  = 0.063 

summary(Mg_conc_ca_model)

Mg_ca_data %>% 
group_by(fert_type) %>% 
  summarise(mean = mean(percent_change_icp),
            std_error = sd(percent_change_icp)/ sqrt(n()))
