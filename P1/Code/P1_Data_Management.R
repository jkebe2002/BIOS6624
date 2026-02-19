#Convert raw hiv data file into file for analysis


#source file
hivdat <- read.csv("~/Downloads/hiv_6624_final.csv")

#convert appropriate variable values to mssing based on code dictionary
hivdat <- hivdat%>%
  mutate(CESD = na_if(CESD, -1))    %>%
  mutate(LIV34 = na_if(LIV34, 9))   %>%
  mutate(FRP = na_if(FRP, 9))       %>%
  mutate(FP = na_if(FP, 9))         %>%
  mutate(HBP = na_if(HBP, 9))       %>%
  mutate(DYSLIP = na_if(DYSLIP, 9)) %>%
  mutate(DIAB = na_if(DIAB, 9))     %>%
  mutate(BMI = na_if(BMI, 999))     %>%
  mutate(income = na_if(income, 9)) %>%
  mutate(HEROPIATE = na_if(HEROPIATE, -9)) 
hivdat$lVload <- log(hivdat$VLOAD, base=10)

#Create variables for differences in outcomes
#only consider observations with all 4 outcomes
hivdat_0 <- hivdat[hivdat$years==0 & complete.cases(hivdat$AGG_MENT, hivdat$lVload, hivdat$AGG_PHYS, hivdat$LEU3N),]
hivdat_2 <- hivdat[hivdat$years==2 & complete.cases(hivdat$AGG_MENT, hivdat$lVload, hivdat$AGG_PHYS, hivdat$LEU3N),]

hivdat_clean <- hivdat_0[hivdat_0$newid %in% hivdat_2$newid,]
#adherence in year 2:
hivdat_clean$ADH <- hivdat_2[hivdat_0$newid %in% hivdat_2$newid,]$ADH

#differences:

hivdat_clean$diff_LEU3N <- hivdat_2[hivdat_0$newid %in% hivdat_2$newid,]$LEU3N - hivdat_0[hivdat_0$newid %in% hivdat_2$newid,]$LEU3N
hivdat_clean$diff_lvload <- hivdat_2[hivdat_0$newid %in% hivdat_2$newid,]$lVload - hivdat_0[hivdat_0$newid %in% hivdat_2$newid,]$lVload
hivdat_clean$diff_AGG_MENT <- hivdat_2[hivdat_0$newid %in% hivdat_2$newid,]$AGG_MENT - hivdat_0[hivdat_0$newid %in% hivdat_2$newid,]$AGG_MENT
hivdat_clean$diff_AGG_PHYS <- hivdat_2[hivdat_0$newid %in% hivdat_2$newid,]$AGG_PHYS - hivdat_0[hivdat_0$newid %in% hivdat_2$newid,]$AGG_PHYS


#education:  5-7 or 1-4 for college educ or higher vs. no college
hivdat_clean$educ_cat <- if_else(hivdat_clean$EDUCBAS>4, "College Degree or Higher", "No College Degree")

# RACE Race 1= White, non-Hispanic
# 2= White, Hispanic
# 3= Black, non-Hispanic
# 4= Black, Hispanic
# 5= American Indian or Alaskan
# Native
# 6= Asian or Pacific Islander
# 7= Other
# 8= Other Hispanic (created
#                    for 2001-03 new recruits)
# Blank= Missin

#create race cat. variable with two options

hivdat_clean$race_cat <- if_else(hivdat_clean$RACE==1, "White, Non-Hispanic", "Other")

# SMOKE Smoking status 1= Never smoked
# 2= Former smoker
# 3= Current smoker
# Blank= Missing
hivdat_clean$smoke_cat <- if_else(hivdat_clean$SMOKE==1, "Never Smoked", 
                                  if_else(hivdat_clean$SMOKE==2, "Former Smoker", 
                                          if_else(!is.na(hivdat_clean$SMOKE), "Current Smoker", NA)))

#remove BMI outliers
hivdat_clean$BMI <- if_else((hivdat_clean$BMI > 0 & hivdat_clean$BMI < 250), hivdat_clean$BMI, NA)




#476 people left in the final dataset...







#Table 1 code:











