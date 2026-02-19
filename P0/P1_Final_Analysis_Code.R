


#Table 1:

library(table1)



t1 <- table1(~|lVload+LEU3N+AGG_MENT+AGG_PHYS+educ_cat+race_cat+smoke_cat+BMI+age, hivdat_clean, caption="Table of descriptive statistics by collection sample.")

t1<- t1flex(t1,tablefn = "qflextable")
t1 %>% flextable::fontsize(size = 7, part = "all") %>% 
  # Reduce font size
  flextable::autofit()
#t1