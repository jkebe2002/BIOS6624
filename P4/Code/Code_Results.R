###################################################################################
# PLEASE READ:                                                                    #
# sim_results is the product of running the simulation in Code_Functions          #
# To produce this dataframe quickly you can just reduce foreach iteration number. #
# Thank you!                                                                      #
###################################################################################

sim_results <- read.csv("~/Documents/GitHub/BIOS6624/P4/sim_results.csv")
library(dplyr)
library(ggplot2)



n_sim <- 10000 #Number of iterations in foreach loop in Code_Functions

sim_results_compiled <- sim_results |>
  group_by(
    case,method,var
  ) |>
  summarize(
    mean_bias <- mean(bias),
    coverage_prop <- sum(coverage) / n_sim,
    
  )


# summarise(
#   n_cases  = n(),
#   mean_age = mean(age_years, na.rm=T),
#   max_age  = max(age_years, na.rm=T),
#   min_age  = min(age_years, na.rm=T),
#   n_males  = sum(gender == "m", na.rm=T))

# Table of 

library(dplyr)
library(tidyr)
library(kableExtra)

#######################################################################################################
# The following code produces the massive table of bias and coverage per variable per method per case #
#######################################################################################################
# -- 1. Prep ---------------------------------------------------------------
sim_prepped <- sim_results |>
  mutate(
    N   = if_else(startsWith(case, "1"), "N=250", "N=500"),
    rho = case_when(
      endsWith(case, "a") ~ "rho=0",
      endsWith(case, "b") ~ "rho=0.35",
      endsWith(case, "c") ~ "rho=0.7"
    ),
    method = factor(method, levels = c("P-Value","AIC","BIC","LASSO1","LASSO2","ENET1","ENET2"))
  ) |>
  summarise(
    Bias     = mean(bias,     na.rm = TRUE),
    Coverage = mean(coverage, na.rm = TRUE),
    .by = c(var, method, N, rho)
  )

# -- 2. Build one table for a given rho ------------------------------------
make_table <- function(data, rho_val) {
  
  wide <- data |>
    filter(rho == rho_val) |>
    pivot_wider(
      names_from  = var,
      values_from = c(Bias, Coverage),
      names_glue  = "{var}_{.value}"
    )
  
  var_labels <- paste0("V0", 1:5)
  var_cols   <- c(rbind(
    paste0(var_labels, "_Bias"),
    paste0(var_labels, "_Coverage")
  ))
  
  header        <- c(2, rep(2, 5))
  names(header) <- c(" ", var_labels)
  
  wide |>
    select(N, method, all_of(var_cols)) |>
    arrange(N, method) |>
    kbl(
      digits    = 3,
      caption   = paste0("Bias and Coverage by Method and Variable (", rho_val, ")"),
      col.names = c("N", "Method", rep(c("Bias", "Coverage"), 5))
    ) |>
    kable_classic(full_width = FALSE) |>
    add_header_above(header) |>
    pack_rows("N = 250", 1, 7) |>
    pack_rows("N = 500", 8, 14)
}



make_combined_table <- function(data) {

  make_wide <- function(rho_val) {
    var_labels <- paste0("V0", 1:5)
    var_cols   <- c(rbind(
      paste0(var_labels, "_Bias"),
      paste0(var_labels, "_Coverage")
    ))

    data |>
      filter(rho == rho_val) |>
      pivot_wider(
        names_from  = var,
        values_from = c(Bias, Coverage),
        names_glue  = "{var}_{.value}"
      ) |>
      arrange(N, method) |>
      select(N, method, all_of(var_cols))
  }

  combined <- bind_rows(
    make_wide("rho=0"),
    make_wide("rho=0.35"),
    make_wide("rho=0.7")
  )

  header        <- c(2, rep(2, 5))
  names(header) <- c(" ", paste0("V0", 1:5))

  combined |>
    kbl(
      digits    = 3,
      caption   = "Bias and Coverage by Method, Variable, and Correlation Structure",
      col.names = c("N", "Method", rep(c("Bias", "Cov."), 5)),  # abbreviated to save space
      booktabs  = TRUE,    # cleaner lines for PDF
      linesep   = ""       # suppress extra space every 5 rows
    ) |>
    kable_styling(
      latex_options = c("hold_position", "scale_down"),  # scale_down fits wide tables to page
      font_size     = 8                                   # reduce font size if needed
    ) |>
    add_header_above(header) |>
    pack_rows("rho = 0",    1,  14) |>
    pack_rows("rho = 0.35", 15, 28) |>
    pack_rows("rho = 0.7",  29, 42) |>
    pack_rows("N = 250",  1,  7,  indent = TRUE) |>
    pack_rows("N = 500",  8,  14, indent = TRUE) |>
    pack_rows("N = 250",  15, 21, indent = TRUE) |>
    pack_rows("N = 500",  22, 28, indent = TRUE) |>
    pack_rows("N = 250",  29, 35, indent = TRUE) |>
    pack_rows("N = 500",  36, 42, indent = TRUE)
}


########################
#Plot of Type 2 errors #
########################

library(dplyr)
library(tidyr)

# 1. Summaries
type2_by_method <- sim_results |>
  group_by(method, var) |>
  summarize(type2_total = sum(type2), .groups = "drop") |>
  mutate(grouping = sprintf("Method = %s", method))

type2_by_case <- sim_results |>
  group_by(case, var) |>
  summarize(type2_total = sum(type2), .groups = "drop") |>
  mutate(grouping = sprintf("Case = %s", case))

# 2. Combine
type2_combined <- bind_rows(type2_by_method, type2_by_case)

# 3. Pivot wider so columns are V01, V02, ..., V05
type2_wide <- type2_combined |>
  select(grouping, var, type2_total) |>
  pivot_wider(
    names_from = var,
    values_from = type2_total
  ) |>
  arrange(grouping)

type2plot<- ggplot(type2_combined, aes(x = var, y = type2_total, color = grouping, group = grouping)) +
  geom_line(size = .8) +
  geom_point(size = 1.2) +
  labs(
    x = "Variable",
    y = "Type II Errors",
    color = "Condition",
    title = "Type II Errors Across Variables"
  ) +
  theme_minimal(base_size = 9)+
  theme(text = element_text(family = "serif"))











############
#Heat Maps##
############

library(dplyr)
library(tidyr)
library(ggplot2)

# -- 1. Summarise: per simulation, did ALL variables meet the criterion? ----
# Assumes one row per simulation x variable x method x case

# change sim_data
# change sim_id (should just be floor(row #) / 42)
sim_results$sim_id <- floor((seq_len(nrow(sim_results)) - 1) / 210)
heatmap_data <- sim_results |>
  mutate(
    N   = if_else(startsWith(case, "1"), "N=250", "N=500"),
    rho = case_when(
      endsWith(case, "a") ~ "rho=0",
      endsWith(case, "b") ~ "rho=0.35",
      endsWith(case, "c") ~ "rho=0.7"
    ),
    method = factor(method, levels = c("P-Value","AIC","BIC","LASSO1","LASSO2","ENET1","ENET2"))
  ) |>
  # For each simulation run, check if ALL variables satisfy the criterion
  summarise(
   all_trueP  = as.integer(all(trueP  == 1, na.rm = TRUE)),
    all_falseP = as.integer(all(falseP == 1, na.rm = TRUE)),
   any_type1 = as.integer(any(type1 > 0, na.rm = TRUE)),
   any_type2 = as.integer(any(type2 > 0, na.rm = TRUE)),
   all_type1 = mean(type1),
   all_type2 = sum(type2),
    trueP=sum(trueP),
    .by = c(sim_id, method, case, N, rho)
  ) |>
  # Then average across simulations to get rates
  summarise(
    rate_trueP  = mean(all_trueP),
    rate_falseP = mean(all_falseP),
    rate_type1_v = mean(all_type1),
    rate_type2_v = mean(all_type2),
    rate_type1_m = mean(any_type1),
    rate_type2_m = mean(any_type2),
    
    .by = c(method, case, N, rho)
  ) |>
  # Nice case label for x-axis
  mutate(case_label = paste0(N, "\n", rho))

# -- 2. Helper to make one heatmap -----------------------------------------
make_heatmap <- function(data, rate_var, title, low, mid, high) {
  data |>
    ggplot(aes(x = case_label, y = method, fill = {{ rate_var }})) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round({{rate_var}}, 3)),
              size = 2.3,        
              family = "serif") +
    scale_fill_gradient2(
      low      = low,
      mid      = mid,
      high     = high,
      midpoint = 0.5,
      #limits   = c(0, 1),
      #labels   = scales::percent,
      name     = "Rate"
    ) +
    scale_x_discrete(limits = unique(data$case_label)) +
    labs(title = title, x = NULL, y = "Method") +
    theme_minimal(base_size = 8) +   # also reduced base_size
    theme(
      panel.grid  = element_blank(),
      text        = element_text(family = "serif"),
      axis.text.x = element_text(hjust = 0.5)
    )
}

# -- 3. Plot both heatmaps -------------------------------------------------
p_trueP <- make_heatmap(
  heatmap_data, rate_trueP,
  low      = "#d73027",
  mid      = "#ffffbf",
  high     = "#1a9850",
  "True Positive Rate: \n All significant variables included in Model"
)

p_falseP <- make_heatmap(
  heatmap_data, rate_falseP,
  low      = "#1a9850",
  mid      = "#ffffbf",
  high     = "#d73027",
  "False Positive Rate: \n Any non-significant variables included"
)

p_type1 <- make_heatmap(
  heatmap_data, rate_type1_v,
  low      = "#1a9850",
  mid      = "#ffffbf",
  high     = "#d73027",
  "Type 1 Error (# Errors per Model)"
)

p_type2 <- make_heatmap(
  heatmap_data, rate_type2_v,
  low      = "#1a9850",
  mid      = "#ffffbf",
  high     = "#d73027",
  "Type 2 Error (# Errors per Model)"
)

p_type1_m <- make_heatmap(
  heatmap_data, rate_type1_m,
  low      = "#1a9850",
  mid      = "#ffffbf",
  high     = "#d73027",
  "Type 1 Error (Model Level)"
)

p_type2_m <- make_heatmap(
  heatmap_data, rate_type2_m,
  low      = "#1a9850",
  mid      = "#ffffbf",
  high     = "#d73027",
  "Type 2 Error (Model Level)"
)

# Print individually or combine with patchwork
library(patchwork)
p_trueP / p_falseP
p_type1 / p_type2
p_type1_m / p_type2_m





library(gt)

ve_coverage_table <- sim_results |>
  mutate(
    case_label = case_when(
      case == "1a" ~ "N=250, rho=0",
      case == "1b" ~ "N=250, rho=0.35",
      case == "1c" ~ "N=250, rho=0.7",
      case == "2a" ~ "N=500, rho=0",
      case == "2b" ~ "N=500, rho=0.35",
      case == "2c" ~ "N=500, rho=0.7"
    )
  ) |>
  distinct(sim_id, method, case_label, .keep_all = TRUE) |>
  summarise(
    VE       = -round(mean(VE, na.rm = TRUE), 3),
    coverage = round(1 - (sum(type1, na.rm = TRUE) / (n() * 15)), 3),
    .by = c(method, case_label)
  ) |>
  mutate(cell = paste0(VE, " (", coverage, ")")) |>
  pivot_wider(
    names_from  = case_label,
    values_from = cell,
    id_cols     = method
  ) |>
  select(method,
         `N=250, rho=0`, `N=250, rho=0.35`, `N=250, rho=0.7`,
         `N=500, rho=0`, `N=500, rho=0.35`, `N=500, rho=0.7`)
library(gt)

###################################################
# Bias and coverage table for variables 6 thru 20 #
###################################################
ve_coverage_table <- ve_coverage_table |>
  gt() |>
  tab_spanner(
    label   = "N = 250",
    columns = c(`N=250, rho=0`, `N=250, rho=0.35`, `N=250, rho=0.7`)
  ) |>
  tab_spanner(
    label   = "N = 500",
    columns = c(`N=500, rho=0`, `N=500, rho=0.35`, `N=500, rho=0.7`)
  ) |>
  cols_label(
    method          = "Method",
    `N=250, rho=0`  = md("rho = 0"),
    `N=250, rho=0.35` = md("rho = 0.35"),
    `N=250, rho=0.7`  = md("rho = 0.70"),
    `N=500, rho=0`  = md("rho = 0"),
    `N=500, rho=0.35` = md("rho = 0.35"),
    `N=500, rho=0.7`  = md("rho = 0.70")
  ) |>
  tab_source_note(
    source_note = "Values reported as Bias (Coverage). Bias is the average coefficient estimate for null variables retained in the model (true value = 0). Coverage is calculated as 1 - (Type I error rate)."
  ) |>
  tab_options(
    table.font.names        = "Times New Roman",
    table.font.size         = 12,
    heading.align           = "left",
    column_labels.font.weight = "bold",
    table.border.top.style  = "solid",
    table.border.top.width  = px(2),
    table.border.bottom.style = "solid",
    table.border.bottom.width = px(2),
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.width = px(1)
  ) |> tab_caption(caption = "Average bias and coverage for null variables")

