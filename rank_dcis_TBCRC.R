

########################################################################################################

# Load necessary libraries
library(dplyr)
library(survival)
library(survminer)
library(readxl)
library(tidyr)
library(pheatmap)
library(ggplot2)

# Set working directory
setwd("C:/Users/DELL/Downloads")

# Load DICS, expression data 
DCIS <- read.csv("TBCRC_reads_matrix_VST_normalized.csv")
print(head(DCIS))

# Read the RANK expression data (skip the first two header rows)
mmc2 <- read_excel("mmc2.xlsx", skip = 2)
print(head(mmc2))

# Convert gene symbols to uppercase for consistency
mmc2$Geneid <- toupper(mmc2$Geneid)

# Remove rows with missing values
mmc2 <- mmc2[complete.cases(mmc2), ]

# Check for duplicate rows and remove them
if (any(duplicated(DCIS))) {
  print("Duplicate rows found and removed.")
  DCIS <- DCIS %>%
    distinct()
} else {
  print("No duplicate rows found.")
}

# Aggregate duplicate gene symbols by taking the mean expression value
if (any(duplicated(DCIS$X))) {
  print("Duplicate gene symbols found. Aggregating data...")
  DCIS <- DCIS %>%
    group_by(X) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
} else {
  print("No duplicate gene symbols found.")
}

# Remove genes with NA values
DCIS <- DCIS[complete.cases(DCIS), ]

# Extract gene IDs from the RANK expression data
rank_genes <- mmc2$Geneid

# Check for missing genes
missing_genes <- setdiff(rank_genes, DCIS$X)
print("Missing genes:")
print(missing_genes)

# Remove missing genes from the RANK signature
rank_genes <- rank_genes[!rank_genes %in% missing_genes]

# Filter DCIS gene expression data to include only RANK signature genes
filtered_DCIS <- DCIS %>%
  filter(X %in% rank_genes)

if (nrow(filtered_DCIS) == 0) {
  stop("Filtered DCIS data is empty. Check the filtering criteria.")
}

# Extract all column names from the expression data
expression_ids <- colnames(filtered_DCIS)[-1]

# View the first few column names
print(head(expression_ids))

# Extract base IDs by removing the last suffix
patient_ids <- sub("_[^_]+$", "", expression_ids)

# View the first few base IDs
print(head(patient_ids))

# Exclude X column as it is not expression values
expression_data <- filtered_DCIS[, -c(1)]

# Transpose the filtered data (genes as columns, samples as rows)
expression_data_transposed <- as.data.frame(t(expression_data))
colnames(expression_data_transposed) <- filtered_DCIS$X

# Set gene names as column names
expression_data_transposed$Patient_ID <- patient_ids

expression_data_transposed$Signature_Score <- rowMeans(expression_data_transposed[, rank_genes], na.rm = TRUE)

# Load DCIS clinical data
DCISClinical <- read.table("DCISClinical.tsv", header = TRUE, sep = "\t")

# Change first column name
colnames(DCISClinical)[1] <- "Patient_ID"

# Merge clinical data with DCIS expression data
merged_data_DCIS <- merge(expression_data_transposed, DCISClinical, by = "Patient_ID")

if (nrow(merged_data_DCIS) == 0) {
  stop("Merged data is empty. Check the merge criteria.")
}

# Define high and low signature score groups
merged_data_DCIS$Signature_Group <- ifelse(
  merged_data_DCIS$Signature_Score > median(merged_data_DCIS$Signature_Score, na.rm = TRUE),
  "High", "Low"
)

# Save merged data to a CSV file
write.csv(merged_data_DCIS, file = "merged_data_DCIS.csv", row.names = FALSE)

survival_data <- merged_data_DCIS %>%
  select(Patient_ID, Age.at.Diagnosis..years., Progression.or.Recurrence, Days.to.Last.Known.Disease.Status, Tumor.Grade, Signature_Score, Signature_Group)

# Convert Progression.or.Recurrence to binary (1 = recurrence, 0 = no recurrence)
survival_data$Event <- ifelse(survival_data$Progression.or.Recurrence == "Yes - Progression or Recurrence", 1, 0)

# Check for missing values
print("Missing values in survival data:")
print(colSums(is.na(survival_data)))

# Remove rows with missing values in critical columns
survival_data <- survival_data %>%
  filter(!is.na(Age.at.Diagnosis..years.) & !is.na(Progression.or.Recurrence) & !is.na(Days.to.Last.Known.Disease.Status) 
         & !is.na(Event) & !is.na(Signature_Score) & !is.na(Signature_Group) & !is.na(Tumor.Grade))

if (nrow(survival_data) == 0) {
  stop("Survival data is empty. Check the data preparation steps.")
}

# Perform Kaplan-Meier survival analysis for DCIS recurrence
fit <- survfit(Surv(Days.to.Last.Known.Disease.Status, Event) ~ Signature_Group, data = survival_data)

# Generate Kaplan-Meier plot for DCIS recurrence
km_plot <- ggsurvplot(
  fit,
  data = survival_data,
  pval = TRUE,                              # Add p-value to the plot
  pval.size = 4,                            # Increase p-value font size
  risk.table = TRUE,                        # Add risk table below the plot
  risk.table.height = 0.35,                 # Increase risk table height
  risk.table.y.text = FALSE,                # Remove y-axis text in risk table
  risk.table.fontsize = 5,                  # Further increase risk table font size
  conf.int = TRUE,                          # Add confidence intervals
  title = "Recurrence-Free Survival Analysis Based on RANK Overexpression Signature", # Updated title
  xlab = "Time (Days)",                     # Label for the x-axis
  ylab = "Recurrence-Free Probability",     # Label for the y-axis
  legend.title = "Signature Group",         # Title for the legend
  legend.labs = c("Low Signature", "High Signature"), # Labels for the groups
  palette = c("blue", "red"),               # Colors for the groups
  break.time.by = 365,                      # Break x-axis by 1 year (365 days)
  ggtheme = theme_classic(),                # Use a classic theme for a clean look
  font.main = c(16, "bold", "black"),       # Customize title font
  font.x = c(14, "bold", "black"),          # Customize x-axis font
  font.y = c(14, "bold", "black"),          # Customize y-axis font
  font.legend = c(12, "plain", "black"),    # Customize legend font
  font.tickslab = c(12, "plain", "black")   # Customize axis tick labels
)

# Customize the theme to ensure a white background
km_plot$plot <- km_plot$plot + theme(
  panel.background = element_rect(fill = "white"),  # White background for the plot area
  plot.background = element_rect(fill = "white"),   # White background for the entire plot
  axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)), # Increase x-axis title size
  axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)), # Increase y-axis title size
  axis.text.x = element_text(size = 12, face = "bold"), # Increase x-axis text size
  axis.text.y = element_text(size = 12, face = "bold")  # Increase y-axis text size
)

# Customize the risk table theme
km_plot$table <- km_plot$table + theme(
  panel.background = element_rect(fill = "white"),  # White background for the table area
  plot.background = element_rect(fill = "white"),   # White background for the entire table
  axis.title.x = element_text(size = 12, face = "bold"), # Customize x-axis title
  axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), # Rotate x-axis text
  axis.text.y = element_blank()                    # Remove y-axis text in risk table
)

# Combine the Kaplan-Meier plot and risk table into a single layout
combined_plot <- arrange_ggsurvplots(
  list(km_plot),                            # List of ggsurvplot objects
  ncol = 1,                                 # Number of columns in the layout
  nrow = 2,                                 # Number of rows in the layout
  print = FALSE                             # Do not print the plot immediately
)

# Save the combined plot (Kaplan-Meier plot + risk table) in high resolution
ggsave(
  filename = "Combined_KM_Plot_DCIS_RANK_Signature.png", # File name for the combined plot
  plot = combined_plot,                      # Save the combined plot
  width = 20,                                # Increase width to provide more horizontal space
  height = 20,                               # Increase height to provide more vertical space
  dpi = 300,                                 # Resolution (dots per inch)
  limitsize = FALSE                          # Allow larger dimensions
)

# Display the combined plot (optional)
print(combined_plot)

# Perform Cox proportional hazards analysis 
cox_model <- coxph(
  Surv(Days.to.Last.Known.Disease.Status, Event) ~ Signature_Score + Age.at.Diagnosis..years., 
  data = survival_data
)

# Summarize the Cox model
cox_summary <- summary(cox_model)
print(cox_summary)

# Save Cox model results to a CSV file
cox_DCIS_TBCRC_results <- data.frame(
  Variable = names(cox_model$coefficients), # Extract variable names
  HR = exp(cox_model$coefficients),         # Hazard ratios (exponentiated coefficients)
  CI_lower = exp(confint(cox_model)[, 1]),  # Lower bound of 95% confidence interval
  CI_upper = exp(confint(cox_model)[, 2]),  # Upper bound of 95% confidence interval
  P_value = cox_summary$coefficients[, "Pr(>|z|)"] # P-values
)

# Write the results to a CSV file
write.csv(cox_DCIS_TBCRC_results, file = "cox_DCIS_TBCRC_results_signature_score.csv", row.names = FALSE)
