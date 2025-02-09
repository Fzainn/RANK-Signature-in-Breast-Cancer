
# Load required libraries
library(dplyr)
library(survival)
library(survminer)
library(readxl)
library(tidyr)
library(pheatmap)



# Set working directory
setwd("C:/Users/DELL/Downloads")

# Read the RANK expression data (skip the first two header rows)
mmc2 <- read_excel("mmc2.xlsx", skip = 2)

# Convert gene symbols to uppercase for consistency
mmc2$Geneid <- toupper(mmc2$Geneid)


mmc2 <- mmc2[complete.cases(mmc2), ]

# Load the METABRIC gene expression data
METBgene_expression <- read.table("data/data_mrna_illumina_microarray.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Check for duplicate rows and remove them
if (any(duplicated(METBgene_expression))) {
  print("Duplicate rows found and removed.")
  METBgene_expression <- METBgene_expression %>%
    distinct()
} else {
  print("No duplicate rows found.")
}

# Aggregate duplicate gene symbols by taking the mean expression value
if (any(duplicated(METBgene_expression$Hugo_Symbol))) {
  print("Duplicate gene symbols found. Aggregating data...")
  METBgene_expression <- METBgene_expression %>%
    group_by(Hugo_Symbol) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
} else {
  print("No duplicate gene symbols found.")
}

# Remove genes with NA values
METBgene_expression <- METBgene_expression[complete.cases(METBgene_expression), ]

# Extract gene IDs from the RANK expression data
rank_genes <- mmc2$Geneid

# Check for missing genes
missing_genes <- setdiff(rank_genes, METBgene_expression$Hugo_Symbol)
print("Missing genes:")
print(missing_genes)

# Remove missing genes from the RANK signature
rank_genes <- rank_genes[!rank_genes %in% missing_genes]

# Filter METABRIC gene expression data to include only RANK signature genes
filtered_metabric <- METBgene_expression %>%
  filter(Hugo_Symbol %in% rank_genes)

# Transpose the filtered data (genes as columns, samples as rows)
expression_data <- filtered_metabric[, -c(1, 2)]  # Exclude Hugo_Symbol and Entrez_Gene_Id
filtered_metabric_transposed <- as.data.frame(t(expression_data))
colnames(filtered_metabric_transposed) <- filtered_metabric$Hugo_Symbol
filtered_metabric_transposed$PATIENT_ID <- rownames(filtered_metabric_transposed)

filtered_metabric_transposed$Signature_Score <- rowMeans(filtered_metabric_transposed[, rank_genes], na.rm = TRUE)




# Load the METABRIC clinical data
METBclinical_data <- read.table("data/data_clinical_patient.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Replace hyphens with dots in PATIENT_ID
METBclinical_data <- METBclinical_data %>%
  mutate(PATIENT_ID = gsub("-", ".", PATIENT_ID))

# Merge gene expression data with clinical data
merged_data <- merge(filtered_metabric_transposed, METBclinical_data, by = "PATIENT_ID")

# Define high and low signature score groups
merged_data$Signature_Group <- ifelse(
  merged_data$Signature_Score > median(merged_data$Signature_Score, na.rm = TRUE),
  "High", "Low"
)



# Save merged data to a CSV file
write.csv(merged_data, file = "merged_data.csv", row.names = FALSE)





# Extract relevant columns for survival analysis
survival_data <- merged_data %>%
  select(PATIENT_ID, OS_MONTHS, OS_STATUS, AGE_AT_DIAGNOSIS, ER_IHC, Signature_Score, Signature_Group) %>%
  mutate(Event_Status = ifelse(OS_STATUS == "1:DECEASED", 1, 0))  # Convert OS_STATUS to binary




survival_data$OS_MONTHS <- as.numeric(survival_data$OS_MONTHS)
survival_data$Event_Status <- as.numeric(survival_data$Event_Status)
survival_data$AGE_AT_DIAGNOSIS <- as.numeric(survival_data$AGE_AT_DIAGNOSIS)
survival_data$ER_IHC <- as.factor(survival_data$ER_IHC)


# Check for missing values
print("Missing values in survival data:")
print(colSums(is.na(survival_data)))

# Remove rows with missing values in critical columns
survival_data <- survival_data %>%
  filter(!is.na(OS_MONTHS) & !is.na(Event_Status) & !is.na(AGE_AT_DIAGNOSIS) & !is.na(ER_IHC))

# Correct the typo in the ER_IHC column
levels(survival_data$ER_IHC)[levels(survival_data$ER_IHC) == "Positve"] <- "Positive"

# Re-convert to factor with correct levels
survival_data$ER_IHC <- factor(survival_data$ER_IHC, levels = c("Negative", "Positive"))

table(survival_data$ER_IHC)


# Perform Kaplan-Meier survival analysis for the signature score
fit <- survfit(Surv(OS_MONTHS, Event_Status) ~ Signature_Group, data = survival_data)

# Generate Kaplan-Meier plot
km_plot <- ggsurvplot(
  fit,
  data = survival_data,
  pval = TRUE,
  risk.table = TRUE,
  title = "Survival Analysis Based on RANK Overexpression Signature"
)

# Save the Kaplan-Meier plot
ggsave(filename = "KM_plot_RANK_signature.png", plot = km_plot$plot, width = 10, height = 8)

# Perform Cox proportional hazards analysis for the signature score
cox_model <- coxph(Surv(OS_MONTHS, Event_Status) ~ Signature_Score + AGE_AT_DIAGNOSIS + ER_IHC, data = survival_data)
cox_summary <- summary(cox_model)
print(cox_summary)

# Save Cox model results to a CSV file
cox_results <- data.frame(
  Variable = names(cox_model$coefficients),
  HR = exp(cox_model$coefficients),
  CI_lower = exp(confint(cox_model)[, 1]),
  CI_upper = exp(confint(cox_model)[, 2]),
  P_value = cox_summary$coefficients[, "Pr(>|z|)"]
)
write.csv(cox_results, file = "cox_results_signature_score.csv", row.names = FALSE)





















