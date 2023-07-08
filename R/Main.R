
prolog <- function(data_path, clinical_data) {
  library(ggplot2)
  data <- read.csv(data_path)
  clinical_data <- read.csv(clinical_data)


  # Convert selected columns in clinical_data to factors
  for (i in 1:ncol(clinical_data)) {
    clinical_data[, i] <- as.factor(clinical_data[, i])
  }

  # Set the significance level
  significance_level <- 0.05  # Adjust as per your requirement

  # Create an empty vector to store column names with multiple significant p-values
  columns_with_multiple_significant_p <- c()

  # Iterate over each column
  for (column_number in 3:ncol(clinical_data)) {
    # Perform the glm analysis and extract p-values
    predictor <- clinical_data[, column_number]
    MainModel <- glm(clinical_data$Group_ID ~ predictor, data = clinical_data, family = binomial)
    coef_summary <- summary(MainModel)$coef
    p_values <- coef_summary[, "Pr(>|z|)"]

    # Filter for significant p-values
    significant_p_values <- p_values[p_values < significance_level]

    # Check if there are more than one significant p-values
    if (length(significant_p_values) > 1) {
      # Save the column name
      column_name <- colnames(clinical_data)[column_number]
      columns_with_multiple_significant_p <- c(columns_with_multiple_significant_p, column_name)
    }
  }

  # Create a new dataframe with only the columns that have multiple significant p-values
  filtered_data <- subset(clinical_data, select = columns_with_multiple_significant_p)

  # Convert columns in filtered_data to factors
  filtered_data <- as.data.frame(lapply(filtered_data, as.factor))

  # Transpose the data
  Tdata <- t(data)

  # Read the common file and set the first row as column names
  # Read the common file
  Common <- as.data.frame(Tdata, stringsAsFactors = FALSE)
  # Remove row names
  row.names(Common) <- NULL

  colnames(Common) <- Common[1, ]
  Common <- Common[-1, ]
  # Add the "Groups" column from clinical_data at the start of Common
  Common <- cbind(clinical_data$Groups, Common)
  # Set the column names
  colnames(Common)[1] <- "Groups"

  #library(tibble)
  #Common <- rownames_to_column(Common, var = "Index")
  #Common <- Common[, -1]

  Common$Groups <- as.factor(Common$Groups)

  # Convert selected columns in clinical_data to factors
  for (i in 1:ncol(Common)) {
    Common[, i] <- as.factor(Common[, i])
  }
  write.csv(Common, "Common_file2.csv", row.names = FALSE)


  Common<-read.csv("Common_file2.csv")
  Common$Groups<-as.factor(Common$Groups)

  # Create an empty matrix to store the models
  models2 <- matrix(NA, nrow = ncol(Common)-1, ncol = 5)
  colnames(models2) <- c("Predictor", "Estimate", "Std.Error", "Pr(>|z|)", "Adjusted p-value")
  finalModels2 <- list()

  for (i in 2:ncol(Common)) {
    predictor2 <- Common[, i]
    model2 <- glm(Groups ~ predictor2, data = Common, family = binomial)
    finalModels2[[i]] <- model2
    coef <- summary(model2)$coef
    models2[i - 1, "Predictor"] <- colnames(Common)[i]  # Store predictor name
    models2[i - 1, "Estimate"] <- coef[2, "Estimate"]  # Store coefficient estimate
    models2[i - 1, "Std.Error"] <- coef[2, "Std. Error"]  # Store standard error
    models2[i - 1, "Pr(>|z|)"] <- coef[2, "Pr(>|z|)"]  # Store p-value
  }

  # Adjust p-values using FDR correction
  p_values <- models2[, "Pr(>|z|)"]
  adjusted_p_values <- p.adjust(p_values, method = "fdr")

  # Add adjusted p-values to the models2 matrix
  models2[, "Adjusted p-value"] <- adjusted_p_values

  # Set the significance level
  significance_level <- 0.05  # Adjust as per your requirement

  # Create an empty vector to store column names with multiple significant p-values
  protien_columns_with_multiple_significant_p <- c()

  # Iterate over each column
  for (i in 3:nrow(models2)) {
    # Check if the adjusted p-value is significant
    adjusted_p_value <- as.numeric(models2[i, "Adjusted p-value"])

    if (adjusted_p_value < significance_level) {
      # Save the column name
      column_name <- models2[i, "Predictor"]
      protien_columns_with_multiple_significant_p <- c(protien_columns_with_multiple_significant_p, column_name)
    }
  }

  # Create a new dataframe with only the columns that have multiple significant p-values
  protien_filtered_data <- Common[, c("Groups", protien_columns_with_multiple_significant_p)]

  # Concatenate the dataframes using cbind
  concatenated_data <- cbind(protien_filtered_data, filtered_data)

  # Create an empty list to store the glm models
  glm_models <- list()

  # Get the index of the "Groups" column
  groups_column_index <- which(colnames(concatenated_data) == "Groups")

  # Iterate over each column in protein_filtered_data
  for (i in 1:ncol(protien_filtered_data)) {
    # Add the ith column of protein_filtered_data to filtered_data
    combined_data <- cbind(filtered_data, protien_filtered_data[, i])

    # Perform glm analysis on the Groups column with the combined data
    glm_model <- glm(concatenated_data[, groups_column_index] ~ ., data = combined_data, family = binomial)

    # Store the glm model in the list
    glm_models[[i]] <- glm_model
    print(summary(glm_model))
  }

  # Return the results
  return(list(filtered_data = filtered_data, protien_filtered_data = protien_filtered_data, concatenated_data = concatenated_data, glm_models = glm_models))
}

