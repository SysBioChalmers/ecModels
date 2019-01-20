#!/usr/bin/env Rscript
# Loading, merging and scaling the proteomic data used in this study
# Benjamin Sanchez

# Load data:
data_1  <- read.csv(file = 'proteomic_data_1.csv', header = TRUE)
data_2  <- read.csv(file = 'proteomic_data_2.csv', header = TRUE)
NTPdata <- read.csv(file = 'int_std_proteinGroups_theoretical_peptides.txt', sep = '\t', header = TRUE)

# Merge both datasets:
all_data <- merge(data_1, data_2, by = c('ï..Protein.IDs','Mol..weight..kDa.'), all.x = TRUE, all.y = TRUE)
names(all_data) <- gsub('ï..','',names(all_data))

# Filter out unneeded variables:
var_pos  <- c(grep('Protein.IDs',names(all_data)),grep('Mol..weight..kDa.',names(all_data)),grep('Intensity.L.',names(all_data)))
all_data <- all_data[,var_pos]
all_data$Intensity.L.x <- NULL
all_data$Intensity.L.y <- NULL

# Change names:
codes  <- read.csv(file = 'proteomic_codes.txt', header = FALSE, sep = '\t')
for(i in 1:length(codes[,1])) {
  names(all_data) <- gsub(paste0(codes[i,1],'\\>'),codes[i,2],names(all_data))
}

# Compute number of theoretical peptides:
intPos  <- names(NTPdata) == 'Intensity'
iBAQpos <- names(NTPdata) == 'iBAQ'
NTP     <- NTPdata[,intPos]/NTPdata[,iBAQpos]
NTPdata <- data.frame(Protein.IDs = NTPdata$Protein.IDs, theo.peptides = NTP)
NTPdata <- NTPdata[!is.nan(NTPdata[,2]),]

# Scale the MS intensities by NTP and MW:
all_data <- merge(NTPdata, all_data, by = 'Protein.IDs', all.x = FALSE, all.y = FALSE)
intPos   <- grep('Intensity.L.',names(all_data))
all_data[,intPos] <- all_data[,intPos]/all_data$theo.peptides
all_data[,intPos] <- all_data[,intPos]*all_data$Mol..weight..kDa.
all_data[,intPos] <- t(t(all_data[,intPos])/colSums(all_data[,intPos], na.rm = TRUE)) #Final units: g/g detected

# Filter out variable values:
names(all_data) <- gsub('\\<Intensity.L.','',names(all_data))
refs      <- cbind(all_data$REF_1,all_data$REF_2,all_data$REF_3,all_data$REF_4)
new_refs  <- cbind(all_data$REF_2_1,all_data$REF_2_2)
means     <- apply(refs,     1, mean, na.rm = TRUE)
new_means <- apply(new_refs, 1, mean, na.rm = TRUE)
to_keep   <- abs(log2(new_means/means)) <= 1
to_keep[is.na(to_keep)] <- TRUE
all_data <- all_data[to_keep,]
all_data$REF_2_1 <- NULL
all_data$REF_2_2 <- NULL

# Save merged dataset:
write.csv(all_data, file = 'merged_proteomic_data.csv', quote = FALSE)
