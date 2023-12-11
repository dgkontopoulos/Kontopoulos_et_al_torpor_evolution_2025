#!/usr/bin/env Rscript

# This script loads MCMCglmm chains and calculates the correlations 
# between all possible pairs of response variables, based on the 
# variance-covariance matrix. 
#
# The script needs to be run from the command line, with the user 
# providing the directory where the model outputs (*.Rda files) are 
# found, e.g.:
#
# Rscript extract_correlations.R ../Results/MCMCglmm_fits/
#
# The resulting correlation estimates are stored in a file called
# 'correlations.csv', inside the user-provided directory.

library(janitor)
library(MCMCglmm)

#############################
# F  U  N  C  T  I  O  N  S #
#############################

# This function calculates the median posterior correlation between two 
# response variables and its 95% Highest Posterior Density interval.
calculate_correlation <- function(VCV, trait1, trait2)
{
	
	# Get the covariance and the variances of the two responses.
	current_covar <- VCV[,
		paste('trait', trait1, ':trait', trait2, '.Species', sep = '')
	] + VCV[,
		paste('trait', trait1, ':trait', trait2, '.units', sep = '')
	]
	
	current_var1 <- VCV[,
		paste('trait', trait1, ':trait', trait1, '.Species', sep = '')
	] + VCV[,
		paste('trait', trait1, ':trait', trait1, '.units', sep = '')
	]
	
	current_var2 <- VCV[,
		paste('trait', trait2, ':trait', trait2, '.Species', sep = '')
	] + VCV[,
		paste('trait', trait2, ':trait', trait2, '.units', sep = '')
	]
	
	# Calculate the correlation.
	current_cor <- current_covar/sqrt(current_var1 * current_var2)
	
	# Return the median estimate and the 95% HPD interval.
	return(
		list(
			median_cor = round_half_up(median(current_cor), digits = 3), 
			HPD_95_cor = paste(
				'[', 
				round_half_up(HPDinterval(mcmc(current_cor))[,1], digits = 3), 
				',',
				round_half_up(HPDinterval(mcmc(current_cor))[,2], digits = 3),
				']',
				sep = ''
			)
		)
	)
}

# This function generates all possible pairs of response variables,
# calls the function above to calculate correlations, and writes the 
# results to an output file.
extract_correlations <- function(VCV, output_file)
{
	
	# Check whether some of the responses are mass-corrected.
	if ( length(grep('mass_corrected', colnames(VCV))) > 0 )
	{
		traits <- c(
			'Dormancy', 'Body_mass_g', 'BMR_Watt_mass_corrected', 'Brain_size_g_mass_corrected',
			'Max_longevity_years_mass_corrected', 'Range_size_km2', 'Absolute_latitude',
			'Mean_temp', 'SD_temp', 'Annual_precip', 'CV_precip',
			'NPP', 'Migratory', 'Carnivory', 'Herbivory',
			'Fossoriality', 'Aquatic_affinity', 'Hemisphere', 
			'Cathemeral', 'Crepuscular', 'Diurnal', 'Nocturnal'
		)
	} else
	{
		traits <- c(
			'Dormancy', 'Body_mass_g', 'BMR_Watt', 'Brain_size_g',
			'Max_longevity_years', 'Range_size_km2', 'Absolute_latitude',
			'Mean_temp', 'SD_temp', 'Annual_precip', 'CV_precip',
			'NPP', 'Migratory', 'Carnivory', 'Herbivory',
			'Fossoriality', 'Aquatic_affinity', 'Hemisphere', 
			'Cathemeral', 'Crepuscular', 'Diurnal', 'Nocturnal'
		)
	}

	# Get all possible combinations of response variables and store 
	# the correlations in a data frame.
	n_unique_combinations <- length(traits) * (length(traits) - 1) / 2
	correlations <- data.frame(
		Var_1 = rep(NA, n_unique_combinations),
		Var_2 = rep(NA, n_unique_combinations),
		Median_cor = rep(NA, n_unique_combinations),
		HPD_95 = rep(NA, n_unique_combinations)
	)

	row_counter <- 1
	for ( i in 1:(length(traits) - 1) )
	{
		for ( j in (i+1):length(traits) )
		{
			current_cor <- calculate_correlation(
				VCV, traits[i], traits[j]
			)
	
			correlations[row_counter,] <- c(
				traits[i], traits[j], current_cor$median_cor,
				current_cor$HPD_95_cor
			)
			row_counter <- row_counter + 1
		}
	}

	# Write the correlations data frame to the output file.
	write.csv(correlations, file = output_file, row.names = FALSE)
}

############################
# M  A  I  N    C  O  D  E #
############################

# Read the working directory provided by the user as a command line 
# argument.
args <- commandArgs(TRUE)
working_dir <- args[1]

# Look for model fit files in the working directory.
Rda_files_in_working_dir <- list.files(
	path = working_dir, pattern = '^\\d+\\.Rda$'
)

# If there are 5 files only, this means that these come from an
# analysis only across birds or mammals.
#
# If there are 30 files, these come from a whole-dataset analysis.
if ( length(Rda_files_in_working_dir) == 5 )
{
	chains <- 1:5
} else if ( length(Rda_files_in_working_dir) == 30 )
{
	chains <- 1:30
}

# Load the first .Rda file and store its variance-covariance matrix.
cat('Now loading 1.Rda ...\n')
load(paste(working_dir, '/1.Rda', sep = ''))

current_VCV <- fit$VCV

# Load all remaining .Rda files and append their variance-covariance 
# matrices to that of the first file.
for ( i in 2:length(chains) )
{
	cat('Now loading ', chains[i], '.Rda ...\n', sep ='')

	load(paste(working_dir, '/', chains[i], '.Rda', sep = ''))
	current_VCV <- rbind(current_VCV, fit$VCV)
}

# Calculate all correlations between pairs of response variables and 
# write the results to an output file.
extract_correlations(
	current_VCV, paste(working_dir, '/correlations.csv', sep = '')
)
