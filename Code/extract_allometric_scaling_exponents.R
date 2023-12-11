#!/usr/bin/env Rscript

# This script loads MCMCglmm chains and calculates the median posterior 
# values of the allometric scaling exponents of basal metabolic rate, 
# brain mass, and maximum longevity.
#
# The script needs to be run from the command line, with the user 
# providing the directory where the model outputs (*.Rda files) are 
# found, e.g.:
#
# Rscript extract_allometric_scaling_exponents.R ../Results/MCMCglmm_fits/
#
# The resulting estimates are written to an output file called 
# 'allometric_scaling_exponents.csv', inside the user-provided directory.

library(MCMCglmm)

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

# Initialise vectors to store the allometric scaling exponents.
current_exponents_for_BMR <- c()
current_exponents_for_Brain_size <- c()
current_exponents_for_Max_longevity <- c()

# Load each .Rda file and calculate the allometric scaling exponents.
for ( i in 1:length(chains) )
{
	cat('Now loading ', chains[i], '.Rda ...\n', sep ='')
	load(paste(working_dir, '/', chains[i], '.Rda', sep = ''))

	current_exponents_for_BMR <- c(
		current_exponents_for_BMR,
		(
			fit$VCV[,'traitBMR_Watt:traitBody_mass_g.Species'] + 
			fit$VCV[,'traitBMR_Watt:traitBody_mass_g.units']
		) / (
			fit$VCV[,'traitBody_mass_g:traitBody_mass_g.Species'] +
			fit$VCV[,'traitBody_mass_g:traitBody_mass_g.units']
		)
	)

	current_exponents_for_Brain_size <- c(
		current_exponents_for_Brain_size,
		(
			fit$VCV[,'traitBrain_size_g:traitBody_mass_g.Species'] + 
			fit$VCV[,'traitBrain_size_g:traitBody_mass_g.units']
		) / (
			fit$VCV[,'traitBody_mass_g:traitBody_mass_g.Species'] +
			fit$VCV[,'traitBody_mass_g:traitBody_mass_g.units']
		)
	)

	current_exponents_for_Max_longevity <- c(
		current_exponents_for_Max_longevity,
		(
			fit$VCV[,'traitMax_longevity_years:traitBody_mass_g.Species'] + 
			fit$VCV[,'traitMax_longevity_years:traitBody_mass_g.units']
		) / (
			fit$VCV[,'traitBody_mass_g:traitBody_mass_g.Species'] +
			fit$VCV[,'traitBody_mass_g:traitBody_mass_g.units']
		)
	)
}

# Calculate the median value for each exponent.
final_estimates <- data.frame(
	Allometric_scaling_exponent_for_BMR = median(current_exponents_for_BMR),
	Allometric_scaling_exponent_for_Brain_size = median(current_exponents_for_Brain_size),
	Allometric_scaling_exponent_for_Max_longevity = median(current_exponents_for_Max_longevity)
)

# Write the results to an output file.
write.csv(
	final_estimates, file = paste(
		working_dir, '/allometric_scaling_exponents.csv', sep = ''
	),
	row.names = FALSE
)
