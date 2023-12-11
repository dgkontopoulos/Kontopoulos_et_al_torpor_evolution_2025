#!/usr/bin/env Rscript

# This script loads MCMCglmm chains and extracts the posterior samples 
# of the liability threshold between daily torpor and hibernation.
#
# The script needs to be run from the command line, with the user 
# providing the directory where the model outputs (*.Rda files) are 
# found, e.g.:
#
# Rscript extract_threshold_between_daily_torpor_and_hibernation.R ../Results/MCMCglmm_fits/
#
# The resulting estimate is written to an output file called 
# 'hibernation_threshold.txt', inside the user-provided directory.

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

# Load the first .Rda file, get the number of samples in the chain, and 
# store the samples of the hibernation threshold in a vector.
cat('Now loading 1.Rda ...\n')
load(paste(working_dir, '/1.Rda', sep = ''))

samples_per_chain <- nrow(fit$Sol)

torpor_to_hibernation_threshold <- rep(NA, samples_per_chain * length(chains))
torpor_to_hibernation_threshold[1:samples_per_chain] <- fit$CP[,1]

# Load all remaining .Rda files, extract their hibernation threshold 
# samples, and append them to the vector.
for ( i in 2:length(chains) )
{
	rm(fit)
	gc()

	cat('Now loading ', chains[i], '.Rda ...\n', sep ='')
	load(paste(working_dir, '/', chains[i], '.Rda', sep = ''))

	torpor_to_hibernation_threshold[
		(((i - 1) * samples_per_chain) + 1):(i * samples_per_chain)
	] <- fit$CP[,1]
}

# Write the median and the 95% HPD interval of the hibernation threshold
# to an output file.
sink(file = paste(working_dir, '/hibernation_threshold.txt', sep = ''))

cat('median: ', median(torpor_to_hibernation_threshold), '\n', sep = '')
cat('95% HPD interval:\n', sep = '')
cat(HPDinterval(mcmc(torpor_to_hibernation_threshold)))

sink()
