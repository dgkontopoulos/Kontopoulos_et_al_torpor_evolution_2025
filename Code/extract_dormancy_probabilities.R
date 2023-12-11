#!/usr/bin/env Rscript

# This script loads MCMCglmm chains and extracts the dormancy 
# probabilities per taxon.
#
# The script needs to be run from the command line, with the user 
# providing the directory where the model outputs (*.Rda files) are 
# found, e.g.:
#
# Rscript extract_dormancy_probabilities.R ../Results/MCMCglmm_fits/
#
# The resulting estimates are written to an output file called 
# 'dormancy_probabilities.csv', inside the user-provided directory.
#
# Note that this script is meant to be run only on the models fitted 
# with MCMCglmm to the entire dataset.

library(MCMCglmm)

#############################
# F  U  N  C  T  I  O  N  S #
#############################

# This function adds the posterior samples of each chain that support 
# no dormancy, daily torpor, or hibernation per taxon to the current 
# sample counts.
#
# To determine the dormancy state that each sample supports, the 
# function checks where the dormancy liability lies relatively to the
# threshold values.
modify_dormancy_counts <- function()
{
	
	# Add the samples for the root node.
	
	# If the dormancy liability is below 0 (the threshold for daily 
	# torpor), this means that this sample supports no dormancy.
	dormancy_estimates[['root']]$NO <<- dormancy_estimates[['root']]$NO + length(
		which(fit$Sol[,'traitDormancy'] < 0)
	)
	
	# If the dormancy liability is >= 0 and below the hibernation 
	# threshold, this means that this sample supports daily torpor.
	dormancy_estimates[['root']]$Torpor <<- dormancy_estimates[['root']]$Torpor + length(
		which(
			fit$Sol[,'traitDormancy'] >= 0 & 
			fit$Sol[,'traitDormancy'] < fit$CP[,'cutpoint.traitDormancy.1']
		)
	)
	
	# If the dormancy liability is >= the hibernation threshold, this 
	# means that this sample supports hibernation.
	dormancy_estimates[['root']]$Hibernation <<- dormancy_estimates[['root']]$Hibernation + length(
		which(
			fit$Sol[,'traitDormancy'] >= fit$CP[,'cutpoint.traitDormancy.1']
		)
	)

	# Do the same thing for all remaining taxa.
	lapply(
		2:length(sp_names),
		function(j)
		{
			dormancy_estimates[[sp_names[j]]]$NO <<- dormancy_estimates[[sp_names[j]]]$NO + length(
				which(
					(
						fit$Sol[,'traitDormancy'] + 
						fit$Sol[,
							paste(
								'traitDormancy.Species.', 
								sp_names[j], sep = ''
							)
						]
					) < 0
				)
			)

			dormancy_estimates[[sp_names[j]]]$Torpor <<- dormancy_estimates[[sp_names[j]]]$Torpor + length(
				which(
					(
						fit$Sol[,'traitDormancy'] + 
						fit$Sol[,
							paste(
								'traitDormancy.Species.', 
								sp_names[j], sep = ''
							)
						]
					) >= 0 &
					(
						fit$Sol[,'traitDormancy'] + 
						fit$Sol[,
							paste(
								'traitDormancy.Species.', 
								sp_names[j], sep = ''
							)
						]
					) < fit$CP[,'cutpoint.traitDormancy.1']
				)
			)

			dormancy_estimates[[sp_names[j]]]$Hibernation <<- dormancy_estimates[[sp_names[j]]]$Hibernation + length(
				which(
					(
						fit$Sol[,'traitDormancy'] + 
						fit$Sol[,
							paste(
								'traitDormancy.Species.', 
								sp_names[j], sep = ''
							)
						]
					) >= fit$CP[,'cutpoint.traitDormancy.1']
				)
			)
		}
	)
}

############################
# M  A  I  N    C  O  D  E #
############################

# Read the working directory provided by the user as a command line 
# argument.
args <- commandArgs(TRUE)
working_dir <- args[1]

chains <- 1:30

# Load the first .Rda file.
cat('Now loading 1.Rda ...\n', sep = '')
load(paste(working_dir, '/', chains[1], '.Rda', sep = ''))

# Get the names of extant species and internal nodes.
sp_names <- gsub(
	'traitDormancy.Species.', '', 
	grep('Dormancy', colnames(fit$Sol), value = TRUE)
)
sp_names[sp_names == 'traitDormancy'] <- 'root'

# Initialise a list to store the number of samples in support of each 
# dormancy state per taxon.
dormancy_estimates <- list()
for ( i in sp_names )
{
	dormancy_estimates[[i]] <- list(NO = 0, Torpor = 0, Hibernation = 0)
}

# Store the sample counts from the first chain.
modify_dormancy_counts()

# Load all remaining chains and process their samples.
for ( j in 2:length(chains) )
{
	cat('Now loading ', chains[j], '.Rda ...\n', sep ='')
	load(paste(working_dir, '/', chains[j], '.Rda', sep = ''))

	modify_dormancy_counts()
}

# Create a data frame to store the dormancy probabilities per taxon.
results <- data.frame(
	Taxon = sp_names,
	Prob_NO = rep(NA, length(sp_names)),
	Prob_Torpor = rep(NA, length(sp_names)),
	Prob_Hibernation = rep(NA, length(sp_names))
)

# Calculate the probabilities.
for ( i in 1:nrow(results) )
{
	results$Prob_NO[i] <- dormancy_estimates[[sp_names[i]]]$NO / (
		dormancy_estimates[[sp_names[i]]]$NO +
		dormancy_estimates[[sp_names[i]]]$Torpor +
		dormancy_estimates[[sp_names[i]]]$Hibernation
	)
	results$Prob_Torpor[i] <- dormancy_estimates[[sp_names[i]]]$Torpor / (
		dormancy_estimates[[sp_names[i]]]$NO +
		dormancy_estimates[[sp_names[i]]]$Torpor +
		dormancy_estimates[[sp_names[i]]]$Hibernation
	)
	results$Prob_Hibernation[i] <- dormancy_estimates[[sp_names[i]]]$Hibernation / (
		dormancy_estimates[[sp_names[i]]]$NO +
		dormancy_estimates[[sp_names[i]]]$Torpor +
		dormancy_estimates[[sp_names[i]]]$Hibernation
	)
}

# Write the results to an output file.
write.csv(
	results, file = paste(
		working_dir, '/dormancy_probabilities.csv', sep = ''
	), row.names = FALSE
)
