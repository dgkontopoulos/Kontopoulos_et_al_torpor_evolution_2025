#!/usr/bin/env Rscript

# This script checks if a model that was fitted with MCMCglmm was 
# executed for a sufficient number of MCMC generations. It does so
# by examining the effective sample size and potential scale reduction 
# factor of each model parameter.
#
# The script needs to be run from the command line, with the user 
# providing the directory where the model outputs (*.Rda files) are 
# found, e.g.:
#
# Rscript check_ESS_PSRF.R ../Results/MCMCglmm_fits/

library(MCMCglmm)

#############################
# F  U  N  C  T  I  O  N  S #
#############################

# This function checks if the effective sample size (ESS) of each model 
# parameter per chain is at least 400. If not, this means that the 
# parameter space was not adequately explored.
check_ESS <- function(model_fits)
{
	
	# Initialise a vector to store the ESS for each model parameter per 
	# chain.
	all_ESS_vals <- c()

	# Store ESS values in the vector.
	for ( i in 1:length(model_fits) )
	{
		all_ESS_vals <- c(
			all_ESS_vals,
			effectiveSize(model_fits[[i]]$Sol[,1:22]),
			effectiveSize(model_fits[[i]]$VCV),
			effectiveSize(model_fits[[i]]$CP)
		)
	}

	# Remove non-positive ESS values (which can occur if a parameter 
	# is constant).
	all_ESS_vals <- all_ESS_vals[all_ESS_vals > 0]

	# Then, check if any parameter has an ESS value below 400.
	low_ESS <- which(all_ESS_vals < 400)

	if ( length(low_ESS) > 0 )
	{
		stop('PROBLEM! The ESS is not big enough!')
	} else
	{
		cat('Everything OK!\n', sep = '')
	}
}

# This function checks if the potential scale reduction factor (PSRF; 
# Gelman & Rubin, Stat. Sci., 1992) value of each model parameter is 
# below 1.1. If not, this means that the chains did not converge on 
# statistically indistinguishable posterior distributions.
check_PSRF <- function(model_fits)
{
	
	# Initialise lists to store the model parameter values per chain.
	Sols <- 'mcmc.list('
	VCVs <- 'mcmc.list('
	CPs  <- 'mcmc.list('

	# Store the values in the lists.
	for ( i in 1:length(model_fits) )
	{
		Sols <- paste(
			Sols, 'model_fits[[', i, ']]$Sol[,1:22], ', sep = ''
		)
		VCVs <- paste(
			VCVs, 'model_fits[[', i, ']]$VCV, ', sep = ''
		)
		CPs <- paste(
			CPs, 'model_fits[[', i, ']]$CP, ', sep = ''
		)
	}

	Sols <- sub(', $', ')', Sols)
	VCVs <- sub(', $', ')', VCVs)
	CPs  <- sub(', $', ')', CPs)

	# Calculate the PSRF for each parameter.
	psrf_vals <- c(
		gelman.diag(
			eval(parse(text = Sols)), multivariate = FALSE
		)$psrf[,1],
		gelman.diag(
			eval(parse(text = VCVs)), multivariate = FALSE
		)$psrf[,1],
		gelman.diag(
			eval(parse(text = CPs)), multivariate = FALSE
		)$psrf[,1]
	)

	# Remove NA values (which can occur if a parameter is constant).
	psrf_vals <- psrf_vals[!is.na(psrf_vals)]

	# Then, check if any parameter has a PSRF value of 1.1 or higher.
	high_psrf <- which(psrf_vals >= 1.1)

	if ( length(high_psrf) > 0 )
	{
			stop("PROBLEM! The chains have not converged!")
	} else
	{
			cat('Everything OK!\n', sep = '')
	}
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

chain_groups <- list()

# If there are 5 files only, this means that these come from an
# analysis only across birds or mammals.
#
# If there are 30 files, these come from a whole-dataset analysis.
if ( length(Rda_files_in_working_dir) == 5 )
{
	chain_groups[[1]] <- 1:5
} else if ( length(Rda_files_in_working_dir) == 30 )
{
	
	# For the whole-dataset analysis, we fitted 6 different models 
	# to account for the uncertainty in the body mass and aquatic 
	# affinity of the last common ancestor of Amniota.
	#
	# Thus, here we split the chains by model, so that we can examine 
	# the ESS and PSRF separately for each model.
	chain_groups[[1]] <- 1:5
	chain_groups[[2]] <- 6:10
	chain_groups[[3]] <- 11:15
	chain_groups[[4]] <- 16:20
	chain_groups[[5]] <- 21:25
	chain_groups[[6]] <- 26:30
}

# For each group of chains...
for ( i in 1:length(chain_groups) )
{
	
	# ... load the corresponding .Rda files...
	current_model_fits <- list()

	for ( j in 1:length(chain_groups[[i]]) )
	{
		cat('Now loading ', chain_groups[[i]][j], '.Rda ...\n', sep ='')

		load(paste(working_dir, '/', chain_groups[[i]][j], '.Rda', sep = ''))
		current_model_fits[[j]] <- fit
	}

	# ... and check the ESS and PSRF.
	cat('Now checking the ESS... ', sep = '')
	check_ESS(current_model_fits)

	cat('Now checking the PSRF... ', sep = '')
	check_PSRF(current_model_fits)
}
