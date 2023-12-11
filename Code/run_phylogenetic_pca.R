#!/usr/bin/env Rscript

# This script first calculates the median posterior value and 95% HPD 
# interval per response variable per species. It then generates a table 
# of those values for dormancy-capable species only, and performs a 
# phylogenetic principal components analysis (pPCA). Last, it clusters 
# the resulting pPCA scores to identify groups of species that occupy 
# similar regions of the ecophysiological parameter space.
#
# The script needs to be run from the command line, with the user 
# providing the directory where the model outputs (*.Rda files) are 
# found, e.g.:
#
# Rscript run_phylogenetic_pca.R ../Results/MCMCglmm_fits_with_mass_corrections/
#
# The script produces three output files inside the user-provided 
# directory:
#	1) median_and_95_HPD_per_taxon_and_response.csv,
#	2) ppca_torpid_sp.Rda,
#	3) mclust_torpid_sp.Rda.
#
# Note that this script is meant to be run only on the models fitted 
# with MCMCglmm to the entire dataset.

library(MCMCglmm)
library(mclust)
library(phytools)
library(robustbase)

#############################
# F  U  N  C  T  I  O  N  S #
#############################

# This function extracts the median posterior value and the 95% Highest 
# Posterior Density (HPD) interval per response per taxon.
extract_median_and_95_HPD_per_taxon_and_response <- function(working_dir)
{
	chains <- 1:30

	# Load the first .Rda file.
	cat('Now loading 1.Rda ...\n')
	load(paste(working_dir, '/', chains[1], '.Rda', sep = ''))
	
	# Initialise a matrix to store the posterior samples of each 
	# response per species.
	all_Sols <- matrix(
		NA, nrow = 2000 * 30, ncol = ncol(fit$Sol)
	)
	
	# Keep only 2,000 randomly-selected samples per response per species.
	set.seed(1)
	all_Sols[1:2000,] <- fit$Sol[sample(1:nrow(fit$Sol), 2000),]
	
	# Load all remaining chains and store 2,000 randomly selected 
	# samples per chain.
	for ( i in 2:length(chains) )
	{
		rm(fit)
		gc()
	
		cat('Now loading ', chains[i], '.Rda ...\n', sep ='')
		load(paste(working_dir, '/', chains[i], '.Rda', sep = ''))
	
		set.seed(i)
		all_Sols[(((i - 1) * 2000) + 1):(i *2000),] <- fit$Sol[
			sample(1:nrow(fit$Sol), 2000),
		]
	}
	
	# The values we have collected so far represent taxon-specific
	# deviations from the global intercept of each response variable. 
	# The following two lines add the value of the corresponding 
	# intercept to each response, so that we get a proper value per 
	# taxon.
	for ( j in 1:22 )
	{
		all_Sols[
			, grep(colnames(fit$Sol)[j], colnames(fit$Sol))[2:2626]
		] <- all_Sols[,j] + all_Sols[
			, grep(colnames(fit$Sol)[j], colnames(fit$Sol))[2:2626]
		]
	}
	
	# Initialise matrices to store the median values and the lower 
	# and upper bounds of the 95% HPD interval per response and taxon.
	median_matrix <- matrix(nrow = 2626, ncol = 22)
	hpd_95_l_matrix <- matrix(nrow = 2626, ncol = 22)
	hpd_95_u_matrix <- matrix(nrow = 2626, ncol = 22)
	
	# Calculate the median and the 95% HPD interval.
	for ( j in 1:22 )
	{
		median_matrix[,j] <- colMedians(
			all_Sols[,grep(colnames(fit$Sol)[j], colnames(fit$Sol))]
		)
	
		current_HPDs <- HPDinterval(
			mcmc(
				all_Sols[,grep(colnames(fit$Sol)[j], colnames(fit$Sol))]
			)
		)
		hpd_95_l_matrix[,j] <- current_HPDs[,1]
		hpd_95_u_matrix[,j] <- current_HPDs[,2]
	}
	
	# Set the row names of the median matrix to the taxon names.
	row.names(median_matrix) <- gsub(
		'^.*?Species[.]', '', perl = TRUE,
		grep(colnames(fit$Sol)[j], colnames(fit$Sol), value = TRUE)
	)
	
	# Set the name of the root of the phylogeny.
	row.names(median_matrix)[1] <- 'root'
	
	# Set the column names of the median matrix to the names of the 
	# response variables.
	colnames(median_matrix) <- gsub(
		'trait', 'median_', colnames(fit$Sol)[1:22]
	)
	
	# Similarly, rename the row names and column names of the matrices 
	# for the lower and upper bounds of the 95% HPD interval.
	row.names(hpd_95_l_matrix) <- row.names(median_matrix)
	
	colnames(hpd_95_l_matrix) <- gsub(
		'trait', 'HPD_95_l_', colnames(fit$Sol)[1:22]
	)
	
	row.names(hpd_95_u_matrix) <- row.names(median_matrix)
	
	colnames(hpd_95_u_matrix) <- gsub(
		'trait', 'HPD_95_u_', colnames(fit$Sol)[1:22]
	)
	
	# Combine all three matrices in one.
	posterior_dataset_all_taxa <- cbind(
		median_matrix, hpd_95_l_matrix, hpd_95_u_matrix
	)
	
	# Write the resulting matrix to an output_file.
	write.csv(
		posterior_dataset_all_taxa, 
		file = paste(
			working_dir, '/median_and_95_HPD_per_taxon_and_response.csv', 
			sep = ''
		)
	)
	
	return(posterior_dataset_all_taxa)
}

# This function performs a phylogenetic principal components analysis 
# (pPCA) of all ecophysiological variables except for dormancy across 
# dormancy-capable species only. It uses the matrix obtained from the 
# previous function above.
run_ppca <- function(working_dir, posterior_dataset_all_taxa)
{
	
	# Read the main dataset and convert spaces in species' names to 
	# underscores.
	dataset <- read.csv('../Data/dataset.csv')
	dataset$Species <- gsub(' ', '_', dataset$Species)

	# Read the phylogeny.
	tree <- read.tree('../Data/time_calibrated_phylogeny.nwk')
	
	# Remove dormancy-incapable species from the posterior matrix and 
	# from the phylogeny.
	torpid_sp <- dataset$Species[
		!is.na(dataset$Species) & 
		dataset$Dormancy %in% c('Torpor', 'Hibernation')
	]
	posterior_dataset_tips <- posterior_dataset_all_taxa[
		row.names(posterior_dataset_all_taxa) %in% torpid_sp,
	]
	posterior_dataset_tips <- subset(
		posterior_dataset_tips, select = -c(
			median_Dormancy, HPD_95_l_Dormancy, HPD_95_u_Dormancy
		)
	)

	tree_torpid_sp <- drop.tip(
		tree, tree$tip.label[!tree$tip.label %in% torpid_sp]
	)

	# Perform a pPCA and save the result to an output file.
	ppca_MCMCglmm_estimates <- phyl.pca(
		tree_torpid_sp, posterior_dataset_tips, method = 'lambda', 
		mode = 'cor', opt = 'REML'
	)
	save(
		ppca_MCMCglmm_estimates, file = paste(
			working_dir, '/ppca_torpid_sp.Rda', sep = ''
		)
	)

	# Perform clustering of pPCA scores based on Gaussian mixture models, 
	# setting the number of clusters between 1 and 5.
	mclust_fit <- Mclust(ppca_MCMCglmm_estimates$S, G = 1:5)
	
	# Save the clustering result to an output file.
	save(
		mclust_fit, file = paste(
			working_dir, '/mclust_torpid_sp.Rda', sep = ''
		)
	)
}

############################
# M  A  I  N    C  O  D  E #
############################

# Read the working directory provided by the user as a command line 
# argument.
args <- commandArgs(TRUE)
working_dir <- args[1]

# Extract a table of the median posterior value and the bounds of the 
# 95% HPD interval per response per taxon.
posterior_dataset_all_taxa <- extract_median_and_95_HPD_per_taxon_and_response(working_dir)

# Apply a phylogenetic PCA to the above table across extant, 
# dormancy-capable species only.
run_ppca(working_dir, posterior_dataset_all_taxa)
