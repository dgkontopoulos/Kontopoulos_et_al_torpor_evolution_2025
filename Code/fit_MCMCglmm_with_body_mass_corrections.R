#!/usr/bin/env Rscript

# This script fits a model with MCMCglmm to estimate the correlation 
# structure among torpor and 21 ecophysiological variables across 
# 1,338 endotherms, accounting for phylogeny. 3 of these variables 
# (BMR, brain mass, and maximum longevity) are corrected for body mass.
#
# 30 different chains have been specified in this script, allowing the 
# model to 
#	a) thoroughly explore the parameter space and 
#	b) take into account the uncertainty in the aquatic affinity and 
#	   body mass of the last common ancestor of Amniota.
#
# The script needs to be run from the command line, with the user 
# providing the chain ID as follows:
#
# Rscript fit_MCMCglmm_with_body_mass_corrections.R 1
# Rscript fit_MCMCglmm_with_body_mass_corrections.R 2
# ...
# Rscript fit_MCMCglmm_with_body_mass_corrections.R 29
# Rscript fit_MCMCglmm_with_body_mass_corrections.R 30
#
# This script will run for a few days depending on the CPU (~7.5 days on 
# an AMD EPYC 7702 CPU core) and will need ~67 GB of memory. In the end, 
# it will produce an .Rda file per chain, in the 
# ../Results/MCMCglmm_fits_with_mass_corrections/ directory.
#
# To ensure that the chains have converged on statistically equivalent 
# posterior distributions and that the parameter space has been adequately 
# explored, the user needs to run the check_ESS_PSRF.R script next.

library(ape)
library(MCMCglmm)

#############################
# F  U  N  C  T  I  O  N  S #
#############################

# This function reads the dataset and the phylogeny, and prepares all 
# the variables needed for model fitting.
prepare_vars_for_MCMCglmm <- function(chain_id)
{
	
	# Read the dataset and convert spaces in species' names to underscores.
	dataset <- read.csv('../Data/dataset.csv')
	dataset$Species <- gsub(' ', '_', dataset$Species)
	
	# Convert categorical variables to binary factors or ordered factors.

	dataset$Torpor <- factor(
		dataset$Torpor, levels = c('NO', 'Torpor', 'Hibernation')
	)
	
	dataset$Migratory <- factor(
		dataset$Migratory, levels = c('NO', 'YES')
	)

	# Split diet information into two binary variables: carnivory and herbivory.
	dataset$Carnivory <- rep(NA, nrow(dataset))
	dataset$Carnivory[dataset$Diet == 'herbivore'] <- 'NO'
	dataset$Carnivory[dataset$Diet %in% c('carnivore', 'omnivore')] <- 'YES'
	dataset$Carnivory <- factor(
		dataset$Carnivory, levels = c('NO', 'YES')
	)
	
	dataset$Herbivory <- rep(NA, nrow(dataset))
	dataset$Herbivory[dataset$Diet == 'carnivore'] <- 'NO'
	dataset$Herbivory[dataset$Diet %in% c('herbivore', 'omnivore')] <- 'YES'
	dataset$Herbivory <- factor(
		dataset$Herbivory, levels = c('NO', 'YES')
	)
	
	dataset$Fossoriality <- factor(
		dataset$Fossoriality, levels = c('nonfossorial', 'semifossorial', 'fossorial')
	)
	
	dataset$Aquatic_affinity <- factor(
		dataset$Aquatic_affinity, levels = c('very_low', 'low', 'moderate', 'high')
	)
	
	# Convert latitude to absolute latitude.
	dataset$Absolute_latitude <- abs(dataset$Mid_range_lat_dd)
	
	# Create a binary hemisphere variable.
	dataset$Hemisphere <- rep(NA, nrow(dataset))
	dataset$Hemisphere[which(dataset$Mid_range_lat_dd > 0)] <- 'northern'
	dataset$Hemisphere[which(dataset$Mid_range_lat_dd < 0)] <- 'southern'
	dataset$Hemisphere <- factor(
		dataset$Hemisphere, levels = c('southern', 'northern')
	)
	
	# For activity patterns, create 4 binary variables: cathemeral, crepuscular, diurnal, nocturnal.
	dataset$Cathemeral <- rep('NO', nrow(dataset))
	dataset$Cathemeral[which(dataset$Daily_activity == 'cathemeral')] <- 'YES' 
	dataset$Cathemeral[which(is.na(dataset$Daily_activity))] <- NA
	dataset$Cathemeral <- factor(
		dataset$Cathemeral, levels = c('NO', 'YES')
	)
	
	dataset$Crepuscular <- rep('NO', nrow(dataset))
	dataset$Crepuscular[which(dataset$Daily_activity == 'crepuscular')] <- 'YES' 
	dataset$Crepuscular[which(is.na(dataset$Daily_activity))] <- NA
	dataset$Crepuscular <- factor(
		dataset$Crepuscular, levels = c('NO', 'YES')
	)
	
	dataset$Diurnal <- rep('NO', nrow(dataset))
	dataset$Diurnal[which(dataset$Daily_activity == 'diurnal')] <- 'YES' 
	dataset$Diurnal[which(is.na(dataset$Daily_activity))] <- NA
	dataset$Diurnal <- factor(
		dataset$Diurnal, levels = c('NO', 'YES')
	)
	
	dataset$Nocturnal <- rep('NO', nrow(dataset))
	dataset$Nocturnal[which(dataset$Daily_activity == 'nocturnal')] <- 'YES' 
	dataset$Nocturnal[which(is.na(dataset$Daily_activity))] <- NA
	dataset$Nocturnal <- factor(
		dataset$Nocturnal, levels = c('NO', 'YES')
	)
	
	# Read the phylogeny and add names to nodes.
	tree <- read.tree('../Data/time_calibrated_phylogeny.nwk')
	tree$node.label <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)
	
	tree$node.label[tree$node.label == getMRCA(tree, c('Tachyglossus_aculeatus', 'Arctocebus_calabarensis'))] <- 'Mammalia'
	tree$node.label[tree$node.label == getMRCA(tree, c('Isoodon_obesulus', 'Arctocebus_calabarensis'))] <- 'Theria'
	tree$node.label[tree$node.label == getMRCA(tree, c('Cabassous_centralis', 'Arctocebus_calabarensis'))] <- 'Placentalia'
	tree$node.label[tree$node.label == getMRCA(tree, c('Struthio_camelus', 'Sitta_carolinensis'))] <- 'Aves'
	tree$node.label[tree$node.label == getMRCA(tree, c('Apteryx_australis', 'Struthio_camelus'))] <- 'Palaeognathae'
	tree$node.label[tree$node.label == getMRCA(tree, c('Cygnus_olor', 'Sitta_carolinensis'))] <- 'Neognathae'
	tree$node.label[tree$node.label == getMRCA(tree, c('Dendragapus_obscurus', 'Cygnus_olor'))] <- 'Galloanserae'
	tree$node.label[tree$node.label == getMRCA(tree, c('Caprimulgus_guttatus', 'Sitta_carolinensis'))] <- 'Neoaves'
	
	# Add palaeobiological information for specific deep tree nodes.
	
	# Amniota (root of the tree)
	dataset[nrow(dataset) + 1,] <- NA
	dataset$Carnivory[nrow(dataset)] <- 'YES'
	dataset$Herbivory[nrow(dataset)] <- 'NO'
	Amniota_row <- nrow(dataset)
	
	# Mammalia
	dataset[nrow(dataset) + 1,] <- NA
	dataset$Species[nrow(dataset)] <- 'Mammalia'
	dataset$Carnivory[nrow(dataset)] <- 'YES'
	dataset$Herbivory[nrow(dataset)] <- 'NO'
	
	# Theria
	dataset[nrow(dataset) + 1,] <- NA
	dataset$Species[nrow(dataset)] <- 'Theria'
	dataset$Carnivory[nrow(dataset)] <- 'YES'
	dataset$Herbivory[nrow(dataset)] <- 'NO'
	
	# Placentalia
	dataset[nrow(dataset) + 1,] <- NA
	dataset$Species[nrow(dataset)] <- 'Placentalia'
	dataset$Body_mass_g[nrow(dataset)] <- 170
	dataset$Carnivory[nrow(dataset)] <- 'YES'
	dataset$Herbivory[nrow(dataset)] <- 'NO'
	dataset$Brain_size_g[nrow(dataset)] <- 1.036 * 1.17
	dataset$BMR_Watt[nrow(dataset)] <- 7.68 * 170 * 20.1 / 3600 / 3.53
	
	# Aves
	dataset[nrow(dataset) + 1,] <- NA
	dataset$Species[nrow(dataset)] <- 'Aves'
	dataset$Brain_size_g[nrow(dataset)] <- 1.036 * 7
	
	# Palaeognathae
	dataset[nrow(dataset) + 1,] <- NA
	dataset$Species[nrow(dataset)] <- 'Palaeognathae'
	dataset$Body_mass_g[nrow(dataset)] <- 15700
	dataset$Brain_size_g[nrow(dataset)] <- 1.036 * 9.6
	
	# Neognathae
	dataset[nrow(dataset) + 1,] <- NA
	dataset$Species[nrow(dataset)] <- 'Neognathae'
	dataset$Body_mass_g[nrow(dataset)] <- 2900
	dataset$Brain_size_g[nrow(dataset)] <- 1.036 * 5.6
	dataset$BMR_Watt[nrow(dataset)] <- 5.35 * 2900 * 20.1 / 3600 / 3.44
	
	# Galloanserae
	dataset[nrow(dataset) + 1,] <- NA
	dataset$Species[nrow(dataset)] <- 'Galloanserae'
	dataset$Body_mass_g[nrow(dataset)] <- 3050
	dataset$Brain_size_g[nrow(dataset)] <- 1.036 * 5.9
	
	# Neoaves
	dataset[nrow(dataset) + 1,] <- NA
	dataset$Species[nrow(dataset)] <- 'Neoaves'
	dataset$Body_mass_g[nrow(dataset)] <- 1450
	dataset$Brain_size_g[nrow(dataset)] <- 1.036 * 5.4
	
	# Given that there is uncertainty in the body mass and aquatic affinity
	# of the last common ancestor of Amniota, create a data frame to 
	# store all the possible combinations. Also, store a different random 
	# seed for each chain.
	parameters_df <- data.frame(
		putative_body_mass = rep(NA, 30),
		putative_aquatic_affinity = rep(NA, 30),
		seed = rep(NA, 30)
	)	

	current_row <- 1
	for ( putative_body_mass in c(2, 500, 1000) )
	{
		for ( putative_aquatic_affinity in c('low', 'moderate') )
		{
			
			# For all combinations of body mass and aquatic affinity, 
			# run 5 independent chains.
			for ( seed in 1:5 )
			{
				parameters_df[current_row,] <- c(
					putative_body_mass, putative_aquatic_affinity,
					seed
				)
				
				current_row <- current_row + 1
			}
		}
	}
	class(parameters_df$putative_body_mass) <- 'numeric'
	class(parameters_df$seed) <- 'numeric'
	
	# For the last common ancestor of Amniota, assign the appropriate values 
	# for the current chain, based on the user-provided chain ID.
	dataset$Body_mass_g[Amniota_row] <- parameters_df$putative_body_mass[chain_id]
	dataset$BMR_Watt[Amniota_row] <- 2.08 * dataset$Body_mass_g[Amniota_row] * 20.1 / 3600 / 2.8 
	dataset$Aquatic_affinity[Amniota_row] <- parameters_df$putative_aquatic_affinity[chain_id]

	# For 5 species, we have information on BMR, brain mass, or maximum
	# longevity, but not body mass. Therefore, for those species, mass 
	# corrections cannot be readily performed, and the data would have 
	# to be dropped.
	#
	# To avoid this, for those 5 species, we extracted the median posterior 
	# estimate of body mass from the mass-uncorrected MCMCglmm fits 
	# (obtained with the fit_MCMCglmm_without_body_mass_corrections.R script).
	# We then treated those estimates as 'data', allowing us to apply the 
	# mass corrections and to not discard any real data.
	dataset$Body_mass_g[!is.na(dataset$Species) & dataset$Species == 'Myotis_brandtii'] <- 21.9989588668361
	dataset$Body_mass_g[!is.na(dataset$Species) & dataset$Species == 'Myotis_evotis'] <- 8.83362029392236
	dataset$Body_mass_g[!is.na(dataset$Species) & dataset$Species == 'Myotis_sodalis'] <- 13.1374727596513
	dataset$Body_mass_g[!is.na(dataset$Species) & dataset$Species == 'Tamias_sibiricus'] <- 99.9724005561761
	dataset$Body_mass_g[!is.na(dataset$Species) & dataset$Species == 'Aves'] <- 1824.27433779559

	# Apply body mass corrections using the allometric scaling exponents 
	# estimated from the variance-covariance matrix inferred with the 
	# fit_MCMCglmm_without_body_mass_corrections.R script.
	dataset$BMR_Watt_mass_corrected <- dataset$BMR_Watt / dataset$Body_mass_g^0.710
	dataset$Brain_size_g_mass_corrected <- dataset$Brain_size_g / dataset$Body_mass_g^0.678
	dataset$Max_longevity_years_mass_corrected <- dataset$Max_longevity_years / dataset$Body_mass_g^0.164
	
	# Return all the necessary variables for model fitting.
	return(
		list(
			dataset = dataset, tree = tree, seed = parameters_df$seed[chain_id]
		)
	)
}

############################
# M  A  I  N    C  O  D  E #
############################

# Read the chain ID provided by the user as a command line argument.
args <- commandArgs(TRUE)
chain_id <- as.numeric(args[1])

# Prepare all the variables needed for fitting the model.
vars_for_MCMCglmm <- prepare_vars_for_MCMCglmm(chain_id)

# Set the random seed and fit the model.
set.seed(vars_for_MCMCglmm$seed)	
fit <- MCMCglmm(

	# Define the response variables, apply any needed transformations, 
	# and specify a distinct intercept per response.
	cbind(
		log(Body_mass_g), log(BMR_Watt_mass_corrected), log(Brain_size_g_mass_corrected), 
		log(Max_longevity_years_mass_corrected), I(Range_size_km2^(1/5)), sqrt(Absolute_latitude),
		Mean_temp, I(SD_temp^(1/3)), log(Annual_precip), I(CV_precip^(1/4)), I(NPP^(1/3)),
		Torpor, Migratory, Carnivory, Herbivory, Fossoriality, Aquatic_affinity, 
		Hemisphere, Cathemeral, Crepuscular, Diurnal, Nocturnal
	) ~ trait - 1,
	
	# Specify a phylogenetic random effect on each intercept.
	random =~ us(trait):Species,

	# Set the distribution for each response variable.
	family = c(
		'gaussian', 'gaussian', 'gaussian',
		'gaussian', 'gaussian', 'gaussian',
		'gaussian', 'gaussian', 'gaussian', 'gaussian', 'gaussian',
		'threshold', 'threshold', 'threshold', 'threshold', 'threshold', 'threshold',
		'threshold', 'threshold', 'threshold', 'threshold', 'threshold'
	),
	
	# Integrate the phylogenetic variance-covariance matrix into the model.
	ginverse = list(Species=inverseA(vars_for_MCMCglmm$tree, nodes = 'ALL', scale = TRUE)$Ainv),

	# Set relatively uninformative priors.
	prior = list(
		G = list(
			G1 = list(
				V=diag(22), nu=50, alpha.mu=rep(0,22), 
				alpha.V=diag(rep(1000,22))
			)
		),
		R = list(V = diag(22), nu = 50, fix = 22)
	),
	
	# Set the data frame with all the needed data.
	data = vars_for_MCMCglmm$dataset,
	
	# Allow the model to estimate the covariances among response variables.
	rcov =~ us(trait):units,

	# Set the number of iterations, the burn-in, and the sampling frequency.
	nitt = 1000000,
	burnin = 100000,
	thin = 25,
	verbose = TRUE,

	# Store the posterior distributions of random effects and latent variables.
	pr = TRUE,
	pl = TRUE,
	
	# Force the threshold liabilities to range from -7 to 7 to facilitate 
	# their estimation.
	trunc = TRUE
)

# Wait for 10 seconds after fitting has finished (just in case), create 
# an output directory if missing, and save the resulting fit to an .Rda 
# file.
Sys.sleep(10)
dir.create('../Results/MCMCglmm_fits_with_mass_corrections/', showWarnings = FALSE)
save(
	fit, file = paste(
		'../Results/MCMCglmm_fits_with_mass_corrections/', chain_id, 
		'.Rda', sep = ''
	)
)
