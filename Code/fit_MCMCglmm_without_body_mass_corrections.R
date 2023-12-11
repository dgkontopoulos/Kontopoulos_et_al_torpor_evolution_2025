#!/usr/bin/env Rscript

# This script fits a model with MCMCglmm to estimate the correlation 
# structure among dormancy and 21 ecophysiological variables across 
# 1,338 endotherms, accounting for phylogeny.
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
# Rscript fit_MCMCglmm_without_body_mass_corrections.R 1
# Rscript fit_MCMCglmm_without_body_mass_corrections.R 2
# ...
# Rscript fit_MCMCglmm_without_body_mass_corrections.R 29
# Rscript fit_MCMCglmm_without_body_mass_corrections.R 30
#
# This script will run for a few days depending on the CPU (~7.5 days on 
# an AMD EPYC 7702 CPU core) and will need ~67 GB of memory. In the end, 
# it will produce an .Rda file per chain, in the ../Results/MCMCglmm_fits/ 
# directory.
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

	dataset$Dormancy <- factor(
		dataset$Dormancy, levels = c('NO', 'Torpor', 'Hibernation')
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
		log(Body_mass_g), log(BMR_Watt), log(Brain_size_g), 
		log(Max_longevity_years), I(Range_size_km2^(1/5)), sqrt(Absolute_latitude),
		Mean_temp, I(SD_temp^(1/3)), log(Annual_precip), I(CV_precip^(1/4)), I(NPP^(1/3)),
		Dormancy, Migratory, Carnivory, Herbivory, Fossoriality, Aquatic_affinity, 
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
dir.create('../Results/MCMCglmm_fits/', showWarnings = FALSE)
save(fit, file = paste('../Results/MCMCglmm_fits/', chain_id, '.Rda', sep = ''))
