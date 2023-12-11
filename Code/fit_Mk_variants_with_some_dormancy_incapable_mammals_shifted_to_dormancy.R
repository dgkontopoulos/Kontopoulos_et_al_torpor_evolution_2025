#!/usr/bin/env Rscript

# This script fits 12 variants of the Mk model to the dormancy data 
# and performs model averaging on the basis of AIC weights.
#
# Based on the result, it then reconstructs the ancestral states of 
# dormancy for all internal nodes.
#
# Furthermore, in this script, some dormancy-incapable mammals are 
# randomly shifted to either daily torpor or hibernation. 
#
# The script needs to be run from the command line, with the user 
# providing the number of dormancy-incapable mammals that are randomly 
# shifted to dormancy, as well as the random seed, e.g.:
#
# Rscript fit_Mk_variants_with_some_dormancy_incapable_mammals_shifted_to_dormancy.R 20 1

library(geiger)
library(phytools)

#############################
# F  U  N  C  T  I  O  N  S #
#############################

# This function determines the optimum value of Pagel's lambda for a
# variant of the Mk model.
#
# Modified from: http://blog.phytools.org/2018/02/how-to-fit-tree-transformation-for.html
determine_optimum_lambda <- function(tree, x, model)
{
	opt <- optimize(
		lk.lambda, c(0,phytools:::maxLambda(tree)), tree = tree,
		x = x, model = model
	)
	
	return(opt$minimum)
}

# This function fits all 12 Mk model variants and performs model averaging.
fit_models <- function(tree, dormancy_vector, random_shifts, user_seed)
{

	# Get all placental mammals in the dataset.	
	placental_mammals <- as.character(
		na.omit(
			tree$tip.label[
				getDescendants(
					tree, getMRCA(
						tree, c(
							'Tamandua_tetradactyla', 
							'Ictidomys_tridecemlineatus'
						)
					)
				)
			]
		)
	)

	# Find the indices of placental mammals that are listed as 
	# dormancy-incapable in our dataset.
	placental_mammals_without_dormancy <- which(
		names(dormancy_vector) %in% placental_mammals &
		dormancy_vector == 'NO'
	)

	# Randomly shift some of those (number set by the user) to either 
	# daily torpor or hibernation.
	set.seed(user_seed)
	dormancy_vector[
		sample(placental_mammals_without_dormancy, random_shifts)
	] <- sample(c('Torpor', 'Hibernation'), random_shifts, replace = TRUE)
	
	# Fit each Mk variant separately and calculate the AIC.
	
	# All transition rates are equal.
	cat('Now fitting the Mk_ER variant...\n')
	set.seed(1)
	fit_Mk_ER <- fitMk(tree, dormancy_vector, model = 'ER', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ER <- AIC(fit_Mk_ER)
	
	# Forward and backward transition rates are equal.
	cat('Now fitting the Mk_SYM variant...\n')
	set.seed(2)
	fit_Mk_SYM <- fitMk(tree, dormancy_vector, model = 'SYM', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_SYM <- AIC(fit_Mk_SYM)
	
	# All transition rates can differ from each other.
	cat('Now fitting the Mk_ARD variant...\n')
	set.seed(3)
	fit_Mk_ARD <- fitMk(tree, dormancy_vector, model = 'ARD', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ARD <- AIC(fit_Mk_ARD)
	
	# For Mk variants with the lambda parameter, we need to first 
	# determine the optimum lambda value for each variant.
	
	# Equal transition rates + lambda transformation.
	cat('Now fitting the Mk_ER_lambda variant...\n')
	set.seed(4)
	lambda_Mk_ER <- determine_optimum_lambda(tree, dormancy_vector, model = 'ER')
	fit_Mk_ER_lambda <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_ER), dormancy_vector, model = 'ER', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ER_lambda <- -2 * fit_Mk_ER_lambda$logLik + 2 * (attr(logLik(fit_Mk_ER_lambda), 'df') + 1)
	
	# Equal forward/backward transition rates + lambda transformation.
	cat('Now fitting the Mk_SYM_lambda variant...\n')
	set.seed(5)
	lambda_Mk_SYM <- determine_optimum_lambda(tree, dormancy_vector, model = 'SYM')
	fit_Mk_SYM_lambda <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_SYM), dormancy_vector, model = 'SYM', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_SYM_lambda <- -2 * fit_Mk_SYM_lambda$logLik + 2 * (attr(logLik(fit_Mk_SYM_lambda), 'df') + 1)
	
	# Different transition rates + lambda transformation.
	cat('Now fitting the Mk_ARD_lambda variant...\n')
	set.seed(6)
	lambda_Mk_ARD <- determine_optimum_lambda(tree, dormancy_vector, model = 'ARD')
	fit_Mk_ARD_lambda <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_ARD), dormancy_vector, model = 'ARD', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ARD_lambda <- -2 * fit_Mk_ARD_lambda$logLik + 2 * (attr(logLik(fit_Mk_ARD_lambda), 'df') + 1)

	# For Mk variants that do not allow direct shifts between no dormancy 
	# and hibernation, we need to provide a custom transition rate matrix.
	
	# Equal transition rates + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_ER_ordered variant...\n')
	set.seed(7)
	matrix_ER <-matrix(c(0,0,1,0,0,1,1,1,0),3,3,byrow=TRUE,
	    dimnames=list(c("Hibernation","NO","Torpor"),c("Hibernation","NO","Torpor")))
	fit_Mk_ER_ordered <- fitMk(tree, dormancy_vector, model = matrix_ER, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ER_ordered <- AIC(fit_Mk_ER_ordered)
	
	# Equal forward/backward transition rates + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_SYM_ordered variant...\n')
	set.seed(8)
	matrix_SYM <-matrix(c(0,0,1,0,0,2,1,2,0),3,3,byrow=TRUE,
	    dimnames=list(c("Hibernation","NO","Torpor"),c("Hibernation","NO","Torpor")))
	fit_Mk_SYM_ordered <- fitMk(tree, dormancy_vector, model = matrix_SYM, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_SYM_ordered <- AIC(fit_Mk_SYM_ordered)
	
	# Different transition rates + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_ARD_ordered variant...\n')
	set.seed(9)
	matrix_ARD <-matrix(c(0,0,1,0,0,2,3,4,0),3,3,byrow=TRUE,
	    dimnames=list(c("Hibernation","NO","Torpor"),c("Hibernation","NO","Torpor")))
	fit_Mk_ARD_ordered <- fitMk(tree, dormancy_vector, model = matrix_ARD, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ARD_ordered <- AIC(fit_Mk_ARD_ordered)
	
	# Equal transition rates + lambda transformation + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_ER_lambda_ordered variant...\n')
	set.seed(10)
	lambda_Mk_ER_ordered <- determine_optimum_lambda(tree, dormancy_vector, model = matrix_ER)
	fit_Mk_ER_lambda_ordered <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_ER_ordered), dormancy_vector, model = matrix_ER, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ER_lambda_ordered <- -2 * fit_Mk_ER_lambda_ordered$logLik + 2 * (attr(logLik(fit_Mk_ER_lambda_ordered), 'df') + 1)
	
	# Equal forward/backward transition rates + lambda transformation + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_SYM_lambda_ordered variant...\n')
	set.seed(11)
	lambda_Mk_SYM_ordered <- determine_optimum_lambda(tree, dormancy_vector, model = matrix_SYM)
	fit_Mk_SYM_lambda_ordered <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_SYM_ordered), dormancy_vector, model = matrix_SYM, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_SYM_lambda_ordered <- -2 * fit_Mk_SYM_lambda_ordered$logLik + 2 * (attr(logLik(fit_Mk_SYM_lambda_ordered), 'df') + 1)
	
	# Different transition rates + lambda transformation + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_ARD_lambda_ordered variant...\n')
	set.seed(12)
	lambda_Mk_ARD_ordered <- determine_optimum_lambda(tree, dormancy_vector, model = matrix_ARD)
	fit_Mk_ARD_lambda_ordered <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_ARD_ordered), dormancy_vector, model = matrix_ARD, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ARD_lambda_ordered <- -2 * fit_Mk_ARD_lambda_ordered$logLik + 2 * (attr(logLik(fit_Mk_ARD_lambda_ordered), 'df') + 1)
	
	cat('Calculating AIC weights and doing model averaging...\n')
	
	# Calculate the AIC weight of each Mk variant.
	AIC_weights <- aicw(
		c(
			AIC_fit_Mk_ER,					AIC_fit_Mk_SYM,					AIC_fit_Mk_ARD,
			AIC_fit_Mk_ER_lambda,			AIC_fit_Mk_SYM_lambda,			AIC_fit_Mk_ARD_lambda,
			AIC_fit_Mk_ER_ordered,			AIC_fit_Mk_SYM_ordered,			AIC_fit_Mk_ARD_ordered,
			AIC_fit_Mk_ER_lambda_ordered,	AIC_fit_Mk_SYM_lambda_ordered,	AIC_fit_Mk_ARD_lambda_ordered
		)
	)
	
	# Perform model averaging.
	averaged_Q_matrix <-
		as.Qmatrix(fit_Mk_ER) * AIC_weights$w[1] +
		as.Qmatrix(fit_Mk_SYM) * AIC_weights$w[2] +
		as.Qmatrix(fit_Mk_ARD) * AIC_weights$w[3] +
		as.Qmatrix(fit_Mk_ER_lambda) * AIC_weights$w[4] +
		as.Qmatrix(fit_Mk_SYM_lambda) * AIC_weights$w[5] +
		as.Qmatrix(fit_Mk_ARD_lambda) * AIC_weights$w[6] +
		as.Qmatrix(fit_Mk_ER_ordered) * AIC_weights$w[7] +
		as.Qmatrix(fit_Mk_SYM_ordered) * AIC_weights$w[8] +
		as.Qmatrix(fit_Mk_ARD_ordered) * AIC_weights$w[9] +
		as.Qmatrix(fit_Mk_ER_lambda_ordered) * AIC_weights$w[10] +
		as.Qmatrix(fit_Mk_SYM_lambda_ordered) * AIC_weights$w[11] +
		as.Qmatrix(fit_Mk_ARD_lambda_ordered) * AIC_weights$w[12]
	class(averaged_Q_matrix) <- 'matrix'
	
	averaged_pi <-
		fit_Mk_ER$pi * AIC_weights$w[1] +
		fit_Mk_SYM$pi * AIC_weights$w[2] +
		fit_Mk_ARD$pi * AIC_weights$w[3] +
		fit_Mk_ER_lambda$pi * AIC_weights$w[4] +
		fit_Mk_SYM_lambda$pi * AIC_weights$w[5] +
		fit_Mk_ARD_lambda$pi * AIC_weights$w[6] +
		fit_Mk_ER_ordered$pi * AIC_weights$w[7] +
		fit_Mk_SYM_ordered$pi * AIC_weights$w[8] +
		fit_Mk_ARD_ordered$pi * AIC_weights$w[9] +
		fit_Mk_ER_lambda_ordered$pi * AIC_weights$w[10] +
		fit_Mk_SYM_lambda_ordered$pi * AIC_weights$w[11] +
		fit_Mk_ARD_lambda_ordered$pi * AIC_weights$w[12]
	
	# Store the results to an output file.
	results <- mget(ls())
	save(
		list = 'results', 
		file = paste(
			'../Results/Mk_fits_', random_shifts, '_random_shifts_seed_', 
			user_seed, '.Rda', sep = ''
		)
	)
}

# This function returns the log-likelihood of a given Mk model with 
# a lambda transformation.
#
# Copied from: http://blog.phytools.org/2018/02/how-to-fit-tree-transformation-for.html
lk.lambda <- function(lambda, tree, x, ...) 
{
	-logLik(
		fitMk(
			phytools:::lambdaTree(tree, lambda), pi = 'fitzjohn',
			x, ...
        )
	)
}

# This function reconstructs the ancestral states of dormancy across 
# the phylogeny.
reconstruct_ancestral_states <- function(tree, random_shifts, user_seed)
{
	
	# Load the Mk model fitting results.
	load(
		paste(
			'../Results/Mk_fits_', random_shifts, '_random_shifts_seed_', 
			user_seed, '.Rda', sep = ''
		)
	)
	
	# Perform 10,000 stochastic character mapping simulations based on 
	# the model-averaged parameters.
	simulated_trees <- make.simmap(
		tree, x = results$dormancy_vector, Q = results$averaged_Q_matrix,
		pi = results$averaged_pi, nsim = 10000
	)
	
	# Get the dormancy probabilities for each node and write them to 
	# an output file.
	model_estimates <- describe.simmap(simulated_trees)$ace
	
	write.csv(
		model_estimates,
		file = paste(
			'../Results/dormancy_probabilities_Mk_', random_shifts, 
			'_random_shifts_seed_', user_seed, '.csv', sep = ''
		)
	)
}

############################
# M  A  I  N    C  O  D  E #
############################

# Read the number of random shifts and the seed provided by the user as 
# command line arguments.
args <- commandArgs(TRUE)
random_shifts <- as.numeric(args[1])
user_seed <- as.numeric(args[2])

# Read the dataset, and replace spaces in species' names with underscores.
dataset <- read.csv('../Data/dataset.csv')
dataset$Species <- gsub(' ', '_', dataset$Species)

# Read the phylogeny.
tree <- read.tree('../Data/time_calibrated_phylogeny.nwk')
tree$node.label <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)

# Prepare a dormancy factor for Mk model fitting.
dormancy <- factor(dataset$Dormancy)
names(dormancy) <- dataset$Species

# Fit all Mk model variants.
fit_models(tree, dormancy, random_shifts, user_seed)

# Perform stochastic character mapping simulations to reconstruct ancestral
# states and detect direct shifts between no dormancy and hibernation.
reconstruct_ancestral_states(tree, random_shifts, user_seed)
