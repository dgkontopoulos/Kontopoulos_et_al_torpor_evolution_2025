#!/usr/bin/env Rscript

# This script fits 12 variants of the Mk model to the dormancy data 
# and performs model averaging on the basis of AIC weights.
#
# Based on the result, it then reconstructs the ancestral states of 
# dormancy for all internal nodes and identifies branches in which 
# direct shifts between no dormancy and hibernation have occurred.

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
fit_models <- function(tree, dormancy)
{
	
	# Fit each Mk variant separately and calculate the AIC.
	
	# All transition rates are equal.
	cat('Now fitting the Mk_ER variant...\n')
	set.seed(1)
	fit_Mk_ER <- fitMk(tree, dormancy, model = 'ER', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ER <- AIC(fit_Mk_ER)
	
	# Forward and backward transition rates are equal.
	cat('Now fitting the Mk_SYM variant...\n')
	set.seed(2)
	fit_Mk_SYM <- fitMk(tree, dormancy, model = 'SYM', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_SYM <- AIC(fit_Mk_SYM)
	
	# All transition rates can differ from each other.
	cat('Now fitting the Mk_ARD variant...\n')
	set.seed(3)
	fit_Mk_ARD <- fitMk(tree, dormancy, model = 'ARD', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ARD <- AIC(fit_Mk_ARD)
	
	# For Mk variants with the lambda parameter, we need to first 
	# determine the optimum lambda value for each variant.
	
	# Equal transition rates + lambda transformation.
	cat('Now fitting the Mk_ER_lambda variant...\n')
	set.seed(4)
	lambda_Mk_ER <- determine_optimum_lambda(tree, dormancy, model = 'ER')
	fit_Mk_ER_lambda <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_ER), dormancy, model = 'ER', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ER_lambda <- -2 * fit_Mk_ER_lambda$logLik + 2 * (attr(logLik(fit_Mk_ER_lambda), 'df') + 1)
	
	# Equal forward/backward transition rates + lambda transformation.
	cat('Now fitting the Mk_SYM_lambda variant...\n')
	set.seed(5)
	lambda_Mk_SYM <- determine_optimum_lambda(tree, dormancy, model = 'SYM')
	fit_Mk_SYM_lambda <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_SYM), dormancy, model = 'SYM', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_SYM_lambda <- -2 * fit_Mk_SYM_lambda$logLik + 2 * (attr(logLik(fit_Mk_SYM_lambda), 'df') + 1)
	
	# Different transition rates + lambda transformation.
	cat('Now fitting the Mk_ARD_lambda variant...\n')
	set.seed(6)
	lambda_Mk_ARD <- determine_optimum_lambda(tree, dormancy, model = 'ARD')
	fit_Mk_ARD_lambda <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_ARD), dormancy, model = 'ARD', pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ARD_lambda <- -2 * fit_Mk_ARD_lambda$logLik + 2 * (attr(logLik(fit_Mk_ARD_lambda), 'df') + 1)

	# For Mk variants that do not allow direct shifts between no dormancy 
	# and hibernation, we need to provide a custom transition rate matrix.
	
	# Equal transition rates + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_ER_ordered variant...\n')
	set.seed(7)
	matrix_ER <-matrix(c(0,0,1,0,0,1,1,1,0),3,3,byrow=TRUE,
	    dimnames=list(c("Hibernation","NO","Torpor"),c("Hibernation","NO","Torpor")))
	fit_Mk_ER_ordered <- fitMk(tree, dormancy, model = matrix_ER, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ER_ordered <- AIC(fit_Mk_ER_ordered)
	
	# Equal forward/backward transition rates + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_SYM_ordered variant...\n')
	set.seed(8)
	matrix_SYM <-matrix(c(0,0,1,0,0,2,1,2,0),3,3,byrow=TRUE,
	    dimnames=list(c("Hibernation","NO","Torpor"),c("Hibernation","NO","Torpor")))
	fit_Mk_SYM_ordered <- fitMk(tree, dormancy, model = matrix_SYM, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_SYM_ordered <- AIC(fit_Mk_SYM_ordered)
	
	# Different transition rates + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_ARD_ordered variant...\n')
	set.seed(9)
	matrix_ARD <-matrix(c(0,0,1,0,0,2,3,4,0),3,3,byrow=TRUE,
	    dimnames=list(c("Hibernation","NO","Torpor"),c("Hibernation","NO","Torpor")))
	fit_Mk_ARD_ordered <- fitMk(tree, dormancy, model = matrix_ARD, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ARD_ordered <- AIC(fit_Mk_ARD_ordered)
	
	# Equal transition rates + lambda transformation + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_ER_lambda_ordered variant...\n')
	set.seed(10)
	lambda_Mk_ER_ordered <- determine_optimum_lambda(tree, dormancy, model = matrix_ER)
	fit_Mk_ER_lambda_ordered <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_ER_ordered), dormancy, model = matrix_ER, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_ER_lambda_ordered <- -2 * fit_Mk_ER_lambda_ordered$logLik + 2 * (attr(logLik(fit_Mk_ER_lambda_ordered), 'df') + 1)
	
	# Equal forward/backward transition rates + lambda transformation + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_SYM_lambda_ordered variant...\n')
	set.seed(11)
	lambda_Mk_SYM_ordered <- determine_optimum_lambda(tree, dormancy, model = matrix_SYM)
	fit_Mk_SYM_lambda_ordered <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_SYM_ordered), dormancy, model = matrix_SYM, pi='fitzjohn', lik.func="pruning")
	AIC_fit_Mk_SYM_lambda_ordered <- -2 * fit_Mk_SYM_lambda_ordered$logLik + 2 * (attr(logLik(fit_Mk_SYM_lambda_ordered), 'df') + 1)
	
	# Different transition rates + lambda transformation + lack of shifts between no dormancy and hibernation.
	cat('Now fitting the Mk_ARD_lambda_ordered variant...\n')
	set.seed(12)
	lambda_Mk_ARD_ordered <- determine_optimum_lambda(tree, dormancy, model = matrix_ARD)
	fit_Mk_ARD_lambda_ordered <- fitMk(phytools:::lambdaTree(tree, lambda_Mk_ARD_ordered), dormancy, model = matrix_ARD, pi='fitzjohn', lik.func="pruning")
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
		file = '../Results/Mk_fits.Rda'
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
# the phylogeny and detects apparent direct transitions between no 
# dormancy and hibernation.
reconstruct_ancestral_states <- function(tree, dormancy)
{
	
	# Load the Mk model fitting results.
	load('../Results/Mk_fits.Rda')
	
	# Perform 10,000 stochastic character mapping simulations based on 
	# the model-averaged parameters.
	simulated_trees <- make.simmap(
		tree, x = dormancy, Q = results$averaged_Q_matrix,
		pi = results$averaged_pi, nsim = 10000
	)
	
	# Get the dormancy probabilities for each node and write them to 
	# an output file.
	model_estimates <- describe.simmap(simulated_trees)$ace
	
	write.csv(
		model_estimates,
		file = '../Results/dormancy_probabilities_Mk.csv'
	)
	
	# Initialise a data frame to store the frequency of direct transitions 
	# between no dormancy and hibernation for each branch across all 
	# simulations.
	transitions_between_NO_and_Hibernation <- data.frame(
		edge = 1:2625,
		NO_to_Hibernation = 0,
		Hibernation_to_NO = 0
	)
	
	# For each simulation ...
	for ( i in 1:10000 )
	{
		
		# ... and for each branch ...
		for ( j in 1:2625 )
		{
			
			# ... check if there are at least two observed dormancy states.
			if ( length(simulated_trees[[i]]$maps[[j]]) >= 2 )
			{
				
				# Read the dormancy states along the branch and check if 
				# a direct shift from no dormancy to hibernation has 
				# occurred.
				#
				# If yes, increase the relevant counter in the data frame.
				
				k <- 2
				while ( k <= length(simulated_trees[[i]]$maps[[j]]) )
				{
					if (
						names(simulated_trees[[i]]$maps[[j]])[k] == 'Hibernation' && 
						names(simulated_trees[[i]]$maps[[j]])[k - 1] == 'NO'
					)
					{
						transitions_between_NO_and_Hibernation$NO_to_Hibernation[j] <- transitions_between_NO_and_Hibernation$NO_to_Hibernation[j] + 1
						break
					}
					k <- k + 1
				}
			}
		}
	}
	
	# The lines below do the exact same thing as above, but for direct 
	# shifts from hibernation to lack of dormancy.
	for ( i in 1:10000 )
	{
		for ( j in 1:2625 )
		{
			if ( length(simulated_trees[[i]]$maps[[j]]) >= 2 )
			{
				k <- 2
				while ( k <= length(simulated_trees[[i]]$maps[[j]]) )
				{
					if (
						names(simulated_trees[[i]]$maps[[j]])[k] == 'NO' && 
						names(simulated_trees[[i]]$maps[[j]])[k - 1] == 'Hibernation'
					)
					{
						transitions_between_NO_and_Hibernation$Hibernation_to_NO[j] <- transitions_between_NO_and_Hibernation$Hibernation_to_NO[j] + 1
						break
					}
					k <- k + 1
				}
			}
		}
	}
	
	# Divide the number of direct shifts per branch by the total number 
	# of simulations.	
	transitions_between_NO_and_Hibernation$NO_to_Hibernation <- transitions_between_NO_and_Hibernation$NO_to_Hibernation / 10000
	transitions_between_NO_and_Hibernation$Hibernation_to_NO <- transitions_between_NO_and_Hibernation$Hibernation_to_NO / 10000
	
	# Find branches where a direct shift occurred in at least half of all
	# simulations.
	key_edges_NO_to_Hibernation <- which(
		transitions_between_NO_and_Hibernation$NO_to_Hibernation > 0.5
	)
	key_edges_Hibernation_to_NO <- which(
		transitions_between_NO_and_Hibernation$Hibernation_to_NO > 0.5
	)
	
	transitions_between_NO_and_Hibernation <- transitions_between_NO_and_Hibernation[
		c(key_edges_NO_to_Hibernation, key_edges_Hibernation_to_NO),
	]
	
	# Write the results to an output file.
	write.csv(
		transitions_between_NO_and_Hibernation, 
		file = '../Results/apparent_direct_transitions.csv',
		row.names = FALSE
	)
}

############################
# M  A  I  N    C  O  D  E #
############################

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
fit_models(tree, dormancy)

# Perform stochastic character mapping simulations to reconstruct ancestral
# states and detect direct shifts between no dormancy and hibernation.
reconstruct_ancestral_states(tree, dormancy)
