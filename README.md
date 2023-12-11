This repository contains code that reproduces the main analyses of the following study:

>Dimitrios - Georgios Kontopoulos, Danielle L. Levesque, and Michael Hiller: **Numerous independent gains of torpor and hibernation across endotherms, linked with adaptation to diverse environments**. 2023. [Submitted]

---
 
#### Execution

To run the scripts in this repository, please first create "Data/" and "Results/" directories, outside "Code/". 
Inside Data/, place the files deposited in [this Figshare repository](https://doi.org/10.6084/m9.figshare.24746283).

Most scripts need to be run from the command line with some user-provided arguments. Specific instructions 
are provided at the top of each script. Note that some scripts will take more than a week to finish running and 
will require high amounts of memory (typically around 30-70 GB).

Given that some scripts depend on the outputs of other scripts, you can follow the order below 
as a rough guide (but also read the information at the top of each script):

1. fit_Mk_variants.R
2. fit_MCMCglmm_without_body_mass_corrections.R
3. fit_MCMCglmm_without_body_mass_corrections_only_across_mammals.R
4. fit_MCMCglmm_without_body_mass_corrections_only_across_birds.R
5. extract_allometric_scaling_exponents.R
6. fit_MCMCglmm_with_body_mass_corrections.R
7. fit_MCMCglmm_with_body_mass_corrections_only_across_mammals.R
8. fit_MCMCglmm_with_body_mass_corrections_only_across_birds.R
9. extract_threshold_between_daily_torpor_and_hibernation.R
10. extract_correlations.R
11. run_phylogenetic_pca.R
12. fit_Mk_variants_with_some_dormancy_incapable_mammals_shifted_to_dormancy.R
13. fit_MCMCglmm_with_body_mass_corrections_with_some_dormancy_incapable_mammals_shifted_to_dormancy.R
14. check_ESS_PSRF.R
15. extract_dormancy_probabilities.R
