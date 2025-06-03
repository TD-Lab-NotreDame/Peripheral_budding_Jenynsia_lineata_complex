library(OptM)

#need to create a folder with all output files from treemix
folder <- file.path(path="data/tree")                     #path to files of TreeMix replicates with different migration edges (m) to test       
#shows changes in log likelihood for different m, and suggests optimum m as 'change points'

test.optM = optM(folder, tsv ="Evanno.variance.txt")  #another option is the Evanno method - see optM package description for detailed information on output
#if data is robust and all runs have the same likelihoods, SD will be 0 and this function will give an error as it can't produce the ad hoc statistic. 
plot_optM(test.optM, method = "Evanno")                          
#plot the proportion of variation explained by each migration event. Calculates deltaM, which is a second-order rate of change in likelihood weighted by the standard deviation

ggsave("optm_migration.pdf", plot = plot, bg = "transparent", width = 3.2, height = 3)
