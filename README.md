# AraNetBench
Evaluates an inferred regulatory network against state of the art interaction databases in Arabidopsis thaliana

The data used to validate the network : 

+ ConnecTF
+ AtRegNet
+ Litterature data

The available functions :

+ flatten edges (in case of grouped nodes)
+ evaluateNetwork (gives TP, FP, TPR, precision, recall, accuracy, for all databases or for a specific combination of databases)
+ evaluateWeights (gives AUC, AUPR for the weigts for all databases or for a specific combination of databases)
+ testAgainstRandom (computes null distribution of TPR in randomized networks and test observed TPR to return a pvalue)
+ Plot validated network (with edges colored depending on their validation status

Required input : a dataframe of inferred edges for a network, a matrix of weights for evaluateWeights
