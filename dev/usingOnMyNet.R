library(AraNetBench)

# network from current ongoing work on CO2*N
load("~/Documents/Combi/Results/network_N_CO2.RData")


load("D:/These/Combi/Results/network_N_CO2.RData")
net <- network_N_CO2$network_data$edges

res <- evaluate_network(net, 
                        validation = c("TARGET"))


res <- evaluate_network(net, 
                        validation = c("CHIPSeq"))

res <- evaluate_network(net, 
                        validation = c("DAPSeq"))


draw_evaluated_network(res)


test_validation_rate(net, validation = c("TARGET"))
