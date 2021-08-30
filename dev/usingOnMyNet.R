

# network from current ongoing work on CO2*N
load("~/Documents/Combi/Results/network_N_CO2.RData")

net <- network_N_CO2$network_data$edges

res <- evaluate_network(net)

draw_evaluated_network(res)

test_validation_rate(net)
