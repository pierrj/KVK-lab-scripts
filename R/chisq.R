

observed1 <- 684
total1 <- 55318
expected1 <- 0.0100369
observed1_vector <- c(observed1, total1-observed1)
expected1_vector <- c(expected1,1-expected1)
chisq.test(observed1_vector, p = expected1_vector)

observed2 <- 323
total2 <- 30836
expected2 <- 0.00877053
observed2_vector <- c(observed2, total2-observed2)
expected2_vector <- c(expected2,1-expected2)
chisq.test(observed2_vector, p = expected2_vector)

observed3 <- 301
total3 <- 25650
expected3 <- 0.00806683
observed3_vector <- c(observed3, total3-observed3)
expected3_vector <- c(expected3,1-expected3)
chisq.test(observed3_vector, p = expected3_vector)


### events are the number of eccDNAs
### proportion is those with and those without

### expected is if you distribute them randomly what proportion do you get