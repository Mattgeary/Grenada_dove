f.50 <- future(forest.mod, years = 35, reps = 1000, threshold = 0.66)
save(f.50, file = "f.50.RData")
f.70 <- future(forest.mod, years = 55, reps = 1000, threshold = 0.66)
save(f.70, file = "f.70.RData")
