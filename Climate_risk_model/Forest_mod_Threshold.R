locs <- data.frame(cbind(species.pa[,3], species.pa[,2], species.pa$ALL))
names(locs) <- c("Lat", "Long", "PA")
locs$pred <- extract(predictions[[1]], locs[,1:2])

pred <- prediction(locs$pred, locs$PA)               
perf <- performance(pred, "tpr", "fpr")
fpr <- perf@x.values[[1]]
tpr <- perf@y.values[[1]]
sum.v <- tpr + (1-fpr)
index <- which.max(sum.v)
cutoff <- perf@alpha.values[[1]][[index]]

### Used equal sensitivity and specificity instead
threshold(forest.eval.1)




