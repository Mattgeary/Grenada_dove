vals.clim.50 <- data.frame("Scenario" = c(names(f.50[[1]]), paste(names(f.50$climate), "2050", sep = "_")), "Predicted_area" = c(as.data.frame(f.50[[1]])[,1], as.data.frame(t(f.50$climate))[,1]))

vals.clim.70 <- data.frame("Scenario" = c(names(f.70[[1]]), paste(names(f.70$climate), "2070", sep = "_")), "Predicted_area" = c(as.data.frame(f.70[[1]])[,1], as.data.frame(t(f.70$climate))[,1]))

climate.results <- rbind(vals.clim.50, vals.clim.70[-c(1:3),])

write.csv(climate.results, "Climate_results.csv")
