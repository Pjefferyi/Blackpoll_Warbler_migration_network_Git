#plot solar decliantion 
solar_data <- solar(twl$Twilight)

plot(twl$Twilight, solar_data$sinSolarDec, col = "black")

abline(v = c(fall.equi, spring.equi), lty = 2, col = "orange", lwd = 2)

abline(h = c(0.25, -0.25), lty = 2, col = "blue")
abline(h = c(0.18, -0.18), lty = 2, col = "red")
abline(h = c(0.10, -0.10), lty = 2, col = "green")
abline(h = c(0.05, -0.05), lty = 2, col = "purple")

grid()
legend(legend = c("tol = 0.25", "tol = 0.18", "tol = 0.10","tol = 0.05"),
       col = c("blue", "red", "green", "purple"), lty = 2, "topright")
