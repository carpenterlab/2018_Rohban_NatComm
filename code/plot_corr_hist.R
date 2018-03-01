x1 <- (cr.1.2.mean)
x2 <- (cr.1.2.mad)

x <- setdiff(x1, 1)
y <- setdiff(x2, 1)

D <- rbind(data.frame(method = "cov", x = y), data.frame(method = "med.", x = x))
ggplot(D, aes(x = x, group = method, fill = method)) + geom_histogram(aes(y = ..density..), position="identity", alpha = 0.8, bins = 100)
