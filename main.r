classes = c("numeric", "numeric", "numeric", "numeric")
dioecious <- read.csv("dioecious", header=TRUE, colClasses=classes)

twowayanova <- aov(mating ~ Dmin + Dmax, data = dioecious)

summary(twowayanova)
