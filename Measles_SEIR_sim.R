library("prob")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("graphics")
source("SEIR_func.R")
num.steps = 10000
source("Initial_conditions.R")


source("SEIR_measles.R")

disease.state                 <-      initial.disease.state (demographic.ages  ,  v  , 1 ,  num.comps)
disease.state                 =       reduce.susceptibles (0, 5, disease.state, 0.92, num.comps, susceptible.indices)
disease.state                 =       reduce.susceptibles (6, 20, disease.state, 0.92, num.comps, susceptible.indices)
disease.state                 =       reduce.susceptibles (21, max(demographic.ages[, 1]), disease.state, 0.95, num.comps, susceptible.indices)
pop.by.age = number.of.each.age(demographic.ages,disease.state, num.comps)

plot(disease.state[susceptible.indices]/pop.by.age)
a = Run.Sims(num.steps)

disease.state = unlist(a[1])
infections.per.year = unlist(a[2])
