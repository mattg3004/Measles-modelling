library("prob")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("graphics")
source("SEIR_func.R")
num.steps = 10000
time.step = 5
source("Initial_conditions.R")


source("SEIR_measles.R")
source("Fast_sims.R")


disease.state                 <-      initial.disease.state (demographic.ages  ,  v  , 1 ,  num.comps)
disease.state                 =       reduce.susceptibles (0, 5, disease.state, 0.85, num.comps, susceptible.indices)
disease.state                 =       reduce.susceptibles (6, 20, disease.state, 0.86, num.comps, susceptible.indices)
disease.state                 =       reduce.susceptibles (21, max(demographic.ages[, 1]), disease.state, 0.9, num.comps, susceptible.indices)
pop.by.age = number.of.each.age(demographic.ages,disease.state, num.comps)

plot(disease.state[susceptible.indices]/pop.by.age)
a = Run.Sims(num.steps)

disease.state3 = unlist(a[1])
infections.per.year3 = unlist(a[2])
all.infections = unlist(a[3])

a = Run.Sims.fast(num.steps)

l=0
k=0
for (j in 1 :100){
  disease.state                 <-      initial.disease.state (demographic.ages  ,  v  , 0.4 ,  num.comps)
 # disease.state                 =       reduce.susceptibles (0, 5, disease.state, 0.85, num.comps, susceptible.indices)
 # disease.state                 =       reduce.susceptibles (6, 20, disease.state, 0.86, num.comps, susceptible.indices)
#  disease.state                 =       reduce.susceptibles (21, max(demographic.ages[, 1]), disease.state, 0.95, num.comps, susceptible.indices)
  a = Run.Sims(500)
  all.infections = unlist(a[3])
  l = sum(all.infections) + l
  
  disease.state                 <-      initial.disease.state (demographic.ages  ,  v  , 0.5 ,  num.comps)
  disease.state                 =       reduce.susceptibles (0, 5, disease.state, 0.85, num.comps, susceptible.indices)
  disease.state                 =       reduce.susceptibles (6, 20, disease.state, 0.86, num.comps, susceptible.indices)
  disease.state                 =       reduce.susceptibles (21, max(demographic.ages[, 1]), disease.state, 0.95, num.comps, susceptible.indices)
  a = Run.Sims.fast(5000)
  all.infections1 = unlist(a[3])
  k  = sum(all.infections1) + k
}
l
k
