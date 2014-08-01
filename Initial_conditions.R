  demographic.ages               =       read.csv("Ages.csv")
  demographic.ages               =       round( read.csv("Pop_Zimbabwe.csv") )
  
  #starting.proportion           =       1
  #demographic.ages              =       demographic.ages1
  #demographic.ages[,2]          =       ceiling(demographic.ages1[,2]/starting.proportion)
  contacts                      <-      read.csv("Contacts.csv")
  #rescale.pop.to.initial.level  =       1 
  #rescaled                      =       1
  initial.pop                   =       sum(demographic.ages[,2])
  updated.state                 =       matrix(0, num.comps * length(demographic.ages[,1]), 1)
  
  num.comps                     =       4          # number of compartments in the model
  inf.comp                      =       3          # what compartment keeps track of the number of infectious for each age group
  t                             =       0
  R_0                           =       15
  #time.step                     =       5          # number of days to step forward
  vacc.prop                     =       1
  vacc.success                  =       0.85
  v                             =       vacc.prop*vacc.success

  incubation.period             =       8         # length of incubation period on average
  mu                            =       min(1, time.step/incubation.period)   # probability of moving from exposed to infectious class during a timestep
  infectious.period             =       5           # number of days spent in the infected class on average
  rho                           =       min(1,time.step/infectious.period)    # probability of losing infectiousness. not necessarily recovered from the infection, but no longer infectious.
  max.age                       =       10         # assume that when R_0 was calculated previously, this was the maximum age of the people spreading measles
  infectious.indices            =       seq(inf.comp,length(updated.state), num.comps)
  migrant.indices               =       seq( inf.comp  ,  num.comps*(max.age+1) ,  num.comps )
  susceptible.indices           =       seq( 1  ,  length(updated.state) ,  num.comps )
  exposed.indices               =       seq(2, length(updated.state), num.comps)
  recovered.indices             =       seq(4, length(updated.state), num.comps)
  
  initial.prop.susceptible      =       1
  days.per.extra.vaccination    =       50 * 365     # do an additonal vaccination campaign every this many days
  number.sup.vacs               =       0           # number of supplementary vaccination campaigns so far
  max.age.for.supp.vac          =       3           # max age for supplementary vaccinations
  supp.vac                      =       0.5         # proportion of people who are up to the age of the max age for supplementary vaccination who will move to vaccinated class
  
  mixing.matrix                 <-      full.mixing.matrix(contacts, demographic.ages)      # average number of people met of each age grooup, stratified by age
  mixing.matrix                 =       matrix(1,length(demographic.ages[,1]),length(demographic.ages[,1]))
  beta_0                        =       calibrate.beta(mixing.matrix, disease.state, infectious.indices, max.age, time.step, infectious.period, R_0)
  beta_1                        =       0.8
  beta                          =       beta_0 * (1 + beta_1 * cos(2 * pi * t / 365))
  disease.state                 <-      initial.disease.state(demographic.ages  ,  v  , initial.prop.susceptible ,  num.comps)
  av.migrants.per.age.per.day   =       1/(365 * (max.age + 1))
  
  birth.rate                    =       30/(1000*365)
  death.rate                    =       15/(1000*365)
  prob.survival                 =       1-(death.rate)*(time.step)

  
  infecteds.by.time             =       matrix(0, num.steps, 1)
  susceptibles.by.time          =       matrix(0, num.steps, 1)
  average.infection.age         =       matrix(0, num.steps, 1)
  foi.by.time                   =       matrix(0, num.steps, 6)
  all.infectious.by.time        =       matrix(0, num.steps, 1)
  all.exposed.by.time           =       matrix(0, num.steps, 1)
  difference.from.estimate      =       matrix(0, num.steps, 1)
  pop.each.step                 =       matrix(0, num.steps, 1)
  prop.contacts.sus.by.time     =       matrix(0, num.steps, 1)
  total.prop.susceptible        =       matrix(0, num.steps, 1)
  pop.by.year                   =       matrix(0, ((num.steps*time.step/365)+1), 1)
  years2                        =       0