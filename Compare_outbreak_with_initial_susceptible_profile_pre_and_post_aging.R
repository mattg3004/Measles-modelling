library("prob")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("graphics")
source("SEIR_func.R")
source("Initial_conditions.R")


source("SEIR_measles.R")

num.repeats = 1
num.time.steps.pre.aging = 100
aging.section.steps  =  500
num.time.steps.post.aging = 100
infections.by.age.pre.aging  =  matrix(0, length(demographic.ages[, 1]), num.repeats)
infections.by.age.post.aging  =  matrix(0, length(demographic.ages[, 1]), num.repeats)
mixing.matrix                 =       matrix(1,length(demographic.ages[,1]),length(demographic.ages[,1]))
beta                          =       calibrate.beta(mixing.matrix, disease.state, infectious.indices, max.age, time.step, infectious.period, R_0)
for( k in 1 : num.repeats){
  
  disease.state                 <-      initial.disease.state (demographic.ages  ,  v  , 1 ,  num.comps)
  disease.state                 =       reduce.susceptibles (0, 5, disease.state, 0.98, num.comps, susceptible.indices)
  disease.state                 =       reduce.susceptibles (6, 20, disease.state, 0.97, num.comps, susceptible.indices)
  disease.state                 =       reduce.susceptibles (21, max(demographic.ages[, 1]), disease.state, 0.9, num.comps, susceptible.indices)
  
  disease.state[3]  =  1
  
  for (j in 1 : num.time.steps.pre.aging){
    N                           =        sum(disease.state)
    births.average              =        birth.rate * N * time.step 
    births.total                =        rpois(1, births.average)
    births.vac                  =        rbinom(1, births.total, v)              # number of childrn born who go straight to immune class
    disease.state[1]            =        disease.state[1] + (births.total - births.vac)
    disease.state[num.comps]    =        disease.state[num.comps] + births.vac
    new.infected                =        0
    number.infectious           =        0
    estimate.inf                =        0
    
    #foi.ages                    =      foi.approach3(R_0  ,  infectious.days  ,  inf.comp  ,  mixing.matrix  ,  num.comps  ,  disease.state  ,  max.age  ,  t  ,  time.step  ,  beta_1)
    foi.ages                    =      foi.by.next.gen ( mixing.matrix, disease.state, infectious.indices, time.step , infectious.period, beta, demographic.ages, num.comps)
    updated.state               =       matrix( 0  ,  num.comps*length(demographic.ages[,1])  ,  1)
    
    for ( i in 1:length(demographic.ages[, 1]) ){
      
      age                         =      demographic.ages[i, 1]
      foi.age                     =      min(1,foi.ages[age+1, ])
      estimate.inf                =      foi.age*disease.state[((age*num.comps)+1)]  +  estimate.inf
      change.matrix.by.age        =      stochastic.transmission.matrix.exposed.included  (age , ceiling(disease.state * prob.survival) , foi.age  , demographic.ages , time.step ,   rho, mu, num.comps)
      new.infected                =      new.infected   +    change.matrix.by.age[9]
      infections.by.age.pre.aging[i, k]     =      infections.by.age.pre.aging[i, k] + change.matrix.by.age[9]
      number.infectious           =      number.infectious    +    change.matrix.by.age[10]
      
      if(age == max(demographic.ages[,1])){
        updated.state[seq((((i-1)*num.comps)+1), num.comps*(i))]        <-      change.matrix.by.age[1:num.comps]
      }
      
      else{
        updated.state[seq((((i-1)*num.comps)+1), num.comps*(i+1))]      <-      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))] + change.matrix.by.age[seq(1,2*num.comps)]
      }
      
    }
    disease.state                 =       updated.state
  }  
  # reset to 0 infecteds or exposed individuals as we want to investigate the impact ofaging the population without decreasing susceptibles via infection
  disease.state[infectious.indices]  =  0
  disease.state[infectious.indices - 1]  =  0
  
  for ( l in 1 : aging.section.steps){
    updated.state               =        matrix( 0  ,  num.comps*length(demographic.ages[,1])  ,  1)
    N                           =        sum(disease.state)
    births.average              =        birth.rate * N * time.step 
    births.total                =        rpois(1, births.average)
    births.vac                  =        rbinom(1, births.total, v)              # number of childrn born who go straight to immune class
    disease.state[1]            =        disease.state[1] + (births.total - births.vac)
    disease.state[num.comps]    =        disease.state[num.comps] + births.vac
    
    for ( i in 1:length(demographic.ages[, 1]) ){
      
      age                         =      demographic.ages[i, 1]
      change.matrix.by.age        =      stochastic.transmission.matrix.exposed.included  (age , ceiling(disease.state * prob.survival) , 0  , demographic.ages , time.step ,   rho, mu, num.comps)
      if(age == max(demographic.ages[,1])){
        updated.state[seq((((i-1)*num.comps)+1), num.comps*(i))]        <-      change.matrix.by.age[1:num.comps]
      }
      
      else{
        updated.state[seq((((i-1)*num.comps)+1), num.comps*(i+1))]      <-      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))] + change.matrix.by.age[seq(1,2*num.comps)]
      }
      
    }
    disease.state                 =       updated.state
  }  
  
  
  # after aging the population by a given amount, introduce a single infected 
  disease.state[3]            =  1
  
  for (j in 1 : num.time.steps.post.aging){
    N                           =        sum(disease.state)
    births.average              =        birth.rate * N * time.step 
    births.total                =        rpois(1, births.average)
    births.vac                  =        rbinom(1, births.total, v)              # number of childrn born who go straight to immune class
    disease.state[1]            =        disease.state[1] + (births.total - births.vac)
    disease.state[num.comps]    =        disease.state[num.comps] + births.vac
    
    new.infected                =        0
    number.infectious           =        0
    estimate.inf                =        0
    
    #foi.ages                    =      foi.approach3(R_0  ,  infectious.days  ,  inf.comp  ,  mixing.matrix  ,  num.comps  ,  disease.state  ,  max.age  ,  t  ,  time.step  ,  beta_1)
    foi.ages                    =      foi.by.next.gen ( mixing.matrix, disease.state, infectious.indices, time.step , infectious.period, beta, demographic.ages, num.comps)
    updated.state               =       matrix( 0  ,  num.comps*length(demographic.ages[,1])  ,  1)
    
    for ( i in 1:length(demographic.ages[,1]) ){
      
      age                         =      demographic.ages[i, 1]
      foi.age                     =      min(1,foi.ages[age+1, ])
      estimate.inf                =      foi.age*disease.state[((age*num.comps)+1)]  +  estimate.inf
      if(age == 0)
      {
        foi.by.time[j,1]      =      foi.age
      }
      change.matrix.by.age        =      stochastic.transmission.matrix.exposed.included  (age , ceiling(disease.state * prob.survival) , foi.age  , demographic.ages , time.step ,   rho, mu, num.comps)
      new.infected                =      new.infected   +    change.matrix.by.age[9]
      infections.by.age.post.aging[i, k]     =      infections.by.age.post.aging[i, k] + change.matrix.by.age[9]
      number.infectious           =      number.infectious    +    change.matrix.by.age[10]
      
      if(age == max(demographic.ages[,1])){
        updated.state[seq((((i-1)*num.comps)+1), num.comps*(i))]        <-      change.matrix.by.age[1:num.comps]
      }
      
      else{
        updated.state[seq((((i-1)*num.comps)+1), num.comps*(i+1))]      <-      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))] + change.matrix.by.age[seq(1,2*num.comps)]
      }
      
    }
    disease.state                 =       updated.state
  }  
  
}

par(mfrow=c(1,2))
plot(infections.by.age.pre.aging)
plot(infections.by.age.post.aging)



