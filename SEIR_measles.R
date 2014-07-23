  library("prob")
  library("gplots")
  library("ggplot2")
  library("RColorBrewer")
  library("graphics")
  source("SIR_func.R")
  source("Funcs_foi_approach2.R")
  source("Initial_conditions.R")
  
  Run.Sims <-function(num.steps){
  
  for(j in 1:num.steps){
    #print(paste("j =",j))
    N                           =       sum(disease.state)
    births.average              =       birth.rate*N*time.step
    births.total                =       rpois(1, births.average)
    births.vac                  =       rbinom(1, births.total, v)              # number of childrn born who go straight to immune class
    disease.state[1]            =       disease.state[1] + (births.total - births.vac)
    disease.state[num.comps]    =       disease.state[num.comps] + births.vac
    migrant.infecteds           =       rpois(  max.age + 1  ,  av.migrants.per.age.per.day * time.step)
    #print(paste("infected.migrants =", sum(migrant.infecteds)))
    updated.state               =       matrix( 0  ,  num.comps*length(demographic.ages[,1])  ,  1)
    
    disease.state[migrant.indices]    =     disease.state[migrant.indices] + migrant.infecteds
    foi.by.time[j, 4]           =       sum(disease.state[infectious.indices])
    # print(paste("infecteds =",sum(disease.state[seq(inf.comp,length(updated.state),num.comps)])))
    pre.infected                =       sum(disease.state[seq( inf.comp  ,  length(disease.state)  ,  num.comps )])
    
    new.infected                =        0
    estimate.inf                =        0

    #foi.ages                    =      foi.approach3(R_0  ,  infectious.days  ,  inf.comp  ,  mixing.matrix  ,  num.comps  ,  disease.state  ,  max.age  ,  t  ,  time.step  ,  beta_1)
    foi.ages                    =      foi.by.next.gen ( mixing.matrix, disease.state, infectious.indices, time.step , infectious.period, beta, demographic.ages, num.comps)
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
      if(age == max(demographic.ages[,1])){
        updated.state[seq((((i-1)*num.comps)+1), num.comps*(i))]        <-      change.matrix.by.age[1:num.comps]
      }
      
      else{
        updated.state[seq((((i-1)*num.comps)+1), num.comps*(i+1))]      <-      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))] + change.matrix.by.age[seq(1,2*num.comps)]
      }
      
    }
    number.by.age                 =       number.of.each.age(demographic.ages,updated.state)
    prop.infected.by.age          =       updated.state[infectious.indices]/number.by.age
    prop.susceptible              =       updated.state[susceptible.indices]/number.by.age
    total.infecteds               =       sum(updated.state[infectious.indices])
    total.prop.susceptible[j]     =      sum(updated.state[susceptible.indices])/sum(updated.state)
    difference.from.estimate[j]   =      new.infected - estimate.inf
    foi.by.time[j,2]              =      new.infected
    foi.by.time[j,3]              =      sum(migrant.infecteds)
    foi.by.time[j,5]              =      estimate.inf
    prop.contacts.sus.by.time[j]  =      proportion.of.susceptible.contacts(mixing.matrix, demographic.ages, disease.state, susceptible.indices)
    average.infection.age[j]      =      sum(updated.state[infectious.indices]*demographic.ages[,1])
    infecteds.by.time[j]          =      sum(updated.state[infectious.indices])
    infecteds.by.time[j]          =      new.infected
    susceptibles.by.time[j]       =       sum(updated.state[susceptible.indices])
    disease.state                 =       updated.state
    t                    =       t   +    time.step
    
    if ( t > days.per.extra.vaccination * (number.sup.vacs +1)){
      number.sup.vacs              =       number.sup.vacs   +   1
      disease.state                =       supplementary.vacc(disease.state  ,  supp.vac  ,  max.age.for.supp.vac)
    }
  
    years1   =   ceiling(t/365)
    if (years1 > years2)
    {
      years2                        =          years1
      pop.by.year[years2]           =          sum(disease.state)
    }
    if (j %% 500 == 0)
    {
      print(paste("j =",j))
      refresh.plots(j, infecteds.by.time, difference.from.estimate, foi.by.time, demographic.ages, prop.infected.by.age, prop.susceptible, t, pop.by.year, prop.contacts.sus.by.time, total.prop.susceptible)
      
    }
  }
  years   =   ceiling(t/365)
  infections.per.year        =    matrix(0,years,1)
  for (pp in 1:years)
  {
     infections.per.year[pp]   =   sum(infecteds.by.time[((floor(pp-1)*(365/time.step)) +1):floor(pp*365/time.step)])
  }
  return(list(disease.state, infections.per.year))
  }