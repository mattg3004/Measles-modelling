library("prob")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("graphics")
source("SEIR_func.R")
source("Initial_conditions.R")

Run.Sims.fast <-function(num.steps){
  t = 0
  total.births.vac = 0
  foi.by.time = matrix(0, num.steps, 5)
  initial.vac_0  =  disease.state[4]
  for(j in 1:num.steps){
    #print(paste('t =', t, ', beta =', beta))
    foi.by.time[j, 5]           =       sum(disease.state)
    N                           =       sum(disease.state)
    births.average              =       birth.rate * N * time.step
    births.total                =       rpois(1, births.average)
    births.vac                  =       rbinom(1, births.total, v)              # number of childrn born who go straight to immune class
    total.births.vac = total.births.vac + births.vac
    #print(paste('max births.vac should be =', total.births.vac + initial.vac_0))
    disease.state[1]            =       disease.state[1] + (births.total - births.vac)
    disease.state[num.comps]    =       disease.state[num.comps] + births.vac
    migrant.infecteds           =       rpois(  max.age + 1  ,  av.migrants.per.age.per.day * time.step)
    total.migrants              =       total.migrants   +  sum(migrant.infecteds) 
    #print(paste("total.migrants =", sum(total.migrants)))
    
    beta                          =       beta_0 * (1 + beta_1 * cos(2 * pi * t / 365))
    disease.state[migrant.indices]    =     disease.state[migrant.indices] + migrant.infecteds
    
    foi.by.time[j, 2]           =       sum(disease.state[infectious.indices])
    # print(paste("infecteds =",sum(disease.state[seq(inf.comp,length(updated.state),num.comps)])))
    new.infected                =        0
    number.infectious           =        0
    estimate.inf                =        0
    
    #foi.ages                    =      foi.approach3(R_0  ,  infectious.days  ,  inf.comp  ,  mixing.matrix  ,  num.comps  ,  disease.state  ,  max.age  ,  t  ,  time.step  ,  beta_1)
    foi.ages                    =      foi.by.next.gen ( mixing.matrix, disease.state, infectious.indices, time.step , infectious.period, beta, demographic.ages, num.comps)
    
    foi.by.time[j,1]      =      foi.ages[1]
    
    estimate.inf                =      sum(foi.ages*disease.state[susceptible.indices])
    updated.state               =      matrix( 0  ,  length(disease.state)  ,  1)
    disease.state               =      ceiling(disease.state * prob.survival)
  
    x              =      matrix(0, length(demographic.ages[,1]), 2)
    x[ , 1]        =      foi.ages
    x[ , 2]        =      disease.state[susceptible.indices]
    sus.outs       =      apply(x, 1, draw.sus)
    
    x1             =      matrix(0, length(demographic.ages[,1]), 2)
    x1[ , 1]       =      disease.state[exposed.indices]
    exposed.out    =      apply(x1, 1, draw.exposed)
    
    
    x2             =      matrix(0, length(demographic.ages[,1]), 2)
    x2[ , 1]       =      disease.state[infectious.indices]
    inf.out        =      apply(x2, 1, draw.infecteds)
    
    x3             =      matrix(0, length(demographic.ages[,1]), 2)
    x3[ , 1]       =      disease.state[recovered.indices]
    recovered.out   =     apply(x3, 1, draw.recovered)
    
    new.infected       =   sum(sus.outs[2, ])  +  sum(sus.outs[6, ])
    number.infectious  =   sum(exposed.out[3, ]) + sum(exposed.out[7, ])
    
    for (p in 1 : (length(demographic.ages[ ,1]) - 1)){
      updated.state[seq(((num.comps*( p-1 )) + 1) ,num.comps*(p+1))]   =   updated.state[seq(((num.comps*( p-1 )) + 1) , num.comps*(p+1))]  +  sus.outs[ , p]  +  exposed.out[ , p]  +  inf.out[ , p] + recovered.out[ , p]
    }
    p = length(demographic.ages[ ,1])
    updated.state[seq(((num.comps*( p-1 )) + 1) ,num.comps*(p))]   =   updated.state[seq(((num.comps*( p-1 )) + 1) , num.comps*(p))]  +  sus.outs[1:4 , p]  +  exposed.out[1:4 , p]  +  inf.out[1:4 , p]  +  recovered.out[ 1:4, p]
    
    disease.state                 =       updated.state
    #print(paste('births.vac is =',  disease.state[4]))
    number.by.age                 =       number.of.each.age(demographic.ages,updated.state, num.comps)
    prop.infected.by.age          =       updated.state[infectious.indices]/number.by.age
    prop.susceptible              =       updated.state[susceptible.indices]/number.by.age
    total.infecteds               =       sum(updated.state[infectious.indices])
    total.prop.susceptible[j]     =       sum(updated.state[susceptible.indices])/sum(updated.state)
    difference.from.estimate[j]   =       new.infected - estimate.inf
    foi.by.time[j,4]              =       new.infected
    foi.by.time[j,3]              =       estimate.inf
    prop.contacts.sus.by.time[j]  =       proportion.of.susceptible.contacts(mixing.matrix, demographic.ages, disease.state, susceptible.indices)
    average.infection.age[j]      =       sum(updated.state[infectious.indices]*demographic.ages[,1])
    infecteds.by.time[j]          =       sum(updated.state[infectious.indices])
    infecteds.by.time[j]          =       new.infected
    susceptibles.by.time[j]       =       sum(updated.state[susceptible.indices])
    t                    =       t   +    time.step
    
    if ( t > days.per.extra.vaccination * (number.sup.vacs +1)){
      number.sup.vacs              =       number.sup.vacs   +   1
      disease.state                =       supplementary.vacc(disease.state  ,  supp.vac  ,  max.age.for.supp.vac, num.comps)
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
    if(j %% 399  == 0){
      disease.state.fast.pre.infection = disease.state
    }
  }
  years   =   ceiling(t/365)
  infections.per.year        =    matrix(0,years,1)
  for (pp in 1:years)
  {
    infections.per.year[pp]   =   sum(infecteds.by.time[((floor(pp-1)*(365/time.step)) +1):floor(pp*365/time.step)])
  }
  return(list(disease.state, infections.per.year, infecteds.by.time))
}