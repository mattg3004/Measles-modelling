
full.mixing.matrix  <- function(contacts,demographic.ages){
  polymod.ages                =    contacts[, 1]
  contacts                    <-   contacts[, 2 : 16]
  max.age                     <-   max(demographic.ages[, 1])
  contact.matrix              =    matrix(0, max.age+1, max.age+1)
  colnames(contact.matrix)    =    c(0 : max.age)
  rownames(contact.matrix)    =    c(0 : max.age)
  num.polymod.ages            <-   length(polymod.ages)
  row.start.position          =    1
  
  for(i in 1 : ( max(polymod.ages))){
    for (j in 1: (max(polymod.ages))){
      i1 = floor((i-1) / 5) + 1
      j1 = floor((j-1) / 5) + 1
      contact.matrix[ i, j]  =  contacts[ i1, j1]/5
    }
  }
  
  for(i in 1 : ( max(polymod.ages))){
    for ( j in (max(polymod.ages) + 1) : (max.age+1)){
      i1 = floor((i-1) / 5) + 1
      j1 = floor(max(polymod.ages) / 5) + 1
      contact.matrix[ i, j]  =  contacts[ i1, j1]/5
    }
  }
    
  for ( i in (max(polymod.ages) + 1) : (max.age+1)){
    for ( j in 1 :  ( max(polymod.ages))){
      i1 = floor(max(polymod.ages) / 5) + 1
      j1 = floor((j-1) / 5) + 1
      contact.matrix[ i, j]  =  contacts[ i1, j1]/ 23
    }
  }
  
  last.entry.col  =  length(contacts[1, ])
  last.entry.row  =  length(contacts[ , 1])
  for ( i in (max(polymod.ages) + 1) : (max.age+1)){
    for ( j in (max(polymod.ages) + 1) : (max.age+1)){
      contact.matrix[ i, j]  =  contacts[ last.entry.row, last.entry.col]/ 23
    }
  }
  
    
  
  return(contact.matrix)
}




average.contacts.by.age<- function(mixing.matrix,max.age){
  
  contacts.by.age       =     colSums(mixing.matrix[,seq(1,(max.age+1))])
  av.num.contacts       =     mean(contacts.by.age)
  return(av.num.contacts)
}



deterministic.transmission.matrix<-function(age , disease.state , vacc.immune , foi  , demographic.ages , time.step ,  sigma  ,  rho){
  
  age.disease.state         =    disease.state[((age*3)+1):((age+1)*3)]
  number.of.age             =    sum(age.disease.state)
  u                         =    time.step/365
  change.matrix             =    matrix(0,11,1)
  new.infecteds             =    0
  
  A                         =    matrix(0,11,1)
  A[1]                      =    round((1-u)*(1-foi) * age.disease.state[1])
  A[2]                      =    round((1-u)*(foi) * age.disease.state[1])   +   round((1-u)*(1-sigma) * age.disease.state[2])
  A[3]                      =    round((1-u)*(sigma) * age.disease.state[2])   +   round((1-u)*(1-rho) * age.disease.state[3])
  A[4]                      =    round((1-u)*(rho) * age.disease.state[3])   +   round((1-u) * age.disease.state[4])
  A[5]                      =    round((1-u)   *  age.disease.state[5])
  A[6]                      =    round((u)*(1-foi) * age.disease.state[1])
  A[7]                      =    round((u)*(foi) * age.disease.state[1])   +   round((u)*(1-sigma) * age.disease.state[2])
  A[8]                      =    round((u)*(sigma) * age.disease.state[2])   +   round((u)*(1-rho) * age.disease.state[3])
  A[9]                      =    round((u)*(rho) * age.disease.state[3])   +   round((u) * age.disease.state[4])
  A[10]                     =    round((u)   *  age.disease.state[5])
  
  A[11]                     =    round((1-u)*(foi) * age.disease.state[1])   +   round((u)*(foi) * age.disease.state[1])
  
  return(A)
}


number.of.each.age <- function(demographic.ages,disease.state, num.comps){
  
  number.of.age = matrix(0,length(demographic.ages[,1]),1)
  for (i in 1:length(demographic.ages[,1]))
  {
    age = demographic.ages[i,1]
    number.of.age[i]          =    max(1, sum(disease.state[((age * num.comps)+1):((age+1) * num.comps)]))
  }
  #1
  return(number.of.age)
}



stochastic.transmission.matrix.exposed.included <- function(age , disease.state , foi , demographic.ages , time.step , rho, mu , num.comps){
  
  age.disease.state         =    disease.state[((age*num.comps)+1):((age+1)*num.comps)]
  number.of.age             =    sum(age.disease.state)
  prob.age.change           =    time.step/365
  change.matrix             =    matrix(0,10,1)
  new.infecteds             =    0
  newly.exposed             =    0
  # Assume that no one 1 or older gets vaccinnated
  u          =    prob.age.change
  
  susceptibles     =     infecteds    =     recovered     =     vac    =   exposed   =    matrix(0,8,1)
  
  if (age.disease.state[1] > 0){
    susceptibles     =    rmultinom(1, age.disease.state[1] , c( (1-u)*(1-foi) , (1-u)*foi ,0 , 0, u*(1-foi) , u*foi , 0, 0))
    newly.exposed    =    susceptibles[2]    +    susceptibles[6]
  }
  
  
  if (age.disease.state[2] > 0){
    exposed          =    rmultinom(1, age.disease.state[2] , c(0 , (1-u) * (1- mu), (1-u ) * mu,  0  , 0,  u  * (1 - mu)  ,  u  *  mu, 0))
    new.infecteds     =    exposed[3]   +   exposed[7]
  }
  
  
  if (age.disease.state[3] > 0){
    infecteds          =    rmultinom(1, age.disease.state[3] , c(0 , 0 ,  (1-u) * (1- rho) , (1-u ) * rho , 0 , 0 , u  * (1 - rho) ,  u  *  rho))
    
  }
  
  if (age.disease.state[4] > 0){
    recovered        =    rmultinom(1, age.disease.state[4], c( 0 , 0 , 0, (1-u) ,  0 , 0 , 0 ,  u))
  }
  
  
  change.matrix[seq(1,2*num.comps)]      =     susceptibles  + exposed  +  infecteds +  recovered 
  change.matrix[9]   =   newly.exposed
  change.matrix[10]  =   new.infecteds
  #print(new.infecteds)
  return(change.matrix)
}



stochastic.transmission.matrix.exposed.included.no.aging <- function(age , disease.state , foi , demographic.ages , time.step , rho, mu , num.comps){
  
  age.disease.state         =    disease.state[((age*num.comps)+1):((age+1)*num.comps)]
  number.of.age             =    sum(age.disease.state)
  prob.age.change           =    time.step/365
  change.matrix             =    matrix(0,10,1)
  new.infecteds             =    0
  newly.exposed             =    0
  # Assume that no one 1 or older gets vaccinnated
  u          =    0
  
  susceptibles     =     infecteds    =     recovered     =     vac    =   exposed   =    matrix(0,8,1)
  
  if (age.disease.state[1] > 0){
    susceptibles     =    rmultinom(1, age.disease.state[1] , c( (1-u)*(1-foi) , (1-u)*foi ,0 , 0, u*(1-foi) , u*foi , 0, 0))
    newly.exposed    =    susceptibles[2]    +    susceptibles[6]
  }
  
  
  if (age.disease.state[2] > 0){
    exposed          =    rmultinom(1, age.disease.state[2] , c(0 , (1-u) * (1- mu), (1-u ) * mu,  0  , 0,  u  * (1 - mu)  ,  u  *  mu, 0))
    new.infecteds     =    exposed[3]   +   exposed[7]
  }
  
  
  if (age.disease.state[3] > 0){
    infecteds          =    rmultinom(1, age.disease.state[3] , c(0 , 0 ,  (1-u) * (1- rho) , (1-u ) * rho , 0 , 0 , u  * (1 - rho) ,  u  *  rho))
    
  }
  
  if (age.disease.state[4] > 0){
    recovered        =    rmultinom(1, age.disease.state[4], c( 0 , 0 , 0, (1-u) ,  0 , 0 , 0 ,  u))
  }
  
  
  change.matrix[seq(1,2*num.comps)]      =     susceptibles  + exposed  +  infecteds +  recovered 
  change.matrix[9]   =   newly.exposed
  change.matrix[10]  =   new.infecteds
  #print(new.infecteds)
  return(change.matrix)
}




stochastic.transmission.matrix.no.aging<-function(age , disease.state , vacc.immune , foi  , demographic.ages , time.step ,   rho){
  
  age.disease.state         =    disease.state[((age*3)+1):((age+1)*3)]
  number.of.age             =    sum(age.disease.state)
  prob.age.change           =    time.step/365
  change.matrix             =    matrix(0,9,1)
  new.infecteds             =    0
  num.spreaders             =    0
  all.exposed               =    0
  # Assume that no one 1 or older gets vaccinnated
  if (age > 0)
  {
    vacc.immune      =    0
  }
  u          =    0
  v          =    vacc.immune
  
  susceptibles     =     infecteds    =     recovered     =     vac    =   exposed   =    matrix(0,6,1)
  
  if (age.disease.state[1] > 0){
    susceptibles     =    rmultinom(1, age.disease.state[1] , c( (1-u)*(1-foi) , (1-u)*foi , 0 ,  u*(1-foi) , u*foi , 0))
    new.infecteds    =    susceptibles[2]    +    susceptibles[5]
  }
  
  if (age.disease.state[2] > 0){
    infecteds          =    rmultinom(1, age.disease.state[2] , c(0 , (1-u) * (1- rho), (1-u ) * rho,  0  ,  u  * (1 - rho)  ,  u  *  rho))
    
  }
  
  if (age.disease.state[3] > 0){
    #number.recover   =    rbinom(age.disease.state[4],1,min(1,365*time.step/15))
    #prop.recover     =    number.recover/age.disease.state[4]
    recovered        =    rmultinom(1,age.disease.state[3],c( 0 , 0 , (1-u) ,  0 , 0 , u))
    
  }
  
  
  change.matrix[seq(1,6)]      =    c( susceptibles  + infecteds +  recovered) 
  change.matrix[7] = new.infecteds
  #print(new.infecteds)
  return(change.matrix)
}



contacts.per.age.group <- function(mixing.matrix  ,  demographic.ages){
  mixing.times.population <- mixing.matrix
  contacts.per.age <- mixing.matrix
  for (i in 1:length(demographic.ages[,1])){
    k<-matrix(demographic.ages[i,2],length(demographic.ages[,1]),1)
    mixing.times.population[,i] <- k*mixing.matrix[,i]
  }
  for (i in 1:length(demographic.ages[,1])){
    contacts.per.age[i,] <- mixing.times.population[i,]/demographic.ages[i,2]
  }
  
  return(rowSums(contacts.per.age))
}




refresh.plots  <- function(j, infecteds.by.time, difference.from.estimate, foi.by.time, demographic.ages, prop.infected.by.age, prop.susceptible.age, t, pop.by.year, prop.sus.time, total.prop.susceptible ){
  par(mfrow=c(3,3))
  plot(1:j,infecteds.by.time[1:j],type="l",ylab = "new inf",col=rgb(runif(1),runif(1),runif(1)))
  plot(1:j,difference.from.estimate[1:j],type="l",ylab = "diff.est",col=rgb(runif(1),runif(1),runif(1)))
  plot(1:j,cumsum(infecteds.by.time[1:j]),type="l",ylab = "cum inf",col=rgb(runif(1),runif(1),runif(1)))
  plot(1:j,total.prop.susceptible[1:j],type="l",ylab = "total prop sus",col=rgb(runif(1),runif(1),runif(1)))
  
  plot(demographic.ages[,1], prop.infected.by.age,type="l",ylab = "inf prop",xlab = "age",col=rgb(runif(1),runif(1),runif(1)))
  plot(1:j,prop.sus.time[1:j], type="b",ylab = "prop contacts sus",xlab = "step",col=rgb(runif(1),runif(1),runif(1)))
  plot(demographic.ages[,1],prop.susceptible.age,type="b",ylab = "prop susceptible",xlab = "age",col=rgb(runif(1),runif(1),runif(1)))
  years   =   ceiling(t/365)
  if (years > 0)
  {
    infections.per.year        =    matrix(0,years,1)
    for (pp in 1:years)
    {
      infections.per.year[pp]   =   sum(infecteds.by.time[((floor(pp-1)*(365/time.step)) +1):floor(pp*365/time.step)])
    }
    plot(1:years, infections.per.year, type="b", ylab = "num infs", xlab = "year", col= rgb(runif(1), runif(1), runif(1)))
    plot(1:years, 100*infections.per.year/pop.by.year[1:years], type="b", ylab = "% of pop infected", xlab = "year", col= rgb(runif(1), runif(1), runif(1)))
  }
  
}



initial.disease.state <- function(demographic.ages  ,  vacc.immune  ,   initial.prop.susceptible  ,  num.comps){
  disease.state                 =      matrix(0,num.comps*length(demographic.ages[,1]),1)
  
  for (i in 1:length(demographic.ages[,1])){
    disease.state[(((i-1)*num.comps)+1):((i)*num.comps)]      =    ceiling(c( demographic.ages[i,2]*initial.prop.susceptible , 0, 0, demographic.ages[i,2]*(1-initial.prop.susceptible)))
  }
  
  return(disease.state)
}



supplementary.vacc <- function(disease.state, supp.vac , max.age, num.comps){
  disease.state
  for (age in 0:max.age)
  {
    age.disease.state                 =    disease.state[((age*num.comps)+1):((age+1)*num.comps)]
    p = disease.state[((age*num.comps)+1)]
    l                                 =    round(disease.state[((age*num.comps)+1)]*(1-supp.vac))
    disease.state[((age*num.comps)+1)]        =    l
    disease.state[(age+1)*num.comps]          =    round(disease.state[((age*num.comps)+1)]*(supp.vac)) + disease.state[((age+1)*num.comps)]
  }
  return(disease.state)
}




foi.by.next.gen <- function ( mixing.matrix, disease.state, infectious.indices, time.step , infectious.period, beta, demographic.ages, num.comps){
  infectious.by.age  =  disease.state[infectious.indices]
  pop.by.age  =  number.of.each.age (demographic.ages,disease.state, num.comps)
  foi.by.age  =  (colSums((infectious.by.age * t( mixing.matrix) * beta) ) / pop.by.age) * min( 1, time.step / infectious.period)
  return(foi.by.age)  
}



calibrate.beta <- function (mixing.matrix, disease.state, infectious.indices, max.age, time.step, infectious.period, R_0, population.by.age){
  num.infectious  =  matrix(0, (max.age + 1), 1)
  average.infectious  =  matrix(0, (max.age + 1), 1)
  for ( i in 1 : (max.age + 1) ){
    num.infectious[i]  =  1
    average.infectious[i]  =  sum(mixing.matrix[, i]) * min(1, time.step / infectious.period)
    num.infectious[i]  =  0
  }
  g = 0
  for (i in 1 : ( max.age + 1)){
    g = g + (population.by.age[i] / sum(population.by.age[1 : (max.age + 1)])) * average.infectious[i]
  }
  beta  = R_0 * min(1, time.step / infectious.period) / g
  return(beta)
}



reduce.susceptibles <- function(min.age, max.age, disease.state, proportion.sus.to.remove, num.comps, susceptible.indices){
  
  susceptibles  =  susceptible.indices[(min.age + 1) : (max.age + 1)]
  disease.state [ (susceptibles + (num.comps - 1) )]  = disease.state [ (susceptibles + (num.comps - 1) )]   +   round( disease.state[ susceptibles  ] * ( proportion.sus.to.remove ) )
  disease.state [ susceptibles ]  =  round(disease.state[ susceptibles ] * ( 1 - proportion.sus.to.remove ) )

  return(disease.state)
}



draw.sus <- function(x){
  u = time.step / 365
  foi = min(1, x[1])
  numbers = x[2]
  rmultinom(1, numbers , c( (1-u)*(1-foi) , (1-u)*foi ,0 , 0, u*(1-foi) , u*foi , 0, 0))
}

draw.exposed <- function(x){
  u = time.step / 365
  numbers = x[1]
  mu = min(1, time.step / incubation.period)
  rmultinom(1, numbers, c(0 , (1-u) * (1- mu), (1-u ) * mu,  0  , 0,  u  * (1 - mu)  ,  u  *  mu, 0))
}


draw.infecteds <- function(x){
  u = time.step / 365
  rho = min(1, time.step / infectious.period)
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 ,  (1-u) * (1- rho) , (1-u ) * rho , 0 , 0 , u  * (1 - rho) ,  u  *  rho))
}

draw.recovered <- function(x){
  u = time.step / 365
  numbers = x[1]
  rmultinom(1, numbers , c(0 , 0 ,  0, (1-u ) , 0 , 0 , 0 ,  u  ))
}
