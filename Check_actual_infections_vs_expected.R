disease.state                 <-      initial.disease.state(demographic.ages  ,  v  , initial.prop.susceptible ,  num.comps)
age = 10
time.step = 8
disease.state[(age*num.comps) + 3] = 1
new.infected                =        0
estimate.inf                =        0
for (i in 1:2000){
  disease.state                 <-      initial.disease.state(demographic.ages  ,  v  , initial.prop.susceptible ,  num.comps)
  age = 0
  disease.state[(age*num.comps) + 3] = 1

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
}
new.infected
estimate.inf
