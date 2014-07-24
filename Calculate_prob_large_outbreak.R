
probability.of.outbreak <- function(age.initial.infected, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, num.comps, inf.comp, num.sims){
  beta = calibrate.beta(mixing.matrix, disease.state, infectious.indices, max.age, time.step, infectious.period, R_0)
  prob.outbreak   =   matrix(0, length(initial.prop), 1)
  for (p in 1 : length(initial.prop)){
    
    for ( k in 1 : num.sims){
      disease.state                 <-      initial.disease.state(demographic.ages  ,  0  , initial.prop[p] ,  num.comps)
      disease.state[((num.comps * age.initial.infected) + inf.comp)]  =  1
      total.infecteds  =  0
      
      
        
        while((total.infecteds < outbreak.threshold) & (sum(disease.state[infectious.indices]) + sum(disease.state[infectious.indices - 1]) > 0)){
          N                           =       sum(disease.state)
          births.average              =       birth.rate * N * time.step
          births.total                =       rpois(1, births.average)
          births.vac                  =       rbinom(1, births.total, v)              # number of childrn born who go straight to immune class
          disease.state[1]            =       disease.state[1] + (births.total - births.vac)
          disease.state[num.comps]    =       disease.state[num.comps] + births.vac
          new.infected                =        0
          number.infectious           =        0
          estimate.inf                =        0
          
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
            number.infectious           =      number.infectious    +    change.matrix.by.age[10]
            
            if(age == max(demographic.ages[,1])){
              updated.state[seq((((i-1)*num.comps)+1), num.comps*(i))]        <-      change.matrix.by.age[1:num.comps]
            }
            
            else{
              updated.state[seq((((i-1)*num.comps)+1), num.comps*(i+1))]      <-      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))] + change.matrix.by.age[seq(1,2*num.comps)]
            }
            
          }
          
          total.infecteds    =    new.infected  +  total.infecteds
          disease.state   =   updated.state
          
        }
      
      
      if ( total.infecteds > outbreak.threshold){
        prob.outbreak[p]  =  prob.outbreak[p]  + 1
      }
      
    }
    

  }
  
  return(prob.outbreak/num.sims)
}
    