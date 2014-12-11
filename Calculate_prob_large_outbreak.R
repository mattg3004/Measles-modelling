
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
    