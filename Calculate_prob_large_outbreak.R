initial.prop = seq(0, 0.25,0.01)
outbreak.threshold = 100000
mixing.matrix                 <-      full.mixing.matrix(contacts, demographic.ages)
p1 = probability.of.outbreak(0, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, 20)
p2 = probability.of.outbreak(10, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, 250)
p3 = probability.of.outbreak(30, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, 250)
p4 = probability.of.outbreak(90, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, 250)

mixing.matrix  =  matrix(1,length(demographic.ages[,1]),length(demographic.ages[,1]))
p5 = probability.of.outbreak(2, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, 250)

setEPS()
postscript("Prob_large_outbreak3.eps",width = 7.6,height = 6.5)
plot(initial.prop, p1, xlab = "initial proportion susceptible", ylab="probability large outbreak", ylim = (0:1.05))
points(initial.prop, p2, col = rgb(0,0,1),pch = 15)
points(initial.prop, p3, col = rgb(1,0,0), pch=16)
points(initial.prop, p4, col = rgb(1,0,1), pch = 17)
points(initial.prop, p5, col = rgb(0,0.8,0.8), pch = 18)
legend("topleft",c("age = 0", "age = 10", "age = 30", "age = 90", "uniform"), pch = c(1,15,16,17,18),col = c(rgb(0,0,0), rgb(0,0,1), rgb(1,0,0), rgb(1,0,1), rgb(0,0.8,0.8)))
dev.off()

probability.of.outbreak <- function(age.initial.infected, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, num.sims){
  beta = calibrate.beta(mixing.matrix, disease.state, infectious.indices, max.age, time.step, infectious.period, R_0)
  prob.outbreak   =   matrix(0, length(initial.prop), 1)
  for (p in 1 : length(initial.prop)){
    
    for ( k in 1 : num.sims){
      disease.state                 <-      initial.disease.state(demographic.ages  ,  0  , initial.prop[p] ,  num.comps)
      disease.state[((3 * age.initial.infected) + 2)]  =  1
      total.infecteds  =  0
      
      
        
        while((total.infecteds < outbreak.threshold) & (sum(disease.state[infectious.indices]) > 0)){
          N                           =       sum(disease.state)
          births.average              =       birth.rate * N * time.step
          births.total                =       rpois(1, births.average)
          births.vac                  =       rbinom(1, births.total, v)              # number of childrn born who go straight to immune class
          disease.state[1]            =       disease.state[1] + (births.total - births.vac)
          disease.state[num.comps]    =       disease.state[num.comps] + births.vac
          new.infected = 0
          estimate.inf = 0
          updated.state               =      matrix( 0  ,  num.comps*length(demographic.ages[,1])  ,  1)
          foi.ages                    =      foi.by.next.gen ( mixing.matrix, disease.state, infectious.indices, time.step , infectious.period, beta, demographic.ages)
          for ( i in 1:length(demographic.ages[,1]) ){
            
            age                         =      demographic.ages[i, 1]
            
            #foi.age                     =      force.of.infection.by.age(age  ,  mixing.matrix  ,  disease.state  ,  beta  ,  gamma  ,  time.step  , inf.comp  ,  num.comps)
            foi.age                     =      min(1,foi.ages[age+1, ])
            estimate.inf                =      foi.age*disease.state[((age*3)+1)]  +  estimate.inf
            #print(paste("age =",age,",foi =",foi.age))
            if(age == 0)
            {
              foi.by.time[j,1]      =      foi.age
            }
            
            change.matrix.by.age        =      stochastic.transmission.matrix2(age , ceiling(disease.state * prob.survival), v , foi.age  , demographic.ages , time.step ,   rho)
            new.infected                =      new.infected   +    change.matrix.by.age[7]
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
    