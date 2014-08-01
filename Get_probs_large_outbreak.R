initial.prop = seq(0, 0.25,0.01)
outbreak.threshold = 100000
mixing.matrix                 <-      full.mixing.matrix(contacts, demographic.ages)
p1 = probability.of.outbreak(0, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, num.comps, inf.comp, 250)
p2 = probability.of.outbreak(10, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, num.comps, inf.comp, 250)
p3 = probability.of.outbreak(30, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, num.comps, inf.comp, 250)
p4 = probability.of.outbreak(90, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, num.comps, inf.comp, 250)

mixing.matrix  =  matrix(1,length(demographic.ages[,1]),length(demographic.ages[,1]))
p5 = probability.of.outbreak(2, demographic.ages, initial.prop, mixing.matrix, infectious.indices, time.step, infectious.period, outbreak.threshold, num.comps, inf.comp, 250)

setEPS()
postscript("Prob_large_outbreak3.eps",width = 7.6,height = 6.5)
plot(initial.prop, p1, xlab = "initial proportion susceptible", ylab="probability large outbreak", ylim = (0:1.05))
points(initial.prop, p2, col = rgb(0,0,1),pch = 15)
points(initial.prop, p3, col = rgb(1,0,0), pch=16)
points(initial.prop, p4, col = rgb(1,0,1), pch = 17)
points(initial.prop, p5, col = rgb(0,0.8,0.8), pch = 18)
legend("topleft",c("age = 0", "age = 10", "age = 30", "age = 90", "uniform"), pch = c(1,15,16,17,18),col = c(rgb(0,0,0), rgb(0,0,1), rgb(1,0,0), rgb(1,0,1), rgb(0,0.8,0.8)))
dev.off()



p1_t1  =  probability.of.outbreak(0, demographic.ages, .1, mixing.matrix, infectious.indices, 5, infectious.period, outbreak.threshold, num.comps, inf.comp, 100)
