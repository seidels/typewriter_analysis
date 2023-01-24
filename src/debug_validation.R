\## ---------------------------
##
## Script name: debug_validation
##
## Purpose of script: Debug/plot validation outputs
##
## Author: Antoine Zwaans
##
## Date Created: 2023-01-11
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac 

setwd("/Users/azwaans/typewriter_analysis/")   # Antoine's working directory 

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(tidyr)
require(data.table)
require(EvoPhylo)
require(ggplot2)
require(cowplot)

logLikelihoodsClock01 <- exp(c(-5.88190164248573, -5.213172146682981, -4.857733878998448, -4.648883453048166, -4.539651507488826, -4.51797788250131, -4.594137864203139, -4.813186504339441, -5.339629672070573))
logLikelihoodsClock1 <- exp(c(-13.742782681155251, -12.697414907005953, -12.144068732084127, -11.821946169652055, -11.64954781734308, -11.598802618337844, -11.671709389145933, -11.908957546641684, -12.473022136291313))
logLikelihoodsClock10 <- exp(c(-140.4608056222875, -139.30351993193423, -138.7478712787002, -138.4643350480603, -138.35961623384998, -138.40733956701092, -138.62010605261264, -139.06332796411823, -139.96138077347288))




insertProbability <- seq(0.9,0.1,by=-0.1)

data <- data.frame(probability=insertProbability,clock01=logLikelihoodsClock01,clock1=logLikelihoodsClock1,clock10=logLikelihoodsClock10)
#data <- pivot_longer(data,!probability,names_to = "clockRate",values_to = "likelihood")

ggplot(data=data,aes(x=probability,y=likelihood,group=clockRate,color=clockRate)) + geom_line() + xlab("Insert 1 probability") + ylab("Likelihood (of unedited barcode at origin)")
clock01 <- ggplot(data=data,aes(x=insertProbability,y=clock01)) + geom_line() + xlab("Insert 1 probability") + ylab("Likelihood (clockRate = 0.1)")
clock1 <- ggplot(data=data,aes(x=insertProbability,y=clock1)) + geom_line() + xlab("Insert 1 probability") + ylab("Likelihood (clockRate = 1)")
clock10 <- ggplot(data=data,aes(x=insertProbability,y=clock10)) + geom_line() + xlab("Insert 1 probability") + ylab("Likelihood (clockRate = 10)")


plot_grid(clock01,clock1,clock10,nrow=1)

logLikelihoodsClock01 <- c(4.65597417732607E-5, 1.4715177646857703E-4, 2.5349192743219707E-4, 3.310914970542982E-4, 3.5925726676898676E-4, 3.310914970542982E-4, 2.5349192743219707E-4, 1.4715177646857703E-4, 4.65597417732607E-5)
logLikelihoodsClock1 <- c(5.745928610564492E-5, 1.8159971904993955E-4, 3.128338910196222E-4, 4.085993678623637E-4, 4.4335868908676625E-4, 4.085993678623637E-4, 3.128338910196222E-4, 1.8159971904993955E-4, 5.745928610564492E-5)
logLikelihoodsClock10 <- c(4.708221157151317E-40, 1.4880303904083176E-39, 2.563364852226827E-39, 3.3480683784187123E-39, 3.632886695332805E-39, 3.3480683784187123E-39, 2.563364852226827E-39, 1.4880303904083176E-39, 4.708221157151317E-40)




insertProbability <- seq(0.9,0.1,by=-0.1)

data <- data.frame(probability=insertProbability,clock01=logLikelihoodsClock01,clock1=logLikelihoodsClock1,clock10=logLikelihoodsClock10)
#data <- pivot_longer(data,!probability,names_to = "clockRate",values_to = "likelihood")

ggplot(data=data,aes(x=probability,y=likelihood,group=clockRate,color=clockRate)) + geom_line() + xlab("Insert 1 probability") 
clock01_0000 <- ggplot(data=data,aes(x=insertProbability,y=clock01)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 0000 root (Rate = 0.1)")
clock1_0000 <- ggplot(data=data,aes(x=insertProbability,y=clock1)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 0000 root (Rate = 1)")
clock10_0000 <- ggplot(data=data,aes(x=insertProbability,y=clock10)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 0000 root (Rate = 10)")


plot_grid(clock01,clock1,clock10,nrow=1)


logLikelihoodsClock01 <- c(9.19698602928606E-4, 0.003678794411714424, 0.008277287426357452, 0.014715177646857695, 0.022992465073215146, 0.03310914970542981, 0.04506523154350168, 0.05886071058743078, 0.07449558683721708)
logLikelihoodsClock1 <- c(1.1349982440621215E-5, 4.539992976248486E-5, 1.0214984196559093E-4, 1.8159971904993944E-4, 2.8374956101553033E-4, 4.085993678623637E-4, 5.561491395904394E-4, 7.263988761997578E-4, 9.193485776903183E-4)
logLikelihoodsClock10 <- c(9.300189940052121E-43, 3.7200759760208484E-42, 8.37017094604691E-42, 1.4880303904083394E-41, 2.3250474850130305E-41, 3.348068378418764E-41, 4.5570930706255385E-41, 5.952121561633358E-41, 7.533153851442219E-41)




insertProbability <- seq(0.9,0.1,by=-0.1)

data <- data.frame(probability=insertProbability,clock01=logLikelihoodsClock01,clock1=logLikelihoodsClock1,clock10=logLikelihoodsClock10)
#data <- pivot_longer(data,!probability,names_to = "clockRate",values_to = "likelihood")

ggplot(data=data,aes(x=probability,y=likelihood,group=clockRate,color=clockRate)) + geom_line() + xlab("Insert 1 probability") + ylab("Likelihood (of unedited barcode at origin)")
clock01_1000 <- ggplot(data=data,aes(x=insertProbability,y=clock01)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 1000 root (Rate = 0.1)") + theme(axis.text.y = element_text(size = 6)) 
clock1_1000 <- ggplot(data=data,aes(x=insertProbability,y=clock1)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 1000 root (Rate = 1)") + theme(axis.text.x = element_text(size = 6)) 
clock10_1000 <- ggplot(data=data,aes(x=insertProbability,y=clock10)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 1000 root Rate = 10)") + theme(axis.text.x = element_text(size = 6)) 


plot_grid(clock01,clock1,clock10,clock01_0000,clock1_0000,clock10_0000,clock01_1000,clock1_1000,clock10_1000,nrow=3,ncol=3)

logLikelihoodsClock01 <- c(0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233)
logLikelihoodsClock1 <- c(4.539992976248485E-5, 4.539992976248485E-5, 4.539992976248485E-5, 4.539992976248485E-5, 4.539992976248485E-5, 4.539992976248485E-5, 4.539992976248485E-5, 4.539992976248485E-5, 4.539992976248485E-5)
logLikelihoodsClock10 <- c(3.720075976020836E-44, 3.720075976020836E-44, 3.720075976020836E-44, 3.720075976020836E-44, 3.720075976020836E-44, 3.720075976020836E-44, 3.720075976020836E-44, 3.720075976020836E-44, 3.720075976020836E-44)




insertProbability <- seq(0.9,0.1,by=-0.1)

data <- data.frame(probability=insertProbability,clock01=logLikelihoodsClock01,clock1=logLikelihoodsClock1,clock10=logLikelihoodsClock10)
#data <- pivot_longer(data,!probability,names_to = "clockRate",values_to = "likelihood")

ggplot(data=data,aes(x=probability,y=likelihood,group=clockRate,color=clockRate)) + geom_line() + xlab("Insert 1 probability") + ylab("Likelihood (of unedited barcode at origin)")
clock01_1200 <- ggplot(data=data,aes(x=insertProbability,y=clock01)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 1200 root (Rate = 0.1)") + theme(axis.text.y = element_text(size = 6)) 
clock1_1200 <- ggplot(data=data,aes(x=insertProbability,y=clock1)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 1200 root (Rate = 1)") + theme(axis.text.x = element_text(size = 6)) 
clock10_1200 <- ggplot(data=data,aes(x=insertProbability,y=clock10)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 1200 root Rate = 10)") + theme(axis.text.x = element_text(size = 6)) 

plot_grid(clock01,clock1,clock10,clock01_0000,clock1_0000,clock10_0000,clock01_1000,clock1_1000,clock10_1000,clock01_1200,clock1_1200,clock10_1200,nrow=4,ncol=3)



#with sequences 2200 2200
likelihoodTree <- exp(c(-8.088165829399674, -6.611837799605418, -5.716406269603338, -5.0613661338149765, -4.539651507488826, -4.103352404994482, -3.726767645326548, -3.3944568674663973, -3.096390587666701)) 
likelihood2200 <- c(0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233)
likelihood2000 <- c(9.19698602928606E-4, 0.003678794411714424, 0.008277287426357452, 0.014715177646857695, 0.022992465073215146, 0.03310914970542981, 0.04506523154350168, 0.05886071058743078, 0.07449558683721708)
likelihood0000 <- c(5.74811626830379E-7, 9.196986029286064E-6, 4.655974177326069E-5, 1.4715177646857703E-4, 3.5925726676898676E-4, 7.44955868372171E-4, 0.001380122716019739, 0.0023544284234972325, 0.0037713390836341155)





insertProbability <- seq(0.9,0.1,by=-0.1)

data <- data.frame(probability=insertProbability,likelihood=likelihoodTree,internal0=likelihood0000,internal2=likelihood2000,internal22=likelihood2200)
#data <- pivot_longer(data,!probability,names_to = "clockRate",values_to = "likelihood")

likelihood <- ggplot(data=data,aes(x=insertProbability,y=likelihood)) + geom_line() + xlab("Insert 1 probability") + ylab("Likelihood (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")
state0 <- ggplot(data=data,aes(x=insertProbability,y=internal0)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 0000 (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")
state2 <- ggplot(data=data,aes(x=insertProbability,y=internal2)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 2000 (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")
state22 <- ggplot(data=data,aes(x=insertProbability,y=internal22)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 2200 (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")


plot_grid(likelihood,state0,state2,state22,nrow=1)



#with sequences 21100 2100
likelihoodTree <- exp(c(-6.411748421698742, -5.9870908979224895, -5.8833462549899, -5.940503564432857, -6.1205879097907125, -6.425375450145432, -6.889867436932458, -7.613255604161904, -8.93076419355177) )
likelihood2100 <- c(0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233)
likelihood2000 <- c(0.016761507038373844, 0.011772142117486159, 0.007886415520112794, 0.004966372455814472, 0.0028740581341518937, 0.0014715177646857699, 6.20796556976809E-4, 1.8393972058572123E-4, 2.2992465073215154E-5)
likelihood0000 <- c(0.16554574852714907, 0.14715177646857694, 0.1287578044100048, 0.1103638323514327, 0.09196986029286058, 0.07357588823428847, 0.05518191617571635, 0.036787944117144235, 0.018393972058572117)




insertProbability <- seq(0.9,0.1,by=-0.1)

data <- data.frame(probability=insertProbability,likelihood=likelihoodTree,internal0=likelihood0000,internal2=likelihood2000,internal21=likelihood2100)
#data <- pivot_longer(data,!probability,names_to = "clockRate",values_to = "likelihood")

likelihood <- ggplot(data=data,aes(x=insertProbability,y=likelihood)) + geom_line() + xlab("Insert 1 probability") + ylab("Likelihood (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")
state0 <- ggplot(data=data,aes(x=insertProbability,y=internal0)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 0000 (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")
state2 <- ggplot(data=data,aes(x=insertProbability,y=internal2)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 2000 (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")
state21 <- ggplot(data=data,aes(x=insertProbability,y=internal21)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 2100 (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")


plot_grid(likelihood,state0,state2,state21,nrow=1)



#with sequences 21100 2100
likelihoodTree <- exp(c(-6.411748421698742, -5.9870908979224895, -5.8833462549899, -5.940503564432857, -6.1205879097907125, -6.425375450145432, -6.889867436932458, -7.613255604161904, -8.93076419355177) )
likelihood0000 <- c(0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233, 0.36787944117144233)




insertProbability <- seq(0.9,0.1,by=-0.1)

data <- data.frame(probability=insertProbability,likelihood=likelihoodTree,internal0=likelihood0000,internal2=likelihood2000)
#data <- pivot_longer(data,!probability,names_to = "clockRate",values_to = "likelihood")

likelihood <- ggplot(data=data,aes(x=insertProbability,y=likelihood)) + geom_line() + xlab("Insert 1 probability") + ylab("Likelihood (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")
state0 <- ggplot(data=data,aes(x=insertProbability,y=internal0)) + geom_line() + xlab("Insert 1 probability") + ylab("Lik. 0000 (clockRate = 0.1)") + geom_vline(xintercept=0.5,colour="blue")


plot_grid(likelihood,state0,nrow=1)
