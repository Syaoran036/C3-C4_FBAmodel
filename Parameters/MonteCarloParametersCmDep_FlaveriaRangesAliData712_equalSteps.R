# 04/02/2013
# use equal amount of steps for all parameters (ie. 5 steps for xi)


# parameter ranges:
# rang <- list(
#   # max. mesophyll RubisCO activity [?mol / (sqm * s) ]
#   #Vmmax = c(80, 0),
#   # max. bundle sheath RubisCO activity [?mol / (sqm * s) ]
#   #Vsmax = c(0, 60),
#   #fraction of RubisCO expressed in mesophyll
#   beta = c(0.9,0),
#   kcat=c(3.4,8.8),
#   # Max PEPC activity [?mol / (sqm * s)]
#   Vpmax = c(0, 120),
#   # MM constant for PEPC [?bar]
#   # Bauwe 1986
#   Kp = c(2.5*80, 80),
#   # BS CO2 conductance [?mol /(sqm*s)]
#   # chosen (arbitrarily ??) in between literature values
#   gs = c(15*1e-3, 1*1e-3),
#   # Fraction of mesophyll photorespiration occuring in the BS
#   # in Moricandia all GDC in BS, similar results for Flaveria and Panicum
#   xi = c(0, 1)#,
#   #Rd=c(0,0)
# )


# # use Falveria data to estimate ranges:
source("Parameters/Flaveria Parameter Collection.r")
rang <- list(
  # max. mesophyll RubisCO activity [?mol / (sqm * s) ]
  #Vmmax = c(80, 0),
  # max. bundle sheath RubisCO activity [?mol / (sqm * s) ]
  #Vsmax = c(0, 60),
  #fraction of RubisCO expressed in mesophyll
  beta = c(max(beta,na.rm=TRUE),min(beta,na.rm=TRUE)),
 # beta = c(0.9,min(beta,na.rm=TRUE)),
  kcat=c(min(kcatWes,na.rm=TRUE),max(kcatWes,na.rm=TRUE)),
  # Max PEPC activity [?mol / (sqm * s)]
  Vpmax = c(0, 130),
#  Vpmax = c(0, 120),
  # MM constant for PEPC [?bar]
  # Bauwe 1986
  Kp = c(2.5*80, 80),
  # BS CO2 conductance [?mol /(sqm*s)]
  # chosen (arbitrarily ??) in between literature values
  gs = c(15*1e-3, 1*1e-3),
  # Fraction of mesophyll photorespiration occuring in the BS
  # in Moricandia all GDC in BS, similar results for Flaveria and Panicum
  xi = c(min(xiTranscriptome2,na.rm=TRUE), max(xiTranscriptome2,na.rm=TRUE))#,
  #Rd=c(0,0)
)

C<-400
O<-200e3

steps<-rep(5,length(rang))
# only 1 step for xi:
#steps[length(rang)]<-2
# steps for Vpmax:
#steps[3]<-10

#initialise vector of mutational probabilities
#C4prob <- c(beta=0.1,kcat=0.2,Vpmax=0.05,Kp=0.1,gs=0.1,xi=0.8)
C4prob <- c(beta=0.1,kcat=0.2,Vpmax=0.05,Kp=0.1,gs=0.1,xi=3.2)
C4prob <- C4prob/sum(C4prob)

N<-1e5
Etot<-19.35