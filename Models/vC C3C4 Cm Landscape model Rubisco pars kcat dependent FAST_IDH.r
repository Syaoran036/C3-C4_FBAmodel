# C3-C4 Steady State von Caemmerer model as described in von Caemmerer 2000 book
# included dependency of Kc, Ko, Kc/Ko and gammaAst on kcat.

# edit: 23.1.13, added details option

# pars: named vector

#Cm = 300,
#Om = 200*1e3,
#total amount of RubisCO catalytic sites [?mol/sqm]
#Etot = 25.8

AcC3C4lsKin <- function (pars,
                  Cm = NULL,
                  Om = NULL,
                  #total amount of RubisCO catalytic sites [?mol/sqm]
                  Etot = NULL,   # 80 ?mol/(sqm s) / 3.1 1/s = Vcmax(C3) / kcat (F.pringlei, Kubien 2009)
                  details=FALSE       
                  ){
  
                  if(missing(Cm)|missing(Om)|missing(Etot)){stop("parameters missing!")}
                  # carboxylation rate constant [1/s]:
                  kcat = pars["kcat"]

                  # beta: fraction of RubisCO expressed in mesophyll (beween 0 and 1)
                  beta = pars["beta"]
                  # max RubisCO activity for the whole leaf:


                  # max. mesophyll RubisCO activity [?mol / (sqm * s) ]
                  #Vmmax = pars["Vmmax"]
                  Vmmax = beta * Etot * kcat
                  # max. bundle sheath RubisCO activity [?mol / (sqm * s) ]
                  #Vsmax = pars["Vsmax"]
                  Vsmax = (1-beta) * Etot * kcat
                  #Vcmax,needed for Rd:
                  Vcmax <- Vmmax + Vsmax

                  ## MM const of RubsiCO for CO2 [?bar]
#                  # 2dim abs(residuals) fit of Savir data (w/o form II)
#                  Kc = 23.36 * kcat^2.22
#                  # Kc/Ko
#                  # 2dim least squares fit of Savir data (w/o form II)
#                  Kc.Ko = 3e-4  * kcat^1.37
#                  # half the reiprocal rubisco specificity, dimensionless
#                  # (capital) GammaAst/Os = 0.5/Sco = sGammaAst
#                  # Sco was estmated with 2dim least squares fit of Savir data (w/o form II)
#                  sGammaAst = 0.5/ (5011.62* kcat^-0.60)

                  # MM const of RubsiCO for CO2 [?bar]
                  # 2dim least squares fit of Savir data (w/o form II & Synechococcus)
                  Kc = 16.07 * kcat^2.36
                  # Kc/Ko
                  # 2dim least squares fit of Savir data (w/o form II & Synechococcus)
                  Kc.Ko = 0.00037 * kcat^1.16
                  # half the reiprocal rubisco specificity, dimensionless
                  # (capital) GammaAst/Os = 0.5/Sco = sGammaAst
                  # 2dim least squares fit of Savir data (w/o form II & Synechococcus)
                  sGammaAst = 0.5/ (5009.76* kcat^-0.6)

                  # Max PEPC activity [?mol / (sqm * s)]
                  Vpmax = pars["Vpmax"]
                  # MM constant for PEPC [?bar]
                  # Bauwe 1986
                  Kp = pars["Kp"]
                  # BS CO2 conductance [?mol /(sqm*s)]
                  # chosen (arbitrarily ??) in between literature values
                  gs = pars["gs"]
                  # Fraction of mesophyll photorespiration occuring in the BS
                  # in Moricandia all GDC in BS, similar results for Flaveria and Panicum
                  xi = pars["xi"]

                  # leaf mitochondrial respiration   [?mol / (sqm * s) ]
                  # scaled with RubsiCO, based on loose correlation betwen Rd and leaf protein content
                  # rather arbitrarily
                  if("Rd" %in% names(pars) ){
                      Rd <- pars["Rd"]
                      }else{
                      Rd <- 0.01 * Vcmax
                  }

                  # Mesophyll mitochondrial respiration   [mol / (sqm * s) ]
                  # higher than in C4 model
                  # rather arbitrarily
                  if("Rm" %in% names(pars) ){
                      Rm <- pars["Rm"]
                      }else{
                      Rm <- 0.8 * Rd
                  }


                  if("Rs" %in% names(pars)){
                      Rs <- pars["Rs"]
                      }else{
                      Rs <- Rd - Rm
                  }



            # BS O2 conductance [mol /(sqm*s)]
            # related to gs by ratio of diffusivities and solubilities
            # Farquhar 1983 ,Berry and Farquhar 1978
            gO <- 0.047 * gs

            # PEPC reaction
            Vp <- (Cm*Vpmax)/(Cm+Kp)

            # rate of RubisCO carboxylation in the mesophyll
            Vcm <- (Cm * Vmmax) / (Cm +  (Kc + Om*Kc.Ko))

            # rate of RubisCO oxygenations in the mesophyll (see docx)
            Vom <- (2*sGammaAst*Om*Vcm) / Cm
            #Vom <- (2* sGammaAst * Om * Vmmax)/(Cm + Kc * (1+(Om/Ko)))

            ####################################

            a <- 1 - (1/0.047) * Kc.Ko

            b <- -( (Vp + xi * 0.5 * Vom + gs * Cm) + (Vsmax - Rs) + gs * (Kc + Om*Kc.Ko) + (sGammaAst * Vsmax + ((Kc.Ko) * Rs))/0.047 )

            c <- (Vsmax - Rs) * (Vp + xi * 0.5 * Vom + gs*Cm) - (Vsmax*gs*sGammaAst*Om + Rs*gs*(Kc + Om*Kc.Ko))



            #a <- -( -1 + (Kc/(Ko*0.047)) )
            #b <- -( (gs*Cm + Vp + xi * 0.5 * Vom )  + gs*Kc + (gs * Om * Kc/Ko) - Rs + (Rs*Kc)/(0.047*Ko) + Vsmax + (sGammaAst * Vsmax)/0.047 )
            #c <- -( (gs*Cm + Vp + xi * 0.5 * Vom ) * (Rs - Vsmax) + (Kc * gs * Rs) * (1 + Om/Ko) + gs * sGammaAst* Om * Vsmax )

            ####################################

            As <- (-b - sqrt(b^2 - 4*a*c)) / (2*a)

            Am <- ( ( (Cm - sGammaAst * Om) * Vmmax ) / (Cm + (Kc + Om*Kc.Ko) ) ) - Rm

            Ac <- As + Am

            if(!details){
                return(Ac)
            }else{      
                # equation for Cs:
                
                # L = gs (Cs - Cm)
                # As = Vp + xi * 0.5 * Vom - L
                # As = Vp + xi * 0.5 * Vom -  gs (Cs - Cm)
                # -As = -Vp - xi  0.5  Vom +  gs (Cs - Cm)
                #  -As + Vp + xi  0.5  Vom =  gs (Cs - Cm)
                # ((-As + Vp + xi  0.5  Vom)/gs ) + Cm = Cs
                
                
                Cs <- Cm + (Vp + xi * 0.5 * Vom - As) /gs
                #alternativ:
                #Cs <- ( ( -Ac + Vcm - (1-xi)*0.5*Vom + Vp - Rm) / gs) + Cm
                
                # see docx for explanation, assumption is that all Calvin Cycle bpg
                # is reduced in the BS and creates the respective O2 evolution:
                Os <- As/(0.047 * gs) + Om
                
                Vcs <- (Cs*Vsmax)/(Cs+(Kc + Os*Kc.Ko))
                
                Vos <- Vcs*2*sGammaAst*Os/Cs
                
                # leakiness
                L <- gs*(Cs-Cm)    
                
                #GOGAT
                IDH <- pars["IDH"]
                
                result<- list(c())
                  
                  result[[1]] <- Cm
                  result[[2]] <- Cs
                  result[[3]] <- Os
                  result[[4]] <- Ac
                  result[[5]] <- Am
                  result[[6]] <- As
                  result[[7]] <- Vp
                  result[[8]] <- Vcm
                  result[[9]] <- Vom
                  result[[10]] <- Vcs
                  result[[11]] <- Vos
                  result[[12]] <- L
                  result[[13]] <- gs
                  result[[14]] <- xi
                  result[[15]] <- Vsmax
                  result[[16]] <- IDH
                  
                  names(result)<-c("Cm","Cs","Os","Ac","Am","As","Vp","Vcm","Vom","Vcs","Vos","L","gs","xi","Vsmax","IDH")
                  result <- as.data.frame(result)
                                    
                  return(result) 
            }                
}