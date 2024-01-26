
# 29.1.13

#  funcion for coupling vc C3-C4 model with C4GEM-FBA:

# also changed constraints from original: - change fixed biomass and starch production (obvious)
#                            - change possibility of fatty acid uptake
#                           - allow starch production in M
#                           - allow reversibility of GOGAT enzymes

# varFBA  :   should a FVA be conducted and returned?

# 3.4.14:
# 

vC.C3C4.FBA_wo_C2<-function(pars,
                      Cm,
                      Om,
                      Etot,
                      mechModFile="Models/vC C3C4 Cm Landscape model Rubisco pars kcat dependent FAST.r",
                      FBAmodelorg,
                      alg="fba",  #algorithm parameter for optimizeProb
                      costcoefMTF=NULL, # weights for MTF, also used in fluxVarMTF
                      varFBA=FALSE,
                      varReactions=NULL, # passed to the "react" argument of fluxVarMTF
                      ...
                      ){
      
      require(sybilSBML)
      
      ########################################################
      # load vC model
      source(mechModFile)
     
      vCres<-AcC3C4lsKin(pars,Cm=Cm,Om=Om,Etot=Etot,details=TRUE)
      
      C4FBAmodel.const<-FBAmodelorg
      
      #constrain Rubisco Oxigenase reaction, R5BP + O2 -> 3PGA + Pglyc
      lowbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="R03140_p")] <- vCres$Vom
      ##
      lowbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="R03140_pb")] <- vCres$Vos
      
      # constrain Rubisco Carboxylase reactions
      uppbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="R00024_p")] <- vCres$Vcm
      
       uppbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="R00024_pb")] <- vCres$Vcs
#       
      # constrain PEPC flux
      uppbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="R00345_c")] <- -vCres$Vp
      # force decarboxylation in BS by NADPME
      lowbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="R00216_pb")] <- vCres$Vp
      
      #constrain CO2 uptake to net assimilation rate
      uppbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="Ex1")] <- vCres$Ac
      lowbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="Ex1")] <- vCres$Ac
      
      #reduce CO2 diffusion from ms to bs to  -Leakage from bs to ms
      #uppbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="TMB11")] <- -vCres$L
      #lowbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="TMB11")] <- -1000
      
      
#       #constrain Glycine and Serine plasmodesmata exchange 
      #lowbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="TMBgly")] <- vCres$xi * vCres$Vom
      
      #uppbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="TMBser")] <- -0.5 * vCres$xi * vCres$Vom
      
      #glycine decarboxylase:
      lowbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="R01221_mb")] <- 0.5 * vCres$xi * vCres$Vom
      uppbnd(C4FBAmodel.const)[which(react_id(C4FBAmodel.const)=="R01221_m")] <- 0.5 * (1-vCres$xi) * vCres$Vom
      
      
      #######################################################################################
      
      # conduct the FBA:
      
      if(alg=="fba"){
              optSol.inst<-optimizeProb(C4FBAmodel.const,
                                            algorithm=alg,
                                            #costcoefbw=costcoefMTF,# costcoefbw copies to costcoeffw
                                            ...)
      }else if(alg=="mtf"){
            if(is.null(costcoefMTF)){
                costcoefMTF<-rep(1,react_num(C4FBAmodel.const))# equal weights if costcoefbw was not in "..."
            } 
            optSol.inst<-optimizeProb(C4FBAmodel.const,
                                      algorithm=alg,
                                      costcoefbw=costcoefMTF,# costcoefbw copies to costcoeffw
                                      ...)
      }
      
      if(varFBA==TRUE & alg=="fba"){
        
            var<-fluxVar(C4FBAmodel.const,fld=TRUE) #only include fluxes that correspond to MFA
            return(list("lp_ok"=lp_ok(optSol.inst),"Ac_Obj"=c("Ac"=vCres$Ac,"Obj"=lp_obj(optSol.inst)),"vCres"=vCres,"vCpars"=pars,"optSolObj"=optSol.inst,"varFBA"=var))
            
      }else if (varFBA==TRUE & alg=="mtf"){
        
            if(is.null(varReactions)){
                  varReactions<-1:react_num(C4FBAmodel.const) # perform fv for all reactions
            }
            
            # source sysBiolAlg class for mtffv
            source("mtffv/sysBiolAlg_mtffvClass.R")
            # source fluxVarMTF()
            source("mtffv/fluxVarMTF.R")
            
            var<-fluxVarMTF(C4FBAmodel.const,
                          react= varReactions,
                         wtobj = mod_obj(optSol.inst), # the dotproduct of obj coefficients and model fluxes
                         sumAbsFlux = sum(costcoefMTF*abs(getFluxDist(optSol.inst))),
                         costcoefSumAbsFlux = costcoefMTF,
                         fld="fluxes" # argument to optimizer(), return only model fluxes
            )
            
            return(list("lp_ok"=lp_ok(optSol.inst),"Ac_Obj"=c("Ac"=vCres$Ac,"Obj"=lp_obj(optSol.inst)),"vCres"=vCres,"vCpars"=pars,"optSolObj"=optSol.inst,"varFBA"=var))
      }else{
            return(list("lp_ok"=lp_ok(optSol.inst),"Ac_Obj"=c("Ac"=vCres$Ac,"Obj"=lp_obj(optSol.inst)),"vCres"=vCres,"vCpars"=pars,"optSolObj"=optSol.inst))
      }
}
