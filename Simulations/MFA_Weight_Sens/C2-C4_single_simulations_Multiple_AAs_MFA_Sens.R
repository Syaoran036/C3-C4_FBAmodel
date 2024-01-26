setwd("Simulations/MFA_Weight_Sens/")

source("Simulations/MFA_Weight_Sens/vC_FBA_weight_Sens_Function.R")

sens_res<-vC_FBA_weight_Sens(
                  constraint_file="Constraints/reactionConstraintsBMobj_NADPMEC2-C4_Multiple_AAs.tsv",
                  weight_vec=c(1, 1.1, 2, 3, 4, 5,)#c(1, 1.1, 5, 10, 50, 100),
        )

save(sens_res,file="PlasmoWeightSensResult_Multiple_AAs_1to5")

# Analyze results for w=5, here a aKG/Gln Shuttle is used:
raw<-sens_res[[2]]
vCFBAres<-raw[[6]]
vCFBAres.flux <- getFluxDist(vCFBAres$optSolObj)
# function does not return modelOrg Obj, so reload, but make sure it's the same one used earlier!!
load("Models/C4GEM_maize model_H2Oplasmo_C2.rfile")
read.and.set.constraints<-function(file,model){
  constraints<-read.csv(file,sep="\t",strip.white=TRUE)
  #positions of model ids to change
  model.pos<-match(constraints$reactionID,react_id(model))
  if(any(is.na(model.pos))){stop("there are react_ids in the constraint file that don't match the model ids!")}
  uppbnd(model)[model.pos]<-constraints$ub
  lowbnd(model)[model.pos]<-constraints$lb
  return(model)
}

C4FBAmodel.const.NADPME<-read.and.set.constraints("Constraints/reactionConstraintsBMobj_NADPMEC2-C4_Multiple_AAs.tsv",C4FBAmodel)


varMinMax<-cbind(react_name(C4FBAmodel.const.NADPME),apply(getFluxDist(vCFBAres$varFBA),1,min),apply(getFluxDist(vCFBAres$varFBA),1,max) ) 

# return reactionids of reactions using met (String)
get.Met.react_id<- function(mod,met){
  return(react_id(mod)[which(S(mod)[which(met_id(mod)==met),]!=0)])
}
# return reaction names of reactions using met (String)
get.Met.react_name<- function(mod,met){
  return(react_name(mod)[which(S(mod)[which(met_id(mod)==met),]!=0)])
}
# return flux by react_id
get.flux.by.react_id<- function(mod,ids,fluxdist,min=rep(NA,length(fluxdist)),max=rep(NA,length(fluxdist))){
  return( cbind("ids"=ids,
                "names"=strtrim(react_name(mod)[match(ids,react_id(mod))],30) ,
                "fluxes"=round(fluxdist[match(ids,react_id(mod))],digits=4),
                "min"= min[match(ids,react_id(mod))],
                "max"= max[match(ids,react_id(mod))]
  ) )
}


#"S_L_45_Glutamine_cb"
Gln.cb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "L_45_Glutamine[cb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Gln.cb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#L_45_Glutamate_cb
Glu.cb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "L_45_Glutamate[cb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Glu.cb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#NH3_cb
NH3.cb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "NH3[cb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, NH3.cb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )


#"S_L_45_Glutamine_c"
Gln.c.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "L_45_Glutamine[c]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Gln.c.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )
