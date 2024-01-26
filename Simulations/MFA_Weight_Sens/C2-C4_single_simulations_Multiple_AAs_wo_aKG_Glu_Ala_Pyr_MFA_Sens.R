source("Simulations/MFA_Weight_Sens/vC_FBA_weight_Sens_Function.R")



sens_res<-vC_FBA_weight_Sens(
  constraint_file="Constraints/reactionConstraintsBMobj_NADPMEC2-C4_Multiple_AAs_wo_aKG_Glu_Ala_Pyr.tsv",
  weight_vec=c(1, 1.1, 2, 3, 4, 5)#c(1, 1.1, 5, 10, 50, 100)
)

#save(sens_res,file="PlasmoWeightSensResult_wo_aKG_Glu_Ala_Pyr_1to5")

# analyze the Asn/Mal Shutlle in more detail, at 10 it's a unique solution:
#load("PlasmoWeightSensResult_wo_aKG_Glu_Ala_Pyr")

#plasmo.res[["10"]]

vCFBAres<-raw.res[["1.1"]]
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

C4FBAmodel.const.NADPME<-read.and.set.constraints("Constraints/reactionConstraintsBMobj_NADPMEC2-C4_Multiple_AAs_wo_aKG_Glu_Ala_Pyr.tsv",C4FBAmodel)


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


#"40_S_41__45_Malate_cb"
Mal.cb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "40_S_41__45_Malate[cb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Mal.cb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#"40_S_41__45_Malate_pb"
Mal.pb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "40_S_41__45_Malate[pb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Mal.pb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )


#"40_S_41__45_Malate_cb"
Mal.mb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "40_S_41__45_Malate[mb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Mal.mb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )


#Oxaloacetate_pb
OAA.mb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "Oxaloacetate[pb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, OAA.mb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#Oxaloacetate_cb
OAA.cb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "Oxaloacetate[cb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, OAA.cb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#L_45_Aspartate_cb
Asp.cb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "L_45_Aspartate[cb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Asp.cb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#L_45_Asparagine_cb
Asn.cb.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "L_45_Asparagine[cb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Asn.cb.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#L_45_Asparagine_c
Asn.c.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "L_45_Asparagine[c]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Asn.c.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#L_45_Aspartate_c
Asp.c.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "L_45_Aspartate[c]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Asp.c.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#Oxaloacetate_c
OAA.c.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "Oxaloacetate[c]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, OAA.c.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

# PEPC !?

#Phosphoenolpyruvate_cb
PEP.c.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "Phosphoenolpyruvate[c]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, PEP.c.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

