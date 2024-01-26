
###################################################################################################################

# prepare the modelorg object (load and set constraints from file) 

require(sybilSBML)
require(sybil)

# read and set constraints
read.and.set.constraints<-function(file,model){
  constraints<-read.csv(file,sep="\t",strip.white=TRUE)
  #positions of model ids to change
  model.pos<-match(constraints$reactionID,react_id(model))
  if(any(is.na(model.pos))){stop("there are react_ids in the constraint file that don't match the model ids!")}
  uppbnd(model)[model.pos]<-constraints$ub
  lowbnd(model)[model.pos]<-constraints$lb
  return(model)
}


########
#NADPME:

# JBF: adjust to joined biomass function
FBAmodFile="Models/C4GEM_maize model_H2Oplasmo_C2.xml"
constFile.NADPME="Constraints/reactionConstraintsBMobj_NADPMEC2-C4_Multiple_AAs_only_Ser_Gly.tsv"
# read model, remove singletons()
# bndCond is set to TRUE by default, thus identifying external metabolites by 
#C4FBAmodel<-readSBMLmod(FBAmodFile, balanceReact=FALSE, mergeMet=FALSE, singletonMet=TRUE,remMet=TRUE)

# Warning messages:
# 1: 8 metabolites are not used in any reaction and therefore removed from S:
# ...
# 2: 1 reaction is not used and therefore removed from S: “R_BIO_LEAF” 
# 3: removed 1740 singleton metabolites:
# ...
# 4: removed 1352 reactions containing singleton metabolites:
# ...

#save(C4FBAmodel,file="C4GEM_maize model_H2Oplasmo_C2.rfile")
#load("Models/C4GEM_maize model_H2Oplasmo_C2.rfile")

C4FBAmodel.const.NADPME<-read.and.set.constraints(constFile.NADPME,C4FBAmodel)

# set obj function (here maximize mesophyll and bundle sheath biomass)
obj_coef(C4FBAmodel.const.NADPME)[match(c("BIO_m","BIO_b"),react_id(C4FBAmodel.const.NADPME))]<-0.055
# JBF: adjust to joined biomass function
#obj_coef(C4FBAmodel.const.NADPME)[match(c("BIO_m"),react_id(C4FBAmodel.const.NADPME))]<-0.055


###################################################################################################################

# source function for coupled model
# this creates vC.C3C4.FBA(pars,Cm,Om,Etot,mechModFile,FBAmodelorg)
source(file="Models/vC_C4GEM_coupled_NADPME.R")

# source Monte carlo Parameters, those include C,O and Etot
source("Parameters/MonteCarloParametersCmDep_FlaveriaRangesAliData712_equalSteps.R")

parms<-c(beta=rang$beta[1],kcat=rang$kcat[1],Vpmax=rang$Vpmax[1],Kp=rang$Kp[1],gs=rang$gs[1],xi=rang$xi[2])
#parms<-c(beta=rang$beta[1],kcat=rang$kcat[1],Vpmax=6,Kp=rang$Kp[1],gs=rang$gs[1],xi=rang$xi[2])

#plasmodesmata reactions:
plasmoPos<-grep("TMB",react_id(C4FBAmodel.const.NADPME))
# create vector of weights for mtfFBA
mtf.weights<-rep(1,length=length(react_id(C4FBAmodel.const.NADPME)))
# set a higher weight to plasmodesmatal reactions
mtf.weights[plasmoPos]<-1.1

# conduct coupled FBA (save and load file, variability analysis takes a long time)

FVA=TRUE

#       #plasmodesmata reactions:
#plasmoPos<-grep("TMB",react_id(C4FBAmodel.const.NADPME))
# create vector of weights for mtfFBA
#mtf.weights<-rep(1,length=length(react_id(C4FBAmodel.const.NADPME)))
# set a higher weight to plasmodesmatal reactions
#mtf.weights[plasmoPos]<-1.1

vCFBAres<-vC.C3C4.FBA(pars=parms,
                      Cm=C,
                      Om=O,
                      Etot=Etot,
                      mechModFile="Models/vC C3C4 Cm Landscape model Rubisco pars kcat dependent FAST.r",
                      FBAmodelorg=C4FBAmodel.const.NADPME,
                      alg="mtf",
                      varFBA=FVA,
                      #costcoefbw=mtf.weights # costcoefbw copies to costcoeffw
                      costcoefMTF=mtf.weights
)

save(vCFBAres,file="vCFBAres_Multiple_AAs_only_Ser_Gly")

#load("vCFBAres")

if(FVA){
  # create matrix with cols react_name, minflux, maxflux
  varMinMax<-cbind(react_name(C4FBAmodel.const.NADPME),apply(getFluxDist(vCFBAres$varFBA),1,min),apply(getFluxDist(vCFBAres$varFBA),1,max) ) 
}else{
  varMinMax<-cbind(react_name(C4FBAmodel.const.NADPME),rep(NA,length(react_name(C4FBAmodel.const.NADPME))),rep(NA,length(react_name(C4FBAmodel.const.NADPME)))) 
}
vCFBAres

#check solver status
all(vCFBAres$varFBA@lp_ok==0)
all(vCFBAres$varFBA@lp_stat==5)

#####################################################################################################################

# evaluate results:
vCFBAres.flux <- getFluxDist(vCFBAres$optSolObj)

#######################################################################################
# check some fluxes
# exchange reactions:
exchReact<-findExchReact(C4FBAmodel.const.NADPME)
exchNames<-react_id(exchReact)
exchNames<-c(exchNames,"REner01_p","REner01_pb","BIO_m","BIO_b") # add biomass to excahnge reactions, findExchReact() cant find it because multiple products are invloved
exchpos<-c(react_pos(exchReact),match(c("REner01_p","REner01_pb","BIO_m","BIO_b"),react_id(C4FBAmodel.const.NADPME)))

exchFluxes<-cbind(react_id(C4FBAmodel.const.NADPME)[exchpos],vCFBAres.flux[exchpos],varMinMax[,2][exchpos],varMinMax[,3][exchpos],c(met_id(exchReact),NA,NA,NA,NA))

#plasmodesmata reactions:
plasmoPos<-grep("TMB",react_id(C4FBAmodel.const.NADPME))#,"H2O MtoBS"
# mappinmg of plsma react ids to descriptions in matrix form
plasmoNames<-cbind(c("TMBglu","TMBgln","TMBasn","TMBthr","TMBgly", "TMBser", "TMBaKG" ,"TMBh2o", "TMB01",  "TMB02",  "TMB03",  "TMB04",  "TMB05",  "TMB06",  "TMB07",  "TMB08",  "TMB09",  "TMB10",  "TMB11",  "TMB12", "TMB13") ,c("glu MtoBS","gln MtoBS","asn MtoBS","thr MtoBS","gly MtoBS","Ser MtoBs","aKG MtoBS","H2O MtoBS","Malate MtoBs","pyruvate MtoBS","Phosphoglycerate Mto BS","DHAP MtoB","GAP MtoBs","Pi MtoB","Pyrophosphate MtoB","Sucrose MtoB","Aspartate MtoB", "Alanine Mto B","CO2 MtoB","O2 MtoB","PEP M to B"))
plasmoFluxes<-cbind(plasmoNames[match(react_id(C4FBAmodel.const.NADPME)[plasmoPos],plasmoNames[,1]),2],
                    plasmoNames[match(react_id(C4FBAmodel.const.NADPME)[plasmoPos],plasmoNames[,1]),1],
                    round(vCFBAres.flux[plasmoPos],digits=4),
                    varMinMax[,2][plasmoPos],
                    varMinMax[,3][plasmoPos])
colnames(plasmoFluxes)<-c("name","react_id","flux","min","max")

# C4 cycle relevant reactions
C4reacts<-c("Rubisco carb M"="R00024_p","Rubisco carb BS"="R00024_pb","PEPC M"="R00345_c","PEPC BS"="R00345_cb","NAD-ME M"="R00214_m","NAD-ME BS"="R00214_mb",
            "PEPCK M"="R00341_c","PEPCK BS"="R00341_cb","NADPME M"="R00216_p","NADPME BS"="R00216_pb",
            "Malate DH NADP plastid M"="R00343_p","Malate DH NADP mito M"="R00343_m","Malate DH NADP plastid b"="R00343_pb","Malate DH NADP mito b"="R00343_mb",
            "Malate DH NAD plastid M"="R00342_p","Malate DH NAD mito M"="R00342_m","Malate DH NAD plastid b"="R00342_pb","Malate DH NAD mito b"="R00342_mb",
            "PPDK plastid M"="R00206_p","GDC mito M"="R01221_m","GDC mito BS"="R01221_mb",
            "Ala-AT_cytop"="R00258_c","Ala-AT_cytop_bs"="R00258_cb",
            "Asp-AT_cytop"="R00355_c","Asp-AT_mito"="R00355_m","Asp-AT_cytop_bs"="R00355_cb","Asp-AT_mito_bs"="R00355_mb",
            "Glu DH mb"="R00243_mb",
            "GluS_Fd-dep_p"="R00021_p","GluS_Fd-dep_pb"="R00021_pb",
            "GluS_NADPH_dep_c"="R00114_c","GluS_NADPH_dep_cb"="R00114_cb",
            "GluS_NADH_dep_p"="R00093_p","GluS_NADH_dep_pb"="R00093_pb",
            "Gln_Synthetase_pb"="R00253_pb",
            "Glu_NH4_ligase_pb"="R00243_pb"
)

# the model also defines Ala-AT and Asp-AT in plastids, but those are removed in singleton analysis.
C4pos<-match(C4reacts,react_id(C4FBAmodel.const.NADPME))
C4fluxes<-cbind(names(C4reacts),react_id(C4FBAmodel.const.NADPME)[C4pos],vCFBAres.flux[C4pos],varMinMax[,2][C4pos],varMinMax[,3][C4pos])



# objective function reactions:
objPos<-which(obj_coef(C4FBAmodel.const.NADPME)!=0)
objfluxes<-cbind(react_id(C4FBAmodel.const.NADPME)[objPos],vCFBAres.flux[objPos],varMinMax[,2][objPos],varMinMax[,3][objPos])

# fluxes that hit constraints:
# better to a minimization of total flux:
mtfFBA<-optimizeProb(C4FBAmodel.const.NADPME,algorithm="mtf")
mtfvCFBAres.flux<-as.numeric(getFluxDist(mtfFBA))
hitConstPos<-which((abs(mtfvCFBAres.flux-uppbnd(C4FBAmodel.const.NADPME))<1e-4 | (abs(mtfvCFBAres.flux-lowbnd(C4FBAmodel.const.NADPME))<1e-4) ) & mtfvCFBAres.flux!=0)
cbind(react_id(C4FBAmodel.const.NADPME)[hitConstPos],react_name(C4FBAmodel.const.NADPME)[hitConstPos],mtfvCFBAres.flux[hitConstPos],lowbnd(C4FBAmodel.const.NADPME)[hitConstPos],uppbnd(C4FBAmodel.const.NADPME)[hitConstPos])

vCFBAres
exchFluxes
plasmoFluxes
C4fluxes
objfluxes


###############################################################################


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

bsreacts.ids<-grep(".*b$",react_id(C4FBAmodel.const.NADPME),value=TRUE)
bs.flux<-get.flux.by.react_id(C4FBAmodel.const.NADPME,bsreacts.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

bs.flux.relevant=bs.flux[which(abs(as.numeric(bs.flux[,3]))>0.5),]


#Serine
Ser.ids<-get.Met.react_id(C4FBAmodel.const.NADPME,"L_45_Serine[cb]")
#get.Met.react_name(C4FBAmodel.const.NADPME,"L_45_Serine[cb]")
get.flux.by.react_id(C4FBAmodel.const.NADPME, Ser.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

# P-Serine
pSer.ids<-get.Met.react_id(C4FBAmodel.const.NADPME,"O_45_Phospho_45_L_45_serine[cb]")
get.flux.by.react_id(C4FBAmodel.const.NADPME, pSer.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

# phosphonooxypyruvate
pop.ids<-get.Met.react_id(C4FBAmodel.const.NADPME,"3_45_Phosphonooxypyruvate[cb]")
get.flux.by.react_id(C4FBAmodel.const.NADPME, pop.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#3PGA bs cyto
pga.ids<-get.Met.react_id(C4FBAmodel.const.NADPME,"3_45_Phospho_45_D_45_glycerate[cb]")
get.flux.by.react_id(C4FBAmodel.const.NADPME, pga.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#BPG bs cyto
bpg.ids<-get.Met.react_id(C4FBAmodel.const.NADPME,"3_45_Phospho_45_D_45_glyceroyl_32_phosphate[cb]")
get.flux.by.react_id(C4FBAmodel.const.NADPME, bpg.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

#GAP bs cyto
gap.ids<-get.Met.react_id(C4FBAmodel.const.NADPME,"40_2R_41__45_2_45_Hydroxy_45_3_45__40_phosphonooxy_41__45_propanal[cb]")
get.flux.by.react_id(C4FBAmodel.const.NADPME, gap.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )


#FBP bs cyto
fbp.ids<-get.Met.react_id(C4FBAmodel.const.NADPME,"beta_45_D_45_Fructose_32_1_44_6_45_bisphosphate[cb]")
get.flux.by.react_id(C4FBAmodel.const.NADPME, fbp.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )



#F6P bs cyto
f6p.ids<-get.Met.react_id(C4FBAmodel.const.NADPME,"beta_45_D_45_Fructose_32_6_45_phosphate[cb]")
get.flux.by.react_id(C4FBAmodel.const.NADPME, f6p.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )


# cyto Glutamate
Glu.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "L_45_Glutamate[cb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Glu.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )



#all from mito export 
#"L_45_Glutamate_m"
Glu.mito.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "L_45_Glutamate[mb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, Glu.mito.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

# mitochondrial alpha KG
aKG.mito.ids <- get.Met.react_id( C4FBAmodel.const.NADPME, "2_45_Oxoglutarate[mb]" )
get.flux.by.react_id(C4FBAmodel.const.NADPME, aKG.mito.ids, vCFBAres.flux, min= varMinMax[,2], max= varMinMax[,3] )

# results suggest that mitochondrial Glu-DH (R00243_mb) is refixing NH3 to Glu, which is converted to Serine.


##################################################################
# sum of flux over plasmodesmata
sum(abs(as.numeric(plasmoFluxes[,"flux"])))
