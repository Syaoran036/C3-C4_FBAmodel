# sensitivity analysis for PEPC activity

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
#FBAmodFile="Models/C4GEM_maize model_H2Oplasmo_C2.xml"
# read model, remove singletons()
# bndCond is set to TRUE by default, thus identifying external metabolites by 

#C4FBAmodel<-readSBMLmod(FBAmodFile, balanceReact=FALSE, mergeMet=FALSE, singletonMet=TRUE,remMet=TRUE)
#save(C4FBAmodel,file="C4GEM_maize model_H2Oplasmo_C2.rfile")

load("Models/C4GEM_maize model_H2Oplasmo_C2.rfile")

constFile.NADPME="Constraints/reactionConstraints_GS-GOGAT_limited.tsv"
C4FBAmodel.const.NADPME<-read.and.set.constraints(constFile.NADPME,C4FBAmodel)

# set obj function (here maximize mesophyll and bundle sheath biomass)
obj_coef(C4FBAmodel.const.NADPME)[match(c("BIO_m","BIO_b"),react_id(C4FBAmodel.const.NADPME))]<-0.055

##################################################################################################################

# source function for coupled model
# this creates vC.C3C4.FBA(pars,Cm,Om,Etot,mechModFile,FBAmodelorg)
source(file="Models/C2-C4_GS-GOGAT-limited.R")

# source Monte carlo Parameters, those include C,O and Etot
source("Parameters/MonteCarloParametersCmDep_FlaveriaRangesAliData712_equalSteps.R")

GOGAT.range<-seq(0,6,length=7)

res<-list(sol.stat=rep(NA,length=length(GOGAT.range)),C4fluxes=vector(mode="list",length=length(GOGAT.range) ), 
          plasmoFluxes= vector(mode="list",length=length(GOGAT.range) ),
          vCres=vector(mode="list",length=length(GOGAT.range) ),
          obj=rep(NA,length=length(GOGAT.range))
)

#       #plasmodesmata reactions:
plasmoPos<-grep("TMB",react_id(C4FBAmodel.const.NADPME))
# create vector of weights for mtfFBA
mtf.weights<-rep(1,length=length(react_id(C4FBAmodel.const.NADPME)))
# set a higher weight to plasmodesmatal reactions
mtf.weights[plasmoPos]<-1.1


for (i in seq_along(GOGAT.range)){
  
  parms<-c(beta=0.65, kcat=2.77, Vpmax=6, Kp=80, gs=15*1e-3, xi=0.58,GluS=GOGAT.range[i]) #Specific for F.ram Type II C3-C4
  # conduct coupled FBA (save and load file, variability analysis takes a long time)
  
  FVA=FALSE
  
  vCFBAres<-vC.C3C4.FBA(pars=parms,
                        Cm=C,
                        Om=O,
                        Etot=Etot,
                        mechModFile="Models/vC C3C4 Cm Landscape model Rubisco pars kcat dependent FAST.r",
                        FBAmodelorg=C4FBAmodel.const.NADPME,
                        alg="mtf",
                        #varFBA=FVA
                        costcoefMTF=mtf.weights # costcoefbw copies to costcoeffw
  )
  
  # get the solution status from the optSolObj, for glpkAPI, 5 translates to "solution is optimal"  getMeanStatus(5,"GLPKApi")
  res$sol.stat[i]<-lp_stat(vCFBAres$optSolObj)
  res$vCres[[i]]<-vCFBAres
  
  #save(vCFBAres,file="vCFBAres")
  
  #load("vCFBAres")
  
  if(FVA){
    # create matrix with cols react_name, minflux, maxflux
    varMinMax<-cbind(react_name(C4FBAmodel.const.NADPME),apply(getFluxDist(vCFBAres$varFBA),1,min),apply(getFluxDist(vCFBAres$varFBA),1,max) ) 
  }else{
    varMinMax<-cbind(react_name(C4FBAmodel.const.NADPME),rep(NA,length(react_name(C4FBAmodel.const.NADPME))),rep(NA,length(react_name(C4FBAmodel.const.NADPME)))) 
  }
  #vCFBAres
  
  #####################################################################################################################
  
  # evaluate results:
  vCFBAres.flux <- getFluxDist(vCFBAres$optSolObj)
  #       
  #       #######################################################################################
  #       # check some fluxes
  #       # exchange reactions:
  #       exchReact<-findExchReact(C4FBAmodel.const.NADPME)
  #       exchNames<-react_id(exchReact)
  #       exchNames<-c(exchNames,"REner01_p","REner01_pb","BIO_m","BIO_b") # add biomass to excahnge reactions, findExchReact() cant find it because multiple products are invloved
  #       exchpos<-c(react_pos(exchReact),match(c("REner01_p","REner01_pb","BIO_m","BIO_b"),react_id(C4FBAmodel.const.NADPME)))
  #       
  #       exchFluxes<-cbind(react_id(C4FBAmodel.const.NADPME)[exchpos],vCFBAres.flux[exchpos],varMinMax[,2][exchpos],varMinMax[,3][exchpos],c(met_id(exchReact),NA,NA,NA,NA))
  #       
  #,"H2O MtoBS"
  # mappinmg of plsma react ids to descriptions in matrix form
  plasmoNames<-cbind(c("TMBglu","TMBgln","TMBasn","TMBthr","TMBgly", "TMBser", "TMBaKG" ,"TMBh2o", "TMB01",  "TMB02",  "TMB03",  "TMB04",  "TMB05",  "TMB06",  "TMB07",  "TMB08",  "TMB09",  "TMB10",  "TMB11",  "TMB12", "TMB13") ,c("glu MtoBS","gln MtoBS","asn MtoBS","thr MtoBS","gly MtoBS","Ser MtoBs","aKG MtoBS","H2O MtoBS","Malate MtoBs","pyruvate MtoBS","Phosphoglycerate Mto BS","DHAP MtoB","GAP MtoBs","Pi MtoB","Pyrophosphate MtoB","Sucrose MtoB","Aspartate MtoB", "Alanine Mto B","CO2 MtoB","O2 MtoB","PEP M to B"))
  plasmoFluxes<-cbind(plasmoNames[match(react_id(C4FBAmodel.const.NADPME)[plasmoPos],plasmoNames[,1]),2],
                      plasmoNames[match(react_id(C4FBAmodel.const.NADPME)[plasmoPos],plasmoNames[,1]),1],
                      round(vCFBAres.flux[plasmoPos],digits=4),
                      varMinMax[,2][plasmoPos],
                      varMinMax[,3][plasmoPos])
  colnames(plasmoFluxes)<-c("name","react_id","flux","min","max")
  
  res$plasmoFluxes[[i]]<-data.frame(plasmoFluxes, row.names=plasmoFluxes[,1], stringsAsFactors=FALSE)
  
  
  #       # C4 cycle relevant reactions
  C4reacts<-c("Anet"="Ex1","Rubisco carb M"="R00024_p","Rubisco carb BS"="R00024_pb","Rubisco oxyg M"="R03140_p","Rubisco oxyg BS"="R03140_pb","PEPC M"="R00345_c","PEPC BS"="R00345_cb",
              "NADPME M"="R00216_p","NADPME BS"="R00216_pb",
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
              "IDH_MC"="R00709_m","IDH_BSC"="R00709_mb",
              "GDC_cb"="R01221_cb"
  )      
  
  #       # the model also defines Ala-AT and Asp-AT in plastids, but those are removed in singleton analysis.
  C4pos<-match(C4reacts,react_id(C4FBAmodel.const.NADPME))
  C4fluxes<-cbind("react"=names(C4reacts),"react_id"=react_id(C4FBAmodel.const.NADPME)[C4pos],"vCFBAres.flux"=vCFBAres.flux[C4pos],"varMin"=varMinMax[,2][C4pos],"varMax"=varMinMax[,3][C4pos])      
  
  res$C4fluxes[[i]]<-data.frame(C4fluxes,row.names=C4fluxes[,1],stringsAsFactors=FALSE)
  
  #       # objective function reactions:
          objPos<-which(obj_coef(C4FBAmodel.const.NADPME)!=0)
  #       objfluxes<-cbind(react_id(C4FBAmodel.const.NADPME)[objPos],vCFBAres.flux[objPos],varMinMax[,2][objPos],varMinMax[,3][objPos])
  #       
  #       # fluxes that hit constraints:
  #       # better to a minimization of total flux:
  #       mtfFBA<-optimizeProb(C4FBAmodel.const.NADPME,algorithm="mtf")
  #       mtfvCFBAres.flux<-as.numeric(getFluxDist(mtfFBA))
  #    
  
  res$obj[i]<-sum( vCFBAres.flux[which(obj_coef(C4FBAmodel.const.NADPME)!=0)] ) # 
}

# check if all solutions are optimal:
res$sol.stat

res$C4fluxes
res$plasmoFluxes

pdf("PLOT/GOGAT_ME_GDC.pdf",width=8.5*0.3937,height=6.5*0.3937)
# more space for y label
opar<-par(mar=c(5.1, 4.8, 2.1, 2.1),
          cex=0.6)
PEPC_M<-unlist(lapply(res$C4fluxes, function(x){as.numeric(x["PEPC M","vCFBAres.flux"])}))
Rubisco_oxyg_M<-unlist(lapply(res$C4fluxes, function(x){as.numeric(x["Rubisco oxyg M","vCFBAres.flux"])}))
NADPME_BS<-unlist(lapply(res$C4fluxes, function(x){as.numeric(x["NADPME BS","vCFBAres.flux"])}))
GDC_cb<-unlist(lapply(res$C4fluxes, function(x){as.numeric(x["GDC_cb","vCFBAres.flux"])}))
GDC_mito_M<-unlist(lapply(res$C4fluxes, function(x){as.numeric(x["GDC mito M","vCFBAres.flux"])}))
IDH_MC<-unlist(lapply(res$C4fluxes, function(x){as.numeric(x["IDH_MC","vCFBAres.flux"])}))
IDH_BS<-unlist(lapply(res$C4fluxes, function(x){as.numeric(x["IDH_BS","vCFBAres.flux"])}))

Gly_MtoBS<-unlist(lapply(res$plasmoFluxes, function(x){as.numeric(x["gly MtoBS","flux"])}))
Ala_MtoBS<-unlist(lapply(res$plasmoFluxes, function(x){as.numeric(x["Alanine Mto B","flux"])}))
Glu_MtoBS<-unlist(lapply(res$plasmoFluxes, function(x){as.numeric(x["glu MtoBS","flux"])}))
Asp_MtoBS<-unlist(lapply(res$plasmoFluxes, function(x){as.numeric(x["Aspartate MtoB","flux"])}))
Malate_MtoBS<-unlist(lapply(res$plasmoFluxes, function(x){as.numeric(x["Malate MtoBs","flux"])}))

plot(GOGAT.range,
     NADPME_BS, 
     xlim=c(0.5,5.5),
     ylim=c(min(min(NADPME_BS),min(GDC_cb)),max(max(NADPME_BS),max(GDC_cb))),
     xlab=expression("Fd-GOGAT activity [µmol"%*%"m"^-2%*%"s"^-1*"]"),
     ylab=expression("Plasmo. fluxes [µmol"%*%"m"^-2%*%"s"^-1*"]"),     
     #main="Black=AlaAT_m, Grey=Ala_AT_bs, blue=aKG_MtoB, green=glu_MtoBS",
     main=NULL,
     type="l"
)
lines(GOGAT.range, GDC_cb, col="grey")

#lines(GOGAT.range, NADPME_BS, col="grey")
#lines(GOGAT.range, -GDC_mito_BS, col="black",lty=2)
#lines(GOGAT.range, -GDC_mito_M,col="grey",lty=2)
#legend("topright",legend=c("PEPC M","NADPME BS","GDC mito BS","GDC mito M"),ncol=2,cex=0.8,bty="n",col=c("black","grey","black","grey"),lty=c(1,1,2,2))

legend("topleft",legend=c("NADP-ME (BSC)","GDC (BSC)"),ncol=1,cex=0.8,bty="n",col=c("black","grey"),lty=c(1,1))
dev.off()

