


vC_FBA_weight_Sens<-function(constraint_file, 
                             weight_vec,
                             FBAmod_file="Models/C4GEM_maize model_H2Oplasmo_C2.rfile"
                             ){
      
      if(missing(constraint_file)){stop("provide constraint file")}
      if(missing(weight_vec)){stop("provide weight_vec")}
  
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
      constFile.NADPME=constraint_file
      # read model, remove singletons()
      # bndCond is set to TRUE by default, thus identifying external metabolites by 
      
      #save(C4FBAmodel,file="C4GEM_maize model_H2Oplasmo_C2.rfile")
      load(FBAmod_file)
      
      C4FBAmodel.const.NADPME<-read.and.set.constraints(constFile.NADPME,C4FBAmodel)
      
      # set obj function (here maximize mesophyll and bundle sheath biomass)
      obj_coef(C4FBAmodel.const.NADPME)[match(c("BIO_m","BIO_b"),react_id(C4FBAmodel.const.NADPME))]<-0.055
      
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
      
      
      
      # conduct coupled FBA (save and load file, variability analysis takes a long time)
      
      FVA=TRUE
      
      #       #plasmodesmata reactions:
      #plasmoPos<-grep("TMB",react_id(C4FBAmodel.const.NADPME))
      # create vector of weights for mtfFBA
      #mtf.weights<-rep(1,length=length(react_id(C4FBAmodel.const.NADPME)))
      # set a higher weight to plasmodesmatal reactions
      #mtf.weights[plasmoPos]<-1.1
      
      
      plasmo.Weights<-weight_vec
      raw.res<-vector(length(plasmo.Weights), mode="list") # list containing all results
      names(raw.res)<-plasmo.Weights
      plasmo.res<-vector(length(plasmo.Weights), mode="list") # list containing all plasmodesmata fluxes
      names(plasmo.res)<-plasmo.Weights
      C4.res<-vector(length(plasmo.Weights), mode="list") # list containing all C4 relevant fluxes
      names(C4.res)<-plasmo.Weights
      
      for (i in seq_along(plasmo.Weights)){
        # create vector of weights for mtfFBA
        mtf.weights<-rep(1,length=length(react_id(C4FBAmodel.const.NADPME)))
        # set a higher weight to plasmodesmatal reactions
        mtf.weights[plasmoPos]<-plasmo.Weights[i]
        
        vCFBAres<-vC.C3C4.FBA(pars=parms,
                              Cm=C,
                              Om=O,
                              Etot=Etot,
                              mechModFile="Models/vC C3C4 Cm Landscape model Rubisco pars kcat dependent FAST.r",
                              FBAmodelorg=C4FBAmodel.const.NADPME,
                              alg="mtf",
                              varFBA=FVA,
                              #costcoefbw=mtf.weights # costcoefbw copies to costcoeffw
                              costcoefMTF=mtf.weights#,
                              #varReactions=plasmoPos
        )
        
        #save(vCFBAres,file="vCFBAres_Multiple_AAs")
        #load("vCFBAres")
        
        if(FVA){
          # create matrix with cols react_name, minflux, maxflux
          varMinMax<-cbind(react_name(C4FBAmodel.const.NADPME),apply(getFluxDist(vCFBAres$varFBA),1,min),apply(getFluxDist(vCFBAres$varFBA),1,max) ) 
        }else{
          varMinMax<-cbind(react_name(C4FBAmodel.const.NADPME),rep(NA,length(react_name(C4FBAmodel.const.NADPME))),rep(NA,length(react_name(C4FBAmodel.const.NADPME)))) 
        }
        
        #check solver status
        if(!all(vCFBAres$varFBA@lp_ok==0)){stop("!all(vCFBAres$varFBA@lp_ok==0)")}
        if(!all(vCFBAres$varFBA@lp_stat==5)){stop("!all(vCFBAres$varFBA@lp_stat==5)")}
        
        raw.res[[i]]<-vCFBAres
        
        
        # evaluate results:
        vCFBAres.flux <- getFluxDist(vCFBAres$optSolObj)
        
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
        
        plasmo.res[[i]]<-plasmoFluxes
        
        
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
        C4.res[[i]]<-C4fluxes
        
      }
      return(list("plasmo.Weights"=plasmo.Weights,"raw.res"=raw.res,"plasmo.res"=plasmo.res,"C4.res"=C4.res,"modelOrgObj"=C4FBAmodel.const.NADPME))
}

#save(plasmo.Weights,raw.res,plasmo.res,C4.res,file="PlasmoWeightSensResult_Multiple_AAs")
