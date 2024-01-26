# 250414
# write data on sensitivity to xls files


list.of.Mats.to.xls<-function(LoMats,...){
  # Sens returned list of matrices, convert to dfs:
  plasmo.dfs<-lapply(LoMats,as.data.frame,stringsAsFactors=FALSE)
  # round numeic values
  plasmo.dfs.rounded<-lapply(plasmo.dfs,function(x){x$flux=round(as.numeric(x$flux),3);
                                                    x$min=round(as.numeric(x$min),3);
                                                    x$max=round(as.numeric(x$max),3);
                                                    return(x)}
                             )
  
  x<-plasmo.dfs.rounded
  #return(x)
  library(WriteXLS)
  WriteXLS(...,x="x",SheetNames = names(LoMats))
}

for(file in c("PlasmoWeightSensResult_Multiple_AAs_1to5",
              "PlasmoWeightSensResult_wo_aKG_Glu_1to5",
              "PlasmoWeightSensResult_wo_aKG_Glu_Ala_Pyr_1to5",
              "PlasmoWeightSensResult_only_Ser_Gly_1to5",
              "PlasmoWeightSensResult_Multiple_AAs_C4Cycle_active_1to5"              
              )){
      load(file)
      plasmo<-sens_res[["plasmo.res"]]
      rm(sens_res)# safety
      list.of.Mats.to.xls(plasmo,ExcelFileName=paste0(file,".xlsx"))
}