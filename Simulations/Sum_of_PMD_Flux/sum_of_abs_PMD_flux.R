# 250414
# calculate sum of plasmodesmatal flux, use data caclculated for Sens Analysis

setwd("Simulations/MFA_Weight_Sens/results/")



for(file in c("PlasmoWeightSensResult_Multiple_AAs_1to5",
              "PlasmoWeightSensResult_wo_aKG_Glu_1to5",
              "PlasmoWeightSensResult_wo_aKG_Glu_Ala_Pyr_1to5",
              "PlasmoWeightSensResult_only_Ser_Gly_1to5",
              "PlasmoWeightSensResult_Multiple_AAs_C4Cycle_active_1to5"              
)){
  load(file)
  plasmo<-sens_res[["plasmo.res"]][["1.1"]]
  print(file)
  print(sum(abs(as.numeric(plasmo[,"flux"]))))
  rm(sens_res)# safety

}