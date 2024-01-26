setwd("Simulations/MFA_Weight_Sens/")

source("Simulations/MFA_Weight_Sens/vC_FBA_weight_Sens_Function.R")

sens_res<-vC_FBA_weight_Sens(
  constraint_file="Constraints/reactionConstraintsBMobj_NADPMEC2-C4_Multiple_AAs_only_Ser_Gly.tsv",
  weight_vec=c(1, 1.1, 2, 3, 4, 5)#c(1, 1.1, 5, 10, 50, 100)
)

save(sens_res,file="PlasmoWeightSensResult_only_Ser_Gly_1to5")
