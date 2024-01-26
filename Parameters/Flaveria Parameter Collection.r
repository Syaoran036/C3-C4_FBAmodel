# Flaveria Photosynthetic parameters

#tab<-read.csv2("C:/Users/David/Documents/Bio/PhD/Experimental data/Flaveria Parameter Collection.csv")

# [g Chl / m^2] for Conversion:

# whole leaf PEPC-activity in leaf extracts [?mol / (mg Chl * h)]:
PEPCact<- c(NA,40,24,123,123,NA,NA,NA,207,NA,162,1096,553,594,756,600,945,791)

# converted whole leaf PEPC-activity in leaf extracts [?mol/ (sqm * s)]
# assuming 450 mg Chl / m^2
PEPCactSQM<-PEPCact*460/3600
#PEPC act measured with ali on 30/7/12  [?mol/ (sqm * s)]
PEPCactSQMown<-c(9.1,10.2,12.0,33.5,51,NA,9.7,11.2,NA,40.2,40.2,228.2,95.4,NA,542.3,239.2,262.3,NA)

PEPCactSQMall<-ifelse(is.na(PEPCactSQMown),PEPCactSQM,PEPCactSQMown)
#PEPCactSQM<-PEPCactSQMown
#PEPCact<-PEPCactSQMown

PEPCactAuthors <- c(NA,"Holaday A et al. 88", "Bauwe H. 84","Ku M et al. 83", "Ku M et al. 83", NA, NA,NA,"Ku M et al. 83",NA,"Ku M et al. 83","Moore et al. 89","Cheng et al 88","Moore et al. 89","Moore et al. 89","Moore et al. 89","Moore et al. 89","Sudderth et al 2007")
# Which Assay conditions were cited? :
PEPCactAssay <- c(NA,"Hatch 78 pH 8,3","?","Uedan 76 pH 8","Uedan 76 pH 8",NA,NA,NA,"Uedan 76 pH 8",NA,"Uedan 76 pH 8","Uedan 76 pH 8","Uedan 76 pH 8","Uedan 76 pH 8","Uedan 76 pH 8","Uedan 76 pH 8","Uedan 76 pH 8","Pittermann 2000")



SpecNames <- c("Flaveria robusta","Flaveria pringlei","Flaveria cronquistii","Flaveria anomala","Flaveria linearis","Flaveria sonorensis","Flaveria chloraefolia","Flaveria angustifolia",
               "Flaveria pubescens","Flaveria floridana","Flaveria ramosissima","Flaveria vaginata","Flaveria brownii","Flaveria palmeri","Flaveria australasica","Flaveria trinervia","Flaveria bidentis","Flaveria kochiana")
PStypes <- c("C3","C3","C3","C3-C4","C3-C4","C3-C4 type I","C3-C4 type I","C3-C4 type I","C3-C4 type I","C3-C4 type II","C3-C4 type II","C4-like","C4-like","C4-like","C4","C4","C4","C4")
# 13C composition Monson 88, Sudderth 07, Apel 1988
delta13C<-c(NA,-28.8, -29.4, -28.3, -27.9, NA, -27.8, NA, -28.3, -29.9, -28.5, -15.4, -17.4, -16.5,-15.0, -14.3,-15.1,-15.2)
# instantaneous PWUE Vogan & Sgae 2011
PWUE <- c(NA,4.37,4.81,NA,NA,4.11,4.84,4.64,4.62,4.09,4.49,7.13,7.48,9.01,NA,8.98,7.85,9.16)
# CO2 compensation point Gamma (at 30 ?C)  (Vogan & Sage 2011, ?) (originally Ku et al 1991?)
gamma <- c(62.1,62.0,60.4,15.5,27,29.6,29.0,24.1,21.3,9.5,9.0,3.0,6.0,4.7,5.1,3.5,3.2,2.2)
# %14C in C4 acids (Vogan & Sage 2011)
C4Acids<-c(NA,4.1 , 7.7 , NA, NA,3.0 , 11.3, 11.0, 24.9, 36.1, 43.0, 68.0, 68.5, 75.5, NA,     81.4, 72.0, NA)
# Kubien et al 2008
# kcat for CO2 [1/s] ([mol/mol] is Typo in figure 1)
kcat <- c(NA,3.11, 3.13, 3.8,  3.43, 2.69, 3.35, 2.86, NA, 3.19, 2.77, 3.78, 2.58, 3.54, 3.84, 4.42, 4.16, 3.68)
#kcat from Wessinger 1989:
kcatWes<-c(NA,3.8,3.9,NA,5.0,NA,NA,NA,3.8,3.4,NA,8.8,4.0,4.5,3.8,5.4,8.0,NA)
kcatAll<-ifelse(is.na(kcatWes),kcat,kcatWes)
# fraction of rubisco expressed in mesophyll
#from Bauwe 84 (visual observations), Moore 88 (calculation?), Moore 89, Cheng 88 (calculation?)
#beta <- c(NA,1,0.5,NA,NA,NA,NA,0.5,NA,0.52,0.06,0.35,0.068,0.054,0.002,0.008,NA)
#beta <- c(NA,1,0.5,NA,NA,NA,NA,0.5,NA,0.8,0.06,0.35,0.068,0.054,0.002,0.008,NA)
#25.4.12:Cheng 88 removed, Moore 88 calculated based on activity per protoplast:
#beta <- c(NA,NA,1,0.5,NA,NA,NA,NA,0.5,NA,0.65,0.06,NA,0.068,0.054,0.002,0.008,NA)
#5.6.12 corrected Bauwe 84 values for M:BS ratios.
beta <- c(NA,NA,0.95,0.78,NA,NA,NA,NA,0.75,NA,0.65,0.06,0.36,0.068,0.054,0.002,0.008,NA)
xi <-   c(NA,NA,NA,NA,0.97,NA,NA,NA,NA,0.98,NA,NA,NA,NA,NA,0.92,NA,NA)
xiTranscriptome <- c(0.59,0.50,NA,0.86,NA,NA,0.79,0.59,0.90,0.89,0.79,0.98,0.96,0.97,0.98,0.99,0.98,NA)
#new method to calculate xi:
xiTranscriptome2 <- c(0.18,0,NA,0.72,NA,NA,0.59,0.17,0.80,0.79,0.58,0.96,0.93,0.97,0.98,0.98,0.96,NA)

SpecNames <- c("Flaveria robusta","Flaveria pringlei","Flaveria cronquistii","Flaveria anomala","Flaveria linearis","Flaveria sonorensis","Flaveria chloraefolia","Flaveria angustifolia",
               "Flaveria pubescens","Flaveria floridana","Flaveria ramosissima","Flaveria vaginata","Flaveria brownii","Flaveria palmeri","Flaveria australasica","Flaveria trinervia","Flaveria bidentis","Flaveria kochiana")

# AlaAt RPM from Julia Mallmann Data
# Symbols: AlaAT1 | AlaAT1 (ALANINE AMINOTRANSFERAS); ATP binding / L-alanine:2-oxoglutarate aminotransferase | chr1:5922630-5926400 FORWARD
AlaAT.RPM<-c("Flaveria robusta"=202.731306,
             "Flaveria pringlei"=171.307942,
             "Flaveria cronquistii"=NA,
             "Flaveria anomala"=1620.178811,# Fa= anaomala?
             "Flaveria linearis"=NA,
             "Flaveria sonorensis"=NA,
             "Flaveria chloraefolia"=660.8779634,# Fc = chloraefolia
             "Flaveria angustifolia"=210.9280428,
             "Flaveria pubescens"=1068.65008,
             "Flaveria floridana"=NA,
             "Flaveria ramosissima"=1815.026426, 
             "Flaveria vaginata"=4169.832547,# FvE = vaginata
             "Flaveria brownii"=2694.266143,
             "Flaveria palmeri"=2982.725101,# FpaE_rrpm = palmeri?
             "Flaveria australasica"=2782.654058,
             "Flaveria trinervia"=2627.882901,
             "Flaveria bidentis"=2045.220069,
             "Flaveria kochiana"=NA)


#data from Moricandia:
mori.names<-c("M. arvensis")
mori.beta<-NA
mori.kcat<-NA
# Winter 82 (already in [?mol/sqm s])
mori.PEPCactSQM<-c(8.4)
mori.Kp<-NA
mori.gs<-NA
# immunolabeling data from hylton 1988, corrected for M:BS with data from Hattersley 1989
mori.xi<-c(0.97)


# data from Panicum
pani.names<-c("P. milioides","P. hians","P. miliaceum")
# milioides and hians: data from Ku 76, corrected for M:BS area from Wilson 83 and hattersley 84
# miliaceum: edwards 72, using chlM/chlBS fraction in the paper:
pani.beta<-c(0.76,0.83,0.16)
pani.kcat<-rep(NA,times=length(pani.names))
# data from Ku 76, in milloides and hians corrected with chl/sqm from Ku 78, in miliaceum from Usuda 84
pani.PEPCactSQM<-c(10.53,11.73,86.25)
pani.gs<-rep(NA,times=length(pani.names))
#immunolabeling from hylton 88, corrected for % leaf mitochondria in BS from hattersley 89
pani.xi<-c(0.97,NA,NA)
