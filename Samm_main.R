# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014") 
library("ggplot2")
library("deSolve")


climate<-read.csv("climate_file.csv", sep=";")

################################################################################
# BEGIN SAMM FUNCTION

# Functions that are needed in the model

RateModTempExp <- function(Temp){
  rateMod <- 2^((Temp-20)/10)
  return(rateMod)
}

RateModWater <- function(pF){
  corr_M = 1
  if(pF>=2.5){corr_M = (1.625-pF/4) }
  if(pF>=6.5){corr_M = 0}
  if(pF<=1.5){corr_M=(0.6+0.4 * pF/1.5)}
  if(pF<=0){corr_M=0.6}
  return(corr_M)
}

CUE_CN_LMW <- function(LMW_C,LMW_N,CUE_LMW){
  CUE_Sinsabaugh_linear= (((1E-6+LMW_C)/(LMW_N+1E-10))^(-1)) * 13.4 # Formula  extracted from Campbell 2016
  if(is.nan(CUE_Sinsabaugh_linear)){CUE_Sinsabaugh_linear<-0}
  if(is.na(CUE_Sinsabaugh_linear)){CUE_Sinsabaugh_linear<-0}
  if(CUE_Sinsabaugh_linear>1){CUE_Sinsabaugh_linear=1}
  return(CUE_LMW*CUE_Sinsabaugh_linear)
}




# This is SAMM model v1.0

SAMMv1 <- function(Time, State, Pars, Weather){
  with(as.list(c(State, Pars)),{
    
    ############################################################################
    # Set all external driving variables needed to calculate the flows
    ############################################################################
    
    # For the calculation of leaching
    WC_theta = Weather$Simulated_WC_theta[floor(Time)]
    H2O_content=(WC_theta * simulated_depth * 10000) * 1000 #H2O content in kg/ha
    simulated_leaching_m=Weather$Simulated.leaching..mm.[floor(Time)]	/1000	 # leaching converted from mm to m
    H2O_leached=(simulated_leaching_m) * 10000 * 1000 #H2O leached in kg/ha
    
    # Soil temperature and moisture
    temperature=Weather$mean_temp_7.5cm[floor(Time)]
    pF=Weather$pF_from_theta[floor(Time)]
    
    # Calculate the rate modifiers
    corr_F =RateModWater(pF) * RateModTempExp(temperature)
    
    # Calculate the helper variables
    p_LAB = STR_C * pc_STR_LAB/LAB_C # New version - 12.12.2020
    if(p_LAB>1){p_LAB=1	}
    
    a_MIC = max(MIC_C/(K_M_MIC_C+MIC_C),0.05) 
    
    
    
    
    ############################################################################
    # Computing the Carbon flows between pools
    ############################################################################    
    
    STR_C_to_CO2 = STR_C * k_STR * a_MIC * (1-CUE_STR) * corr_F 
    
    STR_C_to_LMW_C = STR_C * k_STR * a_MIC * (CUE_STR) * corr_F
    
    LAB_C_to_CO2 = LAB_C * k_LAB * (1-p_LAB) * a_MIC * (1-CUE_LAB) * corr_F 
    
    LAB_C_to_LMW_C = LAB_C * k_LAB * a_MIC * (1-p_LAB) * (CUE_LAB) * corr_F 
    
    LMW_C_to_CO2 = mu_max * LMW_C * a_MIC * (1-CUE_CN_LMW(LMW_C,LMW_N,CUE_LMW)) * corr_F 
    
    LMW_C_to_MIC_C = mu_max * LMW_C * a_MIC * CUE_CN_LMW(LMW_C,LMW_N,CUE_LMW) * corr_F
    
    STR_C_to_AggSTR_C = min((LMW_C_to_MIC_C+NonMicAgg) * aggfactSTR_C,STR_C) 
    
    LAB_C_to_AggLAB_C = min(STR_C_to_AggSTR_C * pc_STR_LAB,LAB_C)
    
    MAO_C_to_AggMAO_C = min((LMW_C_to_MIC_C+NonMicAgg) * aggfactMAO_C,MAO_C) 
    
    AggSTR_C_to_STR_C = AggSTR_C * kAgg * corr_F
    
    AggLAB_C_to_LAB_C = AggLAB_C * kAgg * corr_F
    
    AggMAO_C_to_MAO_C = AggMAO_C  * kAgg * corr_F
    
    MIC_C_to_CO2 = MIC_C * m_MIC * corr_F
    
    MIC_C_to_MAO_C= f_MICMAOM * MIC_C * k_MIC * corr_F
    
    MIC_C_to_LMW_C= (1-f_MICMAOM) * MIC_C * k_MIC * corr_F
    
    MAO_C_to_LMW_C= MAO_C * k_MAO * a_MIC * corr_F  
    
    
    # leaching from external driver (HYDRUS 1D)
    H2O_LMW_C_concentration=(LMW_C)/H2O_content
    
    LMW_C_leaching= H2O_leached * H2O_LMW_C_concentration
    
    if (LMW_C_leaching>(LMW_C * 0.95)){
      LMW_C_leaching= LMW_C * 0.95	#not more than 95% can be leached
    }
    
    LMW_C_to_MAO_C = LMW_C * K_LMW_MAO * (MAO_C_max-MAO_C)/MAO_C_max * corr_F 
    #MAO_C_to_LMW_CDesorption = MAO_C * KMaoc_to_DocDesorption * (MAO_C/MAO_C_max) * corr_F 
    MAO_C_to_LMW_CDesorption = 0 # Deactivated desorption for now, to reduce equifinality and it is essentially the same as MAO_C_to_LMW_C
    
    
    
    ############################################################################
    # Computing changes in pool sizes dP/dt for Carbon
    ############################################################################  
    
    dSTR_C = -STR_C_to_LMW_C -STR_C_to_CO2 -STR_C_to_AggSTR_C +AggSTR_C_to_STR_C
    
    dLAB_C = -LAB_C_to_LMW_C -LAB_C_to_CO2 -LAB_C_to_AggLAB_C +AggLAB_C_to_LAB_C
    
    dLMW_C = +STR_C_to_LMW_C +LAB_C_to_LMW_C +MIC_C_to_LMW_C +MAO_C_to_LMW_C +
      -LMW_C_to_CO2 -LMW_C_to_MIC_C -LMW_C_to_MAO_C +MAO_C_to_LMW_CDesorption -LMW_C_leaching
    
    dMIC_C = +LMW_C_to_MIC_C -MIC_C_to_CO2 -MIC_C_to_MAO_C -MIC_C_to_LMW_C
    
    dMAO_C = +MIC_C_to_MAO_C +LMW_C_to_MAO_C -MAO_C_to_LMW_CDesorption -MAO_C_to_LMW_C -MAO_C_to_AggMAO_C +AggMAO_C_to_MAO_C
    
    dAggSTR_C = +STR_C_to_AggSTR_C -AggSTR_C_to_STR_C
    
    dAggLAB_C = +LAB_C_to_AggLAB_C -AggLAB_C_to_LAB_C
    
    dAggMAO_C = +MAO_C_to_AggMAO_C -AggMAO_C_to_MAO_C
    
    #adding CO2 for mass balance equation
    dCO2 = +MIC_C_to_CO2 +LMW_C_to_CO2 +LAB_C_to_CO2 +STR_C_to_CO2
    
    dLeachC = +LMW_C_leaching
    
    
    ############################################################################
    # Computing the Nitrogen flows between pools
    ############################################################################
    
    LAB_N_to_LMW_N = LAB_N * k_LAB * a_MIC * (1-p_LAB) * corr_F		
    if (LAB_N_to_LMW_N<0){LAB_N_to_LMW_N= 0}
    
    LMW_N_to_Mic_N = mu_max * LMW_N * a_MIC * corr_F		
    
    MIC_N_to_MAO_N =  f_MICMAOM * MIC_N * k_MIC * corr_F 
    
    MIC_N_to_LMW_N = (1-f_MICMAOM) * MIC_N * k_MIC * corr_F 
    
    MAO_N_to_LMW_N = MAO_N * k_MAO * a_MIC * corr_F 
    
    LMW_N_to_MAO_N = LMW_C_to_MAO_C * (LMW_N/LMW_C)
    
    MAO_to_LMW_NDesorption = MAO_C_to_LMW_CDesorption * (MAO_N/MAO_C) # Desorped at same rate, atm 0
    
    LAB_N_to_AggLAB_N = LAB_C_to_AggLAB_C * (LAB_N/LAB_C)
    
    MAO_N_to_AggMAO_N = MAO_C_to_AggMAO_C * (MAO_N/MAO_C)
    
    AggLAB_N_to_LAB_N = AggLAB_N * kAgg * corr_F
    
    AggMAO_N_to_MAO_N = AggMAO_N * kAgg * corr_F
    
    
    
    H2O_LMW_N_concentration=(LMW_N)/H2O_content
    LMW_N_leaching= H2O_leached * H2O_LMW_N_concentration
    if (LMW_N_leaching>(LMW_N * 0.95)){
      LMW_N_leaching= LMW_N * 0.95	#not more than 95% can be leached
    }
    
    
    # Calculate LMW_N release and immobilization    
    MicLMW_NRelease = 0
    MIC_CN= (MIC_C+1E-10)/(MIC_N+1E-11)
    if (MIC_CN<MicMinCN){
      MicLMW_NRelease=0.5 * (MIC_N-(MIC_C/MicMinCN)) 
    }
    
    MicLMW_NImmobilization = 0
    max_MicLMW_NImmobilization = LMW_N * 0.5
    LMW_N_need = (MIC_C/MicMaxCN)-MIC_N
    if (LMW_N_need>0){
      MicLMW_NImmobilization= max_MicLMW_NImmobilization      
    }
    if (max_MicLMW_NImmobilization>LMW_N_need){
      MicLMW_NImmobilization = LMW_N_need      
    }
    
    
    
    
    ############################################################################
    # Computing changes in pool sizes dP/dt for Nitrogen
    ############################################################################
    dLAB_N = -LAB_N_to_LMW_N -LAB_N_to_AggLAB_N +AggLAB_N_to_LAB_N
    dLMW_N = +LAB_N_to_LMW_N -LMW_N_to_Mic_N +MIC_N_to_LMW_N +MAO_N_to_LMW_N +MicLMW_NRelease -MicLMW_NImmobilization -LMW_N_to_MAO_N +
      + MAO_to_LMW_NDesorption -LMW_N_leaching
    dMIC_N = +LMW_N_to_Mic_N -MIC_N_to_MAO_N -MIC_N_to_LMW_N -MicLMW_NRelease +MicLMW_NImmobilization
    dMAO_N =  +MIC_N_to_MAO_N +LMW_N_to_MAO_N -MAO_to_LMW_NDesorption -MAO_N_to_LMW_N -MAO_N_to_AggMAO_N +AggMAO_N_to_MAO_N
    
    dAggLAB_N = LAB_N_to_AggLAB_N -AggLAB_N_to_LAB_N
    dAggMAO_N = MAO_N_to_AggMAO_N -AggMAO_N_to_MAO_N
    
    
    
    ############################################################################
    # The change of Litterbag pools for calibration
    ############################################################################
    dLitterBagSTR_C= -LitterBagSTR_C * k_STR * a_MIC * corr_F 
    
    # Calculate the helper variables
    LitterBagSTR_CLAB_CProtection = LitterBagSTR_C * pc_STR_LAB/LitterBagLAB_C # New version - 12.12.2020
    if(LitterBagSTR_CLAB_CProtection>1){LitterBagSTR_CLAB_CProtection=1	}
    
    dLitterBagLAB_C= -LitterBagLAB_C * (1-LitterBagSTR_CLAB_CProtection) * k_LAB * a_MIC * corr_F 
    dLitterBagLAB_N= -LitterBagLAB_N * (1-LitterBagSTR_CLAB_CProtection) * k_LAB * a_MIC * corr_F 
    
    
    
    ############################################################################
    # The application of organic resources
    ############################################################################   
    
    #empty Litterbag pool at end of year to make room for new 
    if(floor(Time) %in% (c(382,759,1125,1490,1874,2231,2596,2961,3323,3688,4049,4414,4779,5144,5509,5874,6240,6605,6970,7335,7702,8068,8432,8797,9162)-13)){
      dLitterBagSTR_C= - LitterBagSTR_C * 0.99999
      dLitterBagLAB_C= - LitterBagLAB_C * 0.99999
      dLitterBagLAB_N= - LitterBagLAB_N * 0.99999
    }
    
    # Add litter at the right dates
    if(floor(Time) %in% c(1,382,759,1125,1490,1874,2231,2596,2961,3323,3688,4049,4414,4779,5144,5509,5874,6240,6605,6970,7335,7702,8068,8432,8797,9162)){
      
      dSTR_C = dSTR_C+LitterC * LitterStructuralPercent
      dLAB_C = dLAB_C+ LitterC * (1-LitterStructuralPercent)
      dLAB_N =dLAB_N+LitterN
      
      dLitterBagSTR_C=dLitterBagSTR_C+LitterC * LitterStructuralPercent
      dLitterBagLAB_C=dLitterBagLAB_C+LitterC * (1-LitterStructuralPercent)
      dLitterBagLAB_N=dLitterBagLAB_N+LitterN
      
    }else{
      
      ############################################################################
      # The mass balance is checked on DaysWhere not residue input occur
      ############################################################################     
      
      if(abs(dSTR_C+dLAB_C+dLMW_C+dMIC_C+dMAO_C+dAggSTR_C+dAggLAB_C+dAggMAO_C+dCO2+dLeachC)>0.001){
        # If the mass balance is not closed, an error message is produced
        print(floor(Time))
        stop("Error: the mass balance C is not closed")
      }
      if(abs(dLAB_N+dLMW_N+dMIC_N+dMAO_N+dAggLAB_N+dAggMAO_N+LMW_N_leaching)>0.001){
        # If the mass balance is not closed, an error message is produced
        print(floor(Time))
        stop("Error: the mass balance N is not closed")
      }
      
    }
    
    
    
    ############################################################################
    # The addition of daily C is done
    ############################################################################   
    
    if ((floor(Time)  > 1)){
      dSTR_C = dSTR_C + Daily_Litter_C * Daily_LitterSTR_C
      dLAB_C = dLAB_C + Daily_Litter_C * (1-Daily_LitterSTR_C)
      dLAB_N = dLAB_N + (Daily_Litter_C/Daily_Litter_CN)
    }
    
    
    return(list(c(dSTR_C,
                  dLAB_C,
                  dLMW_C,
                  dMIC_C,
                  dMAO_C,
                  dAggSTR_C,
                  dAggLAB_C,
                  dAggMAO_C,
                  dLAB_N,
                  dLMW_N,
                  dMIC_N,
                  dMAO_N,
                  dAggLAB_N,
                  dAggMAO_N,
                  dLitterBagSTR_C,
                  dLitterBagLAB_C,
                  dLitterBagLAB_N
    )))
  })
}



################################################################################
# END SAMM FUNCTION













# --------------------------
# The parameters are defined
# --------------------------
parameters <- c(
  
  # [Turnover Rates]
  k_STR	= 0.0031,	  # Depolimerisation of Structural litter (d-1 kgC-1)
  k_LAB	= 0.0366,	    # Depolimerisation of Metabolic litter (d-1 kgC-1)
  k_MIC	= 0.007,	          # Microbial death respiration	(d-1)
  k_MAO =	0.0041,	            # MAOM decompositions rate (d-1)
  mu_max	= 0.337,        # Maximum daily uptake by microbes	(d-1 kgC-1) from WANG 2013 MEND MODEL
  kAgg =	0.0323,               #	kAgg decomposition rate (d-1)
  K_M_MIC_C	= 66.8,               #	MIC halt saturation constent (kgC ha-1) estimated half saturation constant
  m_MIC =	0.00064,    # Microbial maintainance respiration	(d-1)
  K_LMW_MAO = 0.018,
  c_SORP = 0.049, # new value based on Georgiou et al. (2025) for 1:1 clay minerals
  
  
  # [CUEs]		
  CUE_STR =	0.65,	
  CUE_LAB	= 0.73,	
  CUE_LMW	= 0.6,	
  
  #[other parameters]		
  MicMinCN = 5.72,
  MicMaxCN = 7.77,
  f_MICMAOM	= 0.21,	
  #BiAC_release_per_growth	= 0.01,	
  pc_STR_LAB	= 2.54,	 # (g g-1)
  aggfactSTR_C = 1.16,                                # (g C g-1)
  aggfactMAO_C =	1.57,	                             # (g C g-1)
  NonMicAgg = 32.5, # in equivalent of soil mic growth in kg C /ha
  #PolyphenolProtectionCapacity= 1,	                   # (g N/g C)
  Daily_Litter_C =	2.28,	                               # kgC ha-1 d-1
  Daily_Litter_CN	= 142.9,	                               # kgC ha-1 d-1
  Daily_LitterSTR_C	= 0.16,	                 # kgC ha-1 d-1

  #KMaoc_to_DocDesorption = 0.001,

  
  #[Added litter values]	# Groundnut	
  LitterC	= 3880,	#kg C/ha
  LitterN	= 228,	#kg N/ha
  LitterStructuralPercent	= 0.13,	#g/g
  PercentPolyphenolsStructural = 0.16,#	g/g
  
  
  #[General values]		
  LitterInAggregateStart = 1400,	# 	kg C/ha 
  InitialSystemC = 4900,		# kg C/ha
  si_cl = 0.176,		# g/g
  
  bulk_density = 1.45,	#t/m^3
  simulated_depth	= 0.15,	#m
  MAO_C_max = NA
)



# ---------------------------------------
# The initial size of the state variables
# ---------------------------------------

yini = c(
  STR_C =	NA,	#kg C/ha
  LAB_C = NA,	#kg C/ha
  LMW_C = 100,	#"kg C/ha	"
  MIC_C = 100,		# kg C/ha
  MAO_C = 1665-575,	# 	kg C/ha
  
  AggSTR_C =	1400 * 0.1 ,	#kg C/ha
  AggLAB_C = 1400 * 0.9 ,	#kg C/ha
  AggMAO_C = 575,	# 	kg C/ha
  
  #[Start values N]		
  LAB_N = NA,#	kg N/ha
  LMW_N	= 5,	#kg N/ha
  MIC_N	= 10,	#kg N/ha
  MAO_N =	NA,	#kg N/ha
  AggLAB_N =	NA ,	#kg C/ha
  AggMAO_N = NA,	# 	kg C/ha
  
  
  # # [LitterBAGValues subtract]	
  LitterBagSTR_C = 0.1,
  LitterBagLAB_C = 0.1,
  LitterBagLAB_N = 0.01
  #LitterBagSTR_Cy2 = 0.1,
  #LitterBagLAB_Cy2 = 0.1,
  #LitterBagLAB_Ny2 = 0.01
  
)

# Calculate start values for pools where it depends on parameter values
InitialLitterC = parameters["InitialSystemC"]-yini["LMW_C"]-yini["MIC_C"]-yini["MAO_C"]-yini["AggMAO_C"]-yini["AggSTR_C"]-yini["AggLAB_C"]
yini["STR_C"] = InitialLitterC  *  parameters["Daily_LitterSTR_C"]
yini["LAB_C"] = InitialLitterC - yini["STR_C"]
yini["LAB_N"] = yini["LAB_C"]/ parameters["Daily_Litter_CN"]
yini["MAO_N"] = yini["MAO_C"]/9 

yini["AggLAB_N"] = yini["AggLAB_C"]/ parameters["Daily_Litter_CN"]
yini["AggMAO_N"] = yini["AggMAO_C"]/9 

yini



# Corrected equation for MAO_C_max
parameters["MAO_C_max"]=as.double(parameters["simulated_depth"]) * (parameters["bulk_density"]*1000 * 10000)* as.double((parameters["si_cl"]))* as.double((parameters["c_SORP"]))




# --------------------------------------
# The times at which output is requested
# --------------------------------------
times <- seq(1,9070,1)



# Plot SOC only the day before residue application
TimesBeforeApplying <-(c(1,382,759,1125,1490,1874,2231,2596,2961,3323,3688,4049,4414,4779,5144,5509,5874,6240,6605,6970,7335,7702,8068,8432,8797,9162))



for(i in 1:1){

  parsGN<-parameters
  
  parsCT<-parameters
  parsCT["LitterC"] = 0
  parsCT["LitterN"] = 0
  parsCT["LitterStructuralPercent"] = 0
  parsCT["PercentPolyphenolsStructural"] = 0
  
  parsDP<-parameters
  parsDP["LitterC"] = 4530
  parsDP["LitterN"] = 57
  parsDP["LitterStructuralPercent"] = 0.32
  parsDP["PercentPolyphenolsStructural"] = 0.27
  
  parsRS<-parameters
  parsRS["LitterC"] = 3670
  parsRS["LitterN"] = 47
  parsRS["LitterStructuralPercent"] = 0.06
  parsRS["PercentPolyphenolsStructural"] = 0.18
  
  parsTM<-parameters
  parsTM["LitterC"] = 4300
  parsTM["LitterN"] = 136
  parsTM["LitterStructuralPercent"] = 0.17
  parsTM["PercentPolyphenolsStructural"] = 0.25
}


ST1<-Sys.time()
outCT <- as.data.frame(rk4(func = SAMMv1, y = yini, parms = parsCT, times = times, Weather = climate));outCT$TRT<-"CT"
outDP <- as.data.frame(rk4(func = SAMMv1, y = yini, parms = parsDP, times = times, Weather = climate));outDP$TRT<-"DP"
outGN <- as.data.frame(rk4(func = SAMMv1, y = yini, parms = parsGN, times = times, Weather = climate));outGN$TRT<-"GN"
outRS <- as.data.frame(rk4(func = SAMMv1, y = yini, parms = parsRS, times = times, Weather = climate));outRS$TRT<-"RS"
outTM <- as.data.frame(rk4(func = SAMMv1, y = yini, parms = parsTM, times = times, Weather = climate));outTM$TRT<-"TM"
outComb<-rbind(outCT,outDP,outGN,outRS,outTM)
print(Sys.time()-ST1)


outComb$AgC<- outComb$AggSTR_C+
  outComb$AggLAB_C+
  outComb$AggMAO_C

outComb$C_N<-(outComb$LMW_C+outComb$MIC_C+outComb$MAO_C+outComb$AggSTR_C+outComb$AggLAB_C+outComb$AggMAO_C+outComb$STR_C+outComb$LAB_C)/
  (outComb$LMW_N+outComb$MIC_N+outComb$MAO_N+outComb$AggLAB_N+outComb$AggMAO_N+outComb$LAB_N)


outComb$AggMAO_C

ggplot() +
  geom_line(data = outComb, aes(x = time/365, y = MIC_C+MAO_C+AggSTR_C+AggLAB_C+AggMAO_C, color = TRT), linetype="dashed") +
  labs(title="Simulated amount of carbon",
       y=expression(paste("Carbon (kg C ha"^{-1}~")")),
       x="Time (year)") +
  theme_classic() +
  coord_cartesian(ylim = c(0, NA))+
  facet_grid(TRT~.)+
  theme(legend.position="NONE", legend.text = element_text(size=10),
        legend.background = element_rect(colour = "black", linetype = "solid"))+
  geom_hline(yintercept=3700,linetype="dashed",size=0.3)+
  geom_hline(yintercept=6300,linetype="dashed",size=0.3)+
  scale_y_continuous(limits = c(0,12000),breaks=c(0,5000,10000))+
  scale_color_viridis_d()






