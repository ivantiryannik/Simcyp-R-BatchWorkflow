# Simcyp-R Batch workflow

#  1.0 Setting up ----
# Load Required Packages 
library("Simcyp")
library("RSQLite")
library("Rcpp")

Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V21\\Screens\\SystemFiles",21,species = SpeciesID$Human, verbose = FALSE)

# Insert simcyp workspace 
path_user <-Simcyp::ScriptLocation()
setwd("-------")

#  2.0 Import Simcyp workspace ----
# import workspace
Simcyp::SetWorkspace("BIDhep0.26Vss31.25FU.wksz") 

# 3.0 Batch process ----
# Update parameters here
Ki<- "idIntRoutesKi11"
Ka<- "idka"
Vss<- "idEnteredVss"
Fg <- "idGutEMMean3A4"

# Define a vector of values   
newVss <- c(0.05,0.25,1.25,6.25,31.25)
# Fg changed via Population CYP3A4 Enzyme Abundance
newFg <- c(231.12, 58.25, 23.8, 9.05, 0)
newKa <- c(0.05,0.1,0.2,0.4,0.8,1,2,3,4,6)
newKi <- c(0.015,0.03,0.06,0.12,0.24,0.48,0.96,1.92,3.84,7.68)

# Run the simulations varying parameters
# For this we will be exporting our results to database. 
DBfilename <- "Batch.db"                               # Create the Database file name, with DB Extension
DBfilepath <- file.path(path_user, DBfilename)                       # File path to save the database results to

# For loop start ----
AUC_CMAX_CL_df<- NULL
start_time <- Sys.time()
for (i in 1:length(newFg)){
  FgValue <- newFg[i] 
  
  for (j in 1:length(newVss)){
    VssValue <- newVss[j] 
    
    for (k in 1:length(newKi)){
      KiValue <- newKi[k]
    
      for (l in 1:length(newKa)){
        KaValue <- newKa[l]
        
        # Set the Simcyp parameters : Fg, Vss, Ki, Ka
        SetPopulationParameter(Fg, PopulationID$GITract, as.numeric(FgValue))  # Fg
        SetCompoundParameter(Vss, CompoundID$Inhibitor1, as.numeric(VssValue))  # Vss
        SetCompoundParameter(Ki, CompoundID$Inhibitor1, as.numeric(KiValue))  # Ki
        SetCompoundParameter(Ka, CompoundID$Inhibitor1, as.numeric(KaValue))  # Ka
        
        capture.output(Simcyp::Simulate(database = DBfilepath), file='NUL')  # suppressing all outputs from Simcyp
        # Connect this R script to your Database file via RSQLite.
        conn <- RSQLite::dbConnect(SQLite(),DBfilepath)
        
        
        # 5.0 Extracting simulation result -----
        # PK summary results 
        AUC_deets_Sub<- GetAUCFrom_DB(ProfileID$CsysFull,CompoundID$Substrate,individual = 1,conn, allDoses = TRUE)
        # LastDose<-nrow(AUC_deets)/2
        
        # AUC inf ratio
        AUCinf<-  AUC_deets_Sub$AUC[AUC_deets_Sub$Inhibition==0][1]  #OVERALL Dose
        AUCinf_DDI<-  AUC_deets_Sub$AUC[AUC_deets_Sub$Inhibition==1][1]  #OVERALL Dose
        AUCinfRatio<- AUCinf_DDI/AUCinf
        
        # Cmax & AUC inhb
        AUC_deets_Inhb <- GetAUCFrom_DB(ProfileID$CsysFull,CompoundID$Inhibitor1,individual = 1,conn, allDoses = TRUE)
        Cmax_inhib <- AUC_deets_Inhb$Cmax[AUC_deets_Inhb$Inhibition==1][1]  #OVERALL Dose
        AUC_deets_Inhb2 <- GetAUCFrom_DB(ProfileID$CsysFull ,CompoundID$Inhibitor1,individual = 1,conn, allDoses = TRUE)
        AUCLastDose <- tail(AUC_deets_Inhb2$AUC, n=1)
        GetProfile_DB(ProfileID$CsysFull,CompoundID$Inhibitor1,individual = 1,conn)
        
        # Setting values for formula
        # Fg values
        if(FgValue == 231.12){
          FgSub <- 0.1
        }else{
          if(FgValue == 58.25){
            FgSub <- 0.3
          }else{
            if(FgValue == 23.8){
              FgSub <- 0.5
            }else{
              if(FgValue == 9.05){
                FgSub <- 0.7
              }else{
                if(FgValue == 0){
                  FgSub <- 1
                }
              }
            }
          }
        }
        
        # Other Values
        fu.p <- GetCompoundResult_DB(individual=1, ResultID$fuAdj, CompoundID$Inhibitor1, conn)   #Fraction unbound in Plasma for individual 1
        Fa <- 1    
        FgInh <- 1
        Ki1 <- KiValue
        Ka1 <- KaValue
        Dose <- GetCompoundParameter(CompoundParameterID$Dose,CompoundID$Inhibitor1)
        Qh <- 101.38
        Qen <- 23.85 
        Rb <- 0.62
        ClearanceInh <- Dose / AUCLastDose
        Css <- Dose / (ClearanceInh * 12)
        MW <- 531.4
        Css2 <- (Css/MW)*1000
        Cmax_inhib2 <- (Cmax_inhib/MW)*1000
        FmSub <- 0.9668
        fumic <- 0.97
        
        # AUCR MSM Formula - Cmaxcss
        Ih <- fu.p * (Cmax_inhib2 + (Fa * FgInh * Ka1 * Dose / MW * 1000) / Qh / Rb)
        Ig <- Fa * Ka1 * Dose / MW * 1000 / Qen
        Ah <- (1 / (1 + (Ih / (Ki1*fumic))))
        Ag <- (1 / (1 + (Ig / (Ki1*fumic))))
        AUCR <- (1 / (Ag * (1 - FgSub) + FgSub)) * (1 / (Ah * FmSub + (1 - FmSub)))
        IMDR <- AUCinfRatio / AUCR
        
        # AUCR MSM Formula - Cavrss
        Iha <- fu.p * (Css2 + (Fa * FgInh * Ka1 * Dose / MW * 1000) / Qh / Rb)
        Iga <- Fa * Ka1 * Dose / MW * 1000 / Qen
        Aha <- (1 / (1 + (Iha / (Ki1*fumic))))
        Aga <- (1 / (1 + (Iga / (Ki1*fumic))))
        AUCRA <- (1 / (Aga * (1 - FgSub) + FgSub)) * (1 / (Aha * FmSub + (1 - FmSub)))
        IMDRA <- AUCinfRatio / AUCRA
        
        
        #Results table
        AUC_CMAX_CL_df<- rbind(AUC_CMAX_CL_df,data.frame(Fg= FgValue,Vss= VssValue,Ki= KiValue, Ka= KaValue, AUCinfRatio, Cmax_inhib, Css, Cmax_inhib2, Css2, AUCLastDose, Iga, AUCR, AUCRA, IMDR, IMDRA, ClearanceInh))
        
        message(sprintf("Done Simcyp simulation; Fg = %5.3f,  Vss =  %5.3f, Ki =  %5.3f, Ka =  %5.3f  AUCinfRatio =  %5.3f,  Cmaxinhb =  %5.3f, Css =  %5.3f, IMDR =  %5.3f, IMDRA =  %5.3f, AUCLastDose =  %5.3f",FgValue,VssValue,KiValue,KaValue,AUCinfRatio,Cmax_inhib,Css,IMDR,IMDRA, AUCLastDose))
        
        
        ## End create matrix ------------------------------------------------
        RSQLite::dbDisconnect(conn)      #Detach database connection 
      
      }
      
    }
    
  }
  
}
End_time <- Sys.time()
BatchTime<- End_time-start_time
BatchTime

# Generating and Outputting matricies 
mydata <- list()

d <- 0# k is the index to access the AUC from the data frame created above

for(w in 1:25){
  
  favssmat <- matrix(data=NA, nrow=length(newKa), ncol=length(newKi))
  
  for(a in 1:length(newKi)){
    
    for(s in 1:length(newKa)){
      
      d <- d +1
      
      favssmat[s,a] <- AUC_CMAX_CL_df[d,]$AUCRA
      
    }
    
  }
  
  mydata[[w]] <- favssmat
  
}




for (i in seq_along(mydata)) {
  tryCatch(
    write.csv(mydata[[i]], file = paste0("matrix_", i, ".csv"), row.names = FALSE),
    error = function(e) print(paste("Error writing matrix", i, ":", e))
  )
}
