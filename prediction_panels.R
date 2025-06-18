### CODE FOR CREATING THE VIRUS PANELS PREDICTION
### INPUT data should have unique ID 
########

### PACKAGES 
library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)

### DIRECTORY 
setwd("...")


# LOAD THE SERUM NEUTRALIZATION DATA
# NEEDS TO BE IN FOLLOWING FORMAT FOR RUNNING THE create_panels FUNCTION: 
# FIRST COLUMN: uniqueID = patient ID (can be a combination of ID and DATE)
# REST OF THE COLUMNS = 1/NT50 measurements on virus panel (ENSURE THAT YOU INVERT THE MEASUREMENTS)
# ENSURE THAT YOU REMOVE CNE40_BC
# VIRUS NAMES SHOULD BE THE SAME IN data.patient AS IN data.mAb.panel

data.patient<-read_xlsx("data_patient.xlsx")

# LOAD OUTCOMES FROM THE create_virus_panels FUNCTION 
load("panelopt_date.Rdata")
paneltest<-read.table("panelopt_date.txt")

# LOAD mAb panel data (should be the same as input for the create_virus_panels FUNCTION)
data.mAb.panel<- read_xlsx("data_mAb_panel.xlsx")


# FUNCTION FOR RUNNING THE VIRUS PANELS METHOD
# ARGUMENTS OF THE FUNCTION:
# neut data from patients
# outcomes from create_virus_panels function
# neut data from mAb
# key between SHCS ID and AT lab ID
prediction_panels<-function(data.patient,panel_opt,paneltest,data.mAb.panel){
  
  colnames(data.patient)[1]<-"uniqueID"
  
  # ID MATCHING FOR PLOT 
  id.match<-as.data.frame(cbind(data.patient$`uniqueID`,seq(1,nrow(data.patient))))
  colnames(id.match)<-c("uniqueID","index")
  
  # PREDICTION AND CORRELATION MATRIX 
  Npan<-length(panel_opt)
  prediction_allpanel<-as.data.frame(matrix(NA,nrow=nrow(data.patient),ncol=(Npan+1)))
  prediction_allpanel$V1<-data.patient$`uniqueID`
  
  prediction_mAb<-as.data.frame(matrix(NA,nrow=nrow(data.patient),ncol=(Npan+1)))
  prediction_mAb$V1<-data.patient$`uniqueID`
  
  cluster.mAb<-data.mAb.panel %>% dplyr::select(mAb,Epitope)
  
  cor_matrix_all <- vector(mode = "list", length = nrow(data.patient))
  for (i in 1:nrow(data.patient)){
    cor_matrix_all[[i]]<-as.data.frame(matrix(NA,nrow=Npan,ncol=nrow(cluster.mAb)))
    colnames(cor_matrix_all[[i]])<-cluster.mAb$mAb
  }
  
  # COMPUTE CORRELATION FOR EACH PANEL 
  for (k in 1:Npan){
    # can remove the print(k) if you don't want to know where you're at in the code
    print(k)
    name.panel.temp<-panel_opt[[k]]
    panel.temp<-data.mAb.panel %>% dplyr::select(mAb,all_of(name.panel.temp))
    data.patient.temp<-data.patient %>% dplyr::select(`uniqueID`,all_of(name.panel.temp))
    cor.matrix<-as.data.frame(matrix(NA,nrow=nrow(data.patient.temp),ncol=nrow(panel.temp)+1))
    colnames(cor.matrix)<-c("uniqueID", panel.temp$mAb)
    cor.matrix$uniqueID<-data.patient$uniqueID
    for (i in 1:nrow(cor.matrix)){
      patient.temp<-cor.matrix$uniqueID[i]
      vec.pat.temp<-data.patient.temp %>% filter(uniqueID==patient.temp) %>% 
        dplyr::select(!c(uniqueID)) %>% as.vector() %>% t() %>% unlist()
      for (j in 1:nrow( panel.temp)){
        mab.temp<- panel.temp$mAb[j]
        vec.mab.temp<- panel.temp%>% filter(mAb==mab.temp) %>% 
          dplyr::select(!c(mAb)) %>% as.vector() %>% t() %>% unlist()
        cor.matrix[i,j+1]<-cor(vec.pat.temp,vec.mab.temp,method="spearman",use="pairwise.complete.obs")
      }
    }
    ## ALL CORRELATIONS
    for (i in 1:nrow(cor.matrix)){
      cor_matrix_all[[i]][k,]<-cor.matrix[i,-1]
    }
    
    ## MAX CORRELATION
    cor.matrix[is.na(cor.matrix)]<-(0-1)
    cor.matrix$maxcolnumber<-max.col(cor.matrix[,-1])
    cor.matrix$maxAb<-colnames(cor.matrix)[cor.matrix$maxcolnumber+1]
    nAb<-nrow(data.mAb.panel)
    for (i in 1:nrow(cor.matrix)){
      if (sum(cor.matrix[i,2:(nAb+1)])==-nAb){
        cor.matrix$maxcor[i]<-0-1
      }
      else{
        cor.matrix$maxcor[i]<-cor.matrix[i,cor.matrix$maxAb[i]]
      }
    }
    cor.matrix.cluster<-merge(cor.matrix,
                              cluster.mAb %>% dplyr::select(maxAb=mAb,Epitope),all.x=TRUE)
    
    cor.matrix.cluster$Epitope[cor.matrix.cluster$maxcor<0.40]<-"no prediction"
    cor.matrix.cluster$maxAb[cor.matrix.cluster$maxcor<0.40]<-"no prediction"
    cor.matrix.cluster<-cor.matrix.cluster %>%
      arrange(factor(uniqueID, levels = prediction_allpanel$V1))
    prediction_allpanel[,k+1]<-cor.matrix.cluster$Epitope
    prediction_mAb[,k+1]<-cor.matrix.cluster$maxAb
  } 
  
  # NEED TO UPDATE HERE IF OTHER EPITOPES 
  prediction_allpanel_group<-prediction_allpanel %>%
    mutate(CD4bsgroup=rowSums(prediction_allpanel == "CD4bs",na.rm = TRUE)) %>%
    mutate(IFFPgroup=rowSums(prediction_allpanel == "IF_FP",na.rm = TRUE)) %>%
    mutate(MPERgroup=rowSums(prediction_allpanel== "MPER",na.rm = TRUE)) %>%
    mutate(V1V2group=rowSums(prediction_allpanel == "V1V2",na.rm = TRUE)) %>%
    mutate(V3Glygroup=rowSums(prediction_allpanel== "V3-Glycan",na.rm = TRUE)) %>%
    mutate(V3GlyT2group=rowSums(prediction_allpanel== "V3-Glycan II",na.rm = TRUE)) %>%
    mutate(SFgroup=rowSums(prediction_allpanel== "SF",na.rm = TRUE)) %>%
    mutate(nopredictiongroup=rowSums(prediction_allpanel== "no prediction",na.rm = TRUE)) %>%
    mutate(sumall=rowSums(cbind(CD4bsgroup,
                                IFFPgroup,
                                MPERgroup,V1V2group,V3Glygroup,
                                V3GlyT2group,
                                nopredictiongroup,SFgroup))) %>%
    mutate(CD4bs_pct=CD4bsgroup/sumall,
           IFFP_pct=IFFPgroup/sumall,
           MPER_pct=MPERgroup/sumall,
           V1V2_pct=V1V2group/sumall,
           V3Gly_pct=V3Glygroup/sumall,
           V3GlyT2_pct=V3GlyT2group/sumall,
           SF_pct=SFgroup/sumall,
           nopred_pct=nopredictiongroup/sumall
    )
  
  data.final<-prediction_allpanel_group
  colnames(data.final)[1]<-"uniqueID"
  
  # CORRELATION AVERAGE VALUES
  data.average.cor<-as.data.frame(matrix(NA,nrow=nrow(data.patient),
                                         ncol=(nrow(cluster.mAb)+1)))
  colnames(data.average.cor)<-c("ID",cluster.mAb$mAb)
  data.average.cor$ID<-data.patient$uniqueID
  for (i in 1:nrow(data.average.cor)){
    for (j in 1:nrow(cluster.mAb)){
      data.average.cor[i,j+1]<-mean(cor_matrix_all[[i]][,j])
    }
  }
  
  # CORRELATION DISTRIBUTION PER PATIENT 
  data_long <- data.average.cor %>%                       
    pivot_longer(colnames(data.average.cor[,-1])) %>% 
    as.data.frame()
  
  cor.data<- vector(mode = "list", length = nrow(data.patient))
  
  for (i in 1:nrow(data.patient)){
    data.temp<-cor_matrix_all[[i]]
    data.temp$Panel<-seq(1:Npan)
    NCOL<-ncol(data.temp)
    
    ## COR DISTRIBUTION 
    data.temp.temp<-data.temp %>%
      pivot_longer(colnames(data.temp)[-NCOL],names_to = "mAb", values_to = "COR") 
    data.temp.temp<-merge(data.temp.temp,
                          cluster.mAb,by="mAb")
    test.temp<-data.temp.temp %>%
      group_by(mAb) %>%
      summarise(meanCOR=mean(COR)) %>%
      filter(meanCOR>0.4)
    data.temp.temp$Epitopecolor<-data.temp.temp$Epitope
    data.temp.temp$Epitopecolor[!data.temp.temp$mAb %in% test.temp$mAb]<-"under threshold"
    cor.data[[i]]<-data.temp.temp
  } 
  
  ### OUTCOME
  list(data.final,cor_matrix_all,data.average.cor,cor.data,id.match)
  
}

# RUN THE FUNCTION (~13 minutes on my mac)
outcome<-prediction_panels(data.patient,panel_opt,paneltest,data.mAb.panel) 

# INDICATE HERE WHERE YOU WANT TO SAVE THE OUTCOME OF THE FUNCTION
data.final<-outcome[[1]]
write.xlsx(data.final,"prediction_date.xlsx")
cor_matrix_all<-outcome[[2]]
save(cor_matrix_all, file="cormatrix_date.RData")
data.average.cor<-outcome[[3]]
write.xlsx(data.average.cor,"averagecor_date.xlsx")
cor.data<-outcome[[4]]
save(cor.data, file="cordata_date.RData")
id.match<-outcome[[5]]
write.xlsx(id.match, file="ID_index_date.xlsx")
