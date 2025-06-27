### CODE FOR CREATING THE VIRUS PANELS 
### INDICATE UPDATES HERE: 

### PACKAGES 
library(readxl)
library(dplyr)
library(tidyr)
library(clevr)
library(foreach)
library(doParallel)
library(doSNOW)

### DIRECTORY 
setwd("...")

# LOAD THE NEUTRALIZATION DATA FOR KNOWN mAb
# NEEDS TO BE IN FOLLOWING FORMAT FOR RUNNING THE create_panels FUNCTION: 
# FIRST COLUMN: "mAb"
# SECOND COLUMN: "Epitope"
# REST OF THE COLUMN = VIRUS PANEL MEASUREMENTS (ONE COLUMN = ONE VIRUS)
# ENSURE THAT YOU REMOVE CNE40_BC

## mAb
data.mAb.panel<- read_xlsx("data_mAb_panel.xlsx")


# FUNCTION FOR CREATING THE SUB VIRUS PANELS
# ARGUMENTS OF THE FUNCTION:
# neut data
# min.panel: minimum size of the virus panels (default=10)
# max.panel: maximum size of the virus panels (default=35)
# n.panel: number panels tested of each size (default=2500)
create_panels <- function(data.mAb.panel,min.panel,max.panel,n.panel) {
  
  data.forcor<-as.data.frame(t(as.matrix((data.mAb.panel)[,-c(1,2)])))
  colnames(data.forcor)<-data.mAb.panel$mAb  
  data.panel<-data.forcor
  
  ## UPDATE HERE IF bnAb CATEGORIES ARE CHANGED 
  mAb_existingcl<-data.mAb.panel %>%
    dplyr::select(mAb,Epitope) %>% 
    mutate(Epitope=dplyr::recode(Epitope,`CD4bs`="1",`IF_FP`="2",`MPER`="3",
                          `SF`="4",`V1V2`="5",`V3-Glycan`="6",`V3-Glycan II`="7"))
  
  num.panel<-max.panel-min.panel+1
  list_panel <- vector(mode = "list", length = num.panel)
  
  set.seed(98)
  seq.length<-ncol(data.mAb.panel)-2
  x <- seq(1, seq.length)
  
  M <- n.panel
  for (j in 1:num.panel){
    N <- min.panel-1+j
    out <- list() 
    while(length(out) < M) {
      out <- c(out,
               unique(replicate(M - length(out), sort(sample(x, N)), simplify = FALSE)))
    }
    list_panel[[j]]<-out
  }
  
  total_length<-num.panel*M
  
  parallel::detectCores()
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  
  NCLUSTER<-length(unique(mAb_existingcl$Epitope))
  
  homogen_vector<-foreach (
    j = 1:total_length,
    .combine = 'rbind',
    .packages='clevr'
  ) %dopar% {
    nlist.temp<-floor((j-1)/M)+1
    klist.temp<-j-(floor((j-1)/M)*M)
    combvirus.temp<-list_panel[[nlist.temp]][klist.temp]
    data.red<-data.panel[unlist(combvirus.temp),]
    cor.matrix<-cor(data.red,method="spearman",use="pairwise.complete.obs")
    if (sum(is.na(cor.matrix))=="0"){
      ncluster<-NCLUSTER
      hier.cl<-hclust(as.dist(1-cor.matrix))
      cl.attr<-as.data.frame(cutree(hier.cl,ncluster))
      colnames(cl.attr)<-"cluster"
      cl.attr$mAb<-rownames(cl.attr)
      
      ## COMPUTE HOMOGENEITY OF CLUSTERS
      data.homogen<-merge(cl.attr,
                          mAb_existingcl,by="mAb")
      c(j,homogeneity(data.homogen$Epitope,data.homogen$cluster),
        completeness(data.homogen$Epitope,data.homogen$cluster),
        v_measure(data.homogen$Epitope,data.homogen$cluster),
        rand_index(data.homogen$Epitope,data.homogen$cluster))
    } 
  } 
  
  homogen_data<-as.data.frame(homogen_vector)
  
  homogen_order<-homogen_data %>%
    arrange(desc(V2),desc(V3),desc(V4))
  
  Npanel<-40
  panel.tokeep<-homogen_order[1:Npanel,]
  
  Npanel<-nrow(panel.tokeep)
  panel_opt <- vector(mode = "list", length = Npanel)
  
  for (i in 1:Npanel){
    j<-panel.tokeep$V1[i]
    nlist.temp<-floor((j-1)/M)+1
    klist.temp<-j-(floor((j-1)/M)*M)
    combvirus.temp<-list_panel[[nlist.temp]][klist.temp]
    panel_opt[[i]]<-rownames(data.panel)[unlist(combvirus.temp)]
  }
  
  ### OUTCOME
  list(homogen_data,panel.tokeep,panel_opt)
}

# RUN THE FUNCTION 
outcome<-create_panels(data.mAb.panel,min.panel=10,max.panel=35,n.panel=2500)

# INDICATE HERE WHERE YOU WANT TO SAVE THE OUTCOME OF THE FUNCTION
homogen_data<-outcome[[1]]
write.table(homogen_data,"allpanels_date.txt",row.names = FALSE)

panel.tokeep<-outcome[[2]]
write.table(panel.tokeep,"panelopt_date.txt",row.names = FALSE)

panel_opt<-outcome[[3]]
save(panel_opt, file="panelopt_date.RData")

