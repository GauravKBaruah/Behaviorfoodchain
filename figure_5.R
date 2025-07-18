
#######################################################################
#  Dynamical stability of the tri-trophic food chain (prey-predator abundance)
#### Figure 5

rm(list=ls())

# Install rEDM v0.6.9
remotes::install_version("rEDM", version = "0.6.9", repos = "https://cloud.r-project.org")

library(rEDM)
packageVersion("rEDM")  # should show ‘0.6.9’


# functions needed for the script

# function to do solve the weighted linear regression   
lm_svdsolve <- function(y, x, ws, subset = seq_along(y)){
  x <- x[subset,]
  y <- y[subset]
  ws <- ws[subset]
  # prepended column of 1s for constant term in linear model
  A <- cbind(1, x) * ws
  A_svd <- svd(A)
  # >>REMOVE SMALL SINGULAR VALUES<<
  s <- A_svd$d
  s_inv <- matrix(0, nrow = dim(x)[2]+1, ncol = dim(x)[2]+1)
  for(i in seq_along(s)){
    if(s[i] >= max(s) * 1e-5)
      s_inv[i,i] <- 1/s[i]
  }
  coeff <- s_inv %*% A_svd$v %*% t(A_svd$u) %*% (ws * y)
  coeff <- t(coeff)
  colnames(coeff) <- c("const",colnames(x))
  return(coeff)
}


#ranking based on EDM - empirical dynamical modelling i.e., S-map as its now called
#timerange : total time 
#tdata : the data subsetted
dom_eigenvalue_EDM<-function(timerange,tdat){
  tdat$ladybugs_all<-(tdat$Total_ladybirds)
  
  timerange <- timerange
  ts_data<-cbind(tdat$Aphids[1:timerange],tdat$ladybugs_all[1:timerange]) # appending two columns of prey and predator abundances
  colnames <- c()
  
  #renaming columns
  for(i in 1:1) {
    colnames <- c(colnames, paste("A", i, sep = ""))
  }
  
  for(i in 1:1) {
    colnames <- c(colnames, paste("P", i, sep = ""))
  }
  
  colnames(ts_data) <- colnames
  data_used<-1:timerange
  targ_colm<-1
  Embedding <- colnames   
  block <- ts_data[,Embedding]
  block <- as.data.frame(apply(block, 2, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
  # from the EDM Deyle et al Proc B paper, this is the theta used for localised linear wieghted regression
  
  test_theta <- block_lnlp(block[data_used,],
                           method = "s-map",
                           num_neighbors = 0,  
                           theta = c(0, 1e-04, 3e-04, 0.001, #we test many different theta and use the one with minimum mae
                                     0.003, 0.01, 0.03, 0.1,
                                     0.3, 0.5, 0.75, 1, 1.5,
                                     2, 3, 4, 6, 8),
                           target_column = 2,
                           silent=T)
  
  #get the best theta value for weigthed localised regression
  best_theta<-test_theta$theta[test_theta$mae==min(test_theta$mae)]
  
  
  jacobians <- list(matrix(0, (2), (2) ))                # to store jacobians for every timepoint
  for(tp in 2:timerange) {                             
    jacobians[[tp]] <- matrix(0, (2), (2))
  }
  
  for (targ_col in 1:(2)) {                          # targ_col =  column of dependent species, animals before plants in ts_data
    #print(targ_col)
    Embedding <- colnames                           # change 'colnames' to specific columns if necessary
    Edim <- length(Embedding)
    
    ts_data <- ts_data[, Embedding]
    
    coeff_names <- sapply(colnames(ts_data), function(x) paste("d", colnames(ts_data)[targ_col], "d", x, sep = ""))
    
    block <- cbind(ts_data[2:dim(ts_data)[1], targ_col], ts_data[1:(dim(ts_data)[1] - 1),])
    norm_consts <- apply(block, 2, function(x) sd(x,na.rm = T))
    block <- as.data.frame(apply(block, 2, function(x) (x - mean(x,na.rm = T)) / sd(x,na.rm = T)))
    
    
    lib <- 1:dim(block)[1]
    pred <- 1:dim(block)[1]
    
    theta <-best_theta
    
    coeff <- array(0, dim = c(length(pred), Edim))
    colnames(coeff) <- coeff_names
    coeff <- as.data.frame(coeff)
    
    for (ipred in 1:length(pred)) {  # partial derivatives for the jacobian matrix
      libs = lib[-pred[ipred]]
      q <- matrix(as.numeric(block[pred[ipred], 2:dim(block)[2]]),
                  ncol = Edim, nrow = length(libs), byrow = TRUE)
      
      distances <- sqrt(rowSums((block[libs, 2:dim(block)[2]] - q)^2))
      dbar <- mean(distances, na.rm=T)
      Ws <- exp(-theta*distances/dbar)
      svd_fit <- lm_svdsolve(y=block[libs, 1], x = block[libs, 2:dim(block)[2]], ws = Ws)  
      jacobians[[ipred  + 1]][targ_col,] <- coeff[ipred, ] <- svd_fit[-1]              # calculate partial derivatives and store in jacobian matrix for time point 'ipred'
    }                                                                                                            # partial derivatives of targ_col with respect to other cols in coeff[ipred, ]
  }   
  dom_eigen<- numeric() 
  
  for (j in 1:timerange) {
    dom_eigen[j] <- Re(eigen(jacobians[[j]], only.values = T)$values[1])  # abs(Re(eigen(jacobians[[j]])$values[,1]))[i]
  }
  
  return(list(dom_eigen=dom_eigen,theta=best_theta))
}

sessionInfo()


#packages needed
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(rEDM)
library(RColorBrewer)

# Read in aphid–ladybird time series dataset
t_dat1 <- read.csv(file = "Data_aphid_ladybird.csv")
str(t_dat1)

t_dat <-subset(t_dat1,t_dat1$Sampling_day>-1)
table(t_dat$Treatment,t_dat$Sampling_day)

# very important to sort by sampling day
t_dat <- t_dat %>% arrange(Cage_ID, Sampling_day)


# excluding certain cages
cage_exclude<-t_dat %>%
  group_by(Cage_ID,Treatment) %>%
  summarise(
    aphid_sd = sd(Aphids, na.rm = TRUE),
    ladybird_sd = sd(Total_ladybirds, na.rm = TRUE)
  ) %>%
  filter(aphid_sd == 0 | ladybird_sd == 0)  # These are problematic
cage_exclude

bad_cages <- cage_exclude$Cage_ID # c(4, 7, 8, 9, 12, 13, 16) #these cages are excluded

# # taking out specific cages of 4, 7, 8, 9, 12, 13, 16 which have state data that does not fit the criteria for S-map
t_dat <- t_dat %>% filter(!Cage_ID %in% bad_cages)

# no ladybug abudnance fluctuations in these cages so excluded
ggplot(t_dat1 %>% filter(Cage_ID %in% bad_cages),
       aes(x = Sampling_day, y = Total_ladybirds, 
           group = Cage_ID, color = as.factor(Cage_ID))) +
  geom_line() +
  facet_grid(Cage_ID~.)+
  labs(title = "Ladybird Time Series in Problematic Cages", color = "Cage ID")

dat_timerange<-expand.grid(timerange=c(24))


################################################## For Mix Treatment ######################################
#for time points 24 we estimate dominant eigenvalue using this for loop
dom_output<-NULL
for(i in 1:nrow(dat_timerange)){  
  
  unique(t_dat$Cage_ID)
  
  timerange<-dat_timerange$timerange[i]
  
  tdat_Mix<-t_dat %>% filter(Treatment == "Mix",Sampling_day>=0)
  
  cages<-unique(tdat_Mix$Cage_ID)
  OUT<-NULL
  for(r in cages){
    tryCatch({ # by passing small errors if arises
      
      temp<-tdat_Mix %>% filter(Cage_ID == r)
      
      d_eigenval<-dom_eigenvalue_EDM(timerange = timerange,tdat = temp)$dom_eigen
      
    }, error=function(e){
      cat("error", r, ":", e$message, "\n")
    })
    
    #storing in dataframe
    temp_dat  <-as.data.frame(cbind(d_eigenval,
                                    (1:timerange),
                                    rep(temp$Position_of_cage[1],each=length(d_eigenval)),
                                    rep(temp$No_of_leaves[1],each=length(d_eigenval)),
                                    rep(temp$Plant_Age[1],each=length(d_eigenval)),
                                    rep(as.character(temp$Treatment[1]),each=length(d_eigenval)),
                                    rep(r,each=length(d_eigenval))))
    
    OUT <- rbind(OUT,temp_dat)  
  }
  
  #storing all that in a dataframe and renaming the columns
  colnames(OUT) <- c("d_eigenvalue", "time","Position_of_cage","No_of_leaves","Plant_Age", "Treatment", "Cage_ID")
  OUT$d_eigenvalue<-as.numeric(OUT$d_eigenvalue)
  OUT$time<-as.numeric(OUT$time)
  OUT$Treatment<-as.factor(OUT$Treatment)
  OUT$Cage_ID<-as.factor(OUT$Cage_ID)
  OUT$Position_of_cage<-as.numeric(OUT$Position_of_cage)
  OUT$No_of_leaves<-as.numeric(OUT$No_of_leaves)
  OUT$Plant_Age<-as.numeric(OUT$Plant_Age)
  
  mean_eigenval_MIX<-OUT %>%select(Cage_ID,d_eigenvalue,Treatment,time,Position_of_cage,No_of_leaves,Plant_Age) %>% 
    group_by(Treatment,Cage_ID,No_of_leaves,Plant_Age,Position_of_cage) %>% 
    summarise(mean_eigen= mean(d_eigenvalue,na.rm = T))
  mean_eigenval_MIX
  OUT %>% ggplot(aes(x=time, d_eigenvalue,color=Cage_ID))+
    geom_point(size=2)+
    theme_classic()+
    ylim(c(-100,100))
  
  ####################################### For Dropper Treatment ################################
  tdat_drp<-t_dat %>% filter(Treatment == "Dropper",Sampling_day>=0)
  cages<-unique(tdat_drp$Cage_ID)
  OUT_drp<-NULL
  for(r in cages){
    tryCatch({
      
      temp<-tdat_drp %>% filter(Cage_ID == r)
      
      d_eigenval<-dom_eigenvalue_EDM(timerange = timerange,tdat = temp)$dom_eigen
      
    }, error=function(e){
      cat("error", r, ":", e$message, "\n")
    })
    
    temp_dat  <-as.data.frame(cbind(d_eigenval,
                                    (1:timerange),
                                    rep(temp$Position_of_cage[1],each=length(d_eigenval)),
                                    rep(temp$No_of_leaves[1],each=length(d_eigenval)),
                                    rep(temp$Plant_Age[1],each=length(d_eigenval)),
                                    rep(as.character(temp$Treatment[1]),each=length(d_eigenval)),
                                    rep(r,each=length(d_eigenval))))
    
    OUT_drp <- rbind(OUT_drp,temp_dat)  
  }
  
  colnames(OUT_drp) <-  c("d_eigenvalue", "time","Position_of_cage","No_of_leaves","Plant_Age", "Treatment", "Cage_ID")
  OUT_drp$d_eigenvalue<-as.numeric(OUT_drp$d_eigenvalue)
  OUT_drp$time<-as.numeric(OUT_drp$time)
  OUT_drp$Treatment<-as.factor(OUT_drp$Treatment)
  OUT_drp$Cage_ID<-as.factor(OUT_drp$Cage_ID)
  OUT_drp$Position_of_cage<-as.numeric(OUT_drp$Position_of_cage)
  OUT_drp$No_of_leaves<-as.numeric(OUT_drp$No_of_leaves)
  OUT_drp$Plant_Age<-as.numeric(OUT_drp$Plant_Age)
  mean_eigenval_drp<-OUT_drp %>% select(Cage_ID,d_eigenvalue,Treatment,time,Position_of_cage,No_of_leaves,Plant_Age) %>% 
    group_by(Treatment,Cage_ID,No_of_leaves,Plant_Age,Position_of_cage) %>% 
    summarise(mean_eigen= mean(d_eigenvalue,na.rm = T))
  mean_eigenval_drp
  OUT_drp %>% ggplot(aes(x=time, d_eigenvalue,color=Cage_ID))+
    geom_point(size=2)+
    theme_classic()+
    ylim(c(-100,100))
  
  
  ################################################ NON_droppers ##################################################
  
  tdat_Ndrp<-t_dat %>% filter(Treatment == "Non-dropper",Sampling_day>=0)
  cages<-unique(tdat_Ndrp$Cage_ID)
  OUT_Ndrp<-NULL
  for(r in cages){
    tryCatch({
      
      temp<-tdat_Ndrp %>% filter(Cage_ID == r)
      
      d_eigenval<-dom_eigenvalue_EDM(timerange = timerange,tdat = temp)$dom_eigen
      
    }, error=function(e){
      cat("error", r, ":", e$message, "\n")
    })
    
    temp_dat  <-as.data.frame(cbind(d_eigenval,
                                    (1:timerange),
                                    rep(temp$Position_of_cage[1],each=length(d_eigenval)),
                                    rep(temp$No_of_leaves[1],each=length(d_eigenval)),
                                    rep(temp$Plant_Age[1],each=length(d_eigenval)),
                                    rep(as.character(temp$Treatment[1]),each=length(d_eigenval)),
                                    rep(r,each=length(d_eigenval))))
    
    OUT_Ndrp <- rbind(OUT_Ndrp,temp_dat)  
  }
  
  colnames(OUT_Ndrp) <-  c("d_eigenvalue", "time","Position_of_cage","No_of_leaves","Plant_Age", "Treatment", "Cage_ID")
  OUT_Ndrp$d_eigenvalue<-as.numeric(OUT_Ndrp$d_eigenvalue)
  OUT_Ndrp$time<-as.numeric(OUT_Ndrp$time)
  OUT_Ndrp$Treatment<-as.factor(OUT_Ndrp$Treatment)
  OUT_Ndrp$Cage_ID<-as.factor(OUT_Ndrp$Cage_ID)
  OUT_Ndrp$Position_of_cage<-as.numeric(OUT_Ndrp$Position_of_cage)
  OUT_Ndrp$No_of_leaves<-as.numeric(OUT_Ndrp$Cage_ID)
  OUT_Ndrp$Plant_Age<-as.numeric(OUT_Ndrp$Plant_Age)
  
  mean_eigenval_Ndrp<-OUT_Ndrp %>% select(Cage_ID,d_eigenvalue,Treatment,time,Position_of_cage,No_of_leaves,Plant_Age) %>% 
    group_by(Treatment,Cage_ID,No_of_leaves,Plant_Age,Position_of_cage) %>% 
    summarise(mean_eigen= mean(d_eigenvalue,na.rm = T))
  mean_eigenval_Ndrp
  OUT_Ndrp %>% ggplot(aes(x=time, d_eigenvalue,color=Cage_ID))+
    geom_point(size=2)+
    theme_classic()+
    ylim(c(-100,100))
  
  mean_stability_local<-rbind(mean_eigenval_drp,mean_eigenval_MIX,mean_eigenval_Ndrp)
  
  mean_stability_local$timerange<-dat_timerange$timerange[i]
  
  dom_output<-rbind(dom_output,mean_stability_local)
}

# Do S-map analysis with the best theta  and create this funciton
#timerange: total time points of timeseries which is 24
# tdat : temporary dataframe
# cages : cage id
# Treatment : Either Mix, non-dropper, or Dropper
#targ_col : target column, either aphid or ladybird 
interaction_str<-function(timerange, tdat, cages,Treatment,targ_col){
  OUT<-NULL
  
  for(r in cages){
    
    
    tdatx<-tdat %>% filter(Cage_ID == r)
    tdatx$ladybugs_all<-(tdatx$Ladybird_larvae+tdatx$Ladybird_adults)
    
    ts_data<-cbind(tdatx$Aphids[1:timerange],tdatx$ladybugs_all[1:timerange])
    colnames <- c()
    
    for(i in 1:1) {
      colnames <- c(colnames, paste("A", i, sep = ""))
    }
    
    for(i in 1:1) {
      colnames <- c(colnames, paste("P", i, sep = ""))
    }
    
    colnames(ts_data) <- colnames
    data_used<-1:timerange
    targ_col <- targ_col
    Embedding <- colnames   
    block <- ts_data[,Embedding]
    block <- as.data.frame(apply(block, 2, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
    
    test_theta <- block_lnlp(block[data_used,],
                             method = "s-map",
                             num_neighbors = 0, # We have to use any value <
                             theta = c(0,  0.001,
                                       0.003, 0.01, 0.03, 0.1,
                                       0.3, 0.5, 0.75, 1, 1.5,
                                       2, 3, 4, 6, 8),
                             target_column = targ_col,
                             silent=T)
    best_theta<- test_theta[which.min(test_theta$mae),"theta"]
    # Define the target column (C1 = column 2)
    best_theta
    
    
    smap_res <-block_lnlp(block[data_used,],
                          method = "s-map",
                          num_neighbors = 0, # we have to use any value < 1
                          theta = best_theta,
                          target_column = targ_col,
                          silent = T,
                          save_smap_coefficients = T)
    
    smap_out <- as.data.frame(smap_res$model_output[[1]])
    
    ## Time series of fluctuating interaction strength
    smap_coef <- as.data.frame(smap_res$smap_coefficients[[1]])
    colnames(smap_coef) <- c("dN_dN","dN_dP","Intercept")
    smap_coef$Treatment <- Treatment
    
    temp_dat  <-as.data.frame(cbind(smap_coef$dN_dN,
                                    smap_coef$dN_dP,
                                    rep(best_theta,each=length(smap_coef$dN_dN)),
                                    rep(r,each=length(smap_coef$dN_dN)),
                                    smap_coef$Intercept,
                                    (1:length(smap_coef$dN_dN)),
                                    rep(as.character(smap_coef$Treatment[1]),each=length(smap_coef$dN_dN)),
                                    rep(r,each=length(smap_coef$dN_dN))))
    
    OUT <- rbind(OUT,temp_dat)  
    
  }
  return(OUT)
}


timerange<-24 #total time points

#mix treatment 
tdat_mix<-t_dat %>% filter(Treatment == "Mix",Sampling_day>=0)
unique(tdat_mix$Cage_ID) # cages 9,12,15 needs to be removed
out_mix<-interaction_str(timerange =timerange,tdat = tdat_mix,cages = unique(tdat_mix$Cage_ID),
                         Treatment = "Mix",targ_col = 1)


# dropper treatment
tdat_drop<-t_dat %>% filter(Treatment == "Dropper",Sampling_day>=0)
unique(tdat_drop$Cage_ID)
out_drop<- interaction_str(timerange =timerange,tdat = tdat_drop,
                           cages = unique(tdat_drop$Cage_ID),Treatment = "Dropper",targ_col = 1)

#non-dropper treatment
tdat_nondrop<-t_dat %>% filter(Treatment == "Non-dropper",Sampling_day>=0)
unique(tdat_nondrop$Cage_ID)
out_nondrop<- interaction_str(timerange =timerange,tdat = tdat_nondrop,cages = unique(tdat_nondrop$Cage_ID),
                              Treatment = "Non-dropper",targ_col=1)



Output<-as.data.frame(rbind(out_nondrop,out_drop,out_mix))
colnames(Output)<- c("dn_dn","dn_dP","theta", "Cage_ID1", "intercept","time","Treatment","Cage_ID")
Output$dn_dn<-as.numeric(Output$dn_dn)
Output$dn_dP<-as.numeric(Output$dn_dP)
Output$intercept<-as.numeric(Output$intercept)
Output$time<-as.numeric(Output$time)
Output$Treatment<-as.factor(Output$Treatment)

Output <- Output %>%
  filter_all(all_vars(!is.nan(.)))


#example dn/dn i.e., change in state densities divided change in state densities, jacobian elements over time
Output %>% ggplot(aes(x=time, y=dn_dP, col=factor(Cage_ID),group=factor(Cage_ID)))+
  geom_line()+
  theme_classic()+
  geom_hline(yintercept = 0,linetype="dashed")+
  facet_wrap(.~Treatment)

#example dn/dn i.e., change in state densities divided change in state densities, jacobian elements over time
Output %>% ggplot(aes(x=time, y=dn_dn, col=factor(Cage_ID),group=factor(Cage_ID)))+
  geom_line()+
  theme_bw()+
  geom_hline(yintercept = 0,linetype="dashed")+
  facet_wrap(.~Treatment)



#mix
tdat_mix<-t_dat %>% filter(Treatment == "Mix",Sampling_day>0)
unique(tdat_mix$Cage_ID) # cages 9,12,15 needs to be removed
out_mix_p<-interaction_str(timerange =timerange,tdat = tdat_mix,cages = unique(tdat_mix$Cage_ID),
                           Treatment = "Mix",targ_col = 2)
#out_mix<-OUT

# dropper
tdat_drop<-t_dat %>% filter(Treatment == "Dropper",Sampling_day>0)
unique(tdat_drop$Cage_ID)
out_drop_p<- interaction_str(timerange =timerange,tdat = tdat_drop,
                             cages = unique(tdat_drop$Cage_ID),Treatment = "Dropper",targ_col = 2)

#non-dropper
tdat_nondrop<-t_dat %>% filter(Treatment == "Non-dropper",Sampling_day>0)
unique(tdat_nondrop$Cage_ID)
out_nondrop_p<- interaction_str(timerange =timerange,tdat = tdat_nondrop,cages = unique(tdat_nondrop$Cage_ID),
                                Treatment = "Non-dropper",targ_col=2)

Output2<-as.data.frame(rbind(out_nondrop_p,out_drop_p,out_mix_p))
colnames(Output2)<- c("dPdn","dP_dP","theta", "Cage_ID1", "intercept","time","Treatment","Cage_ID")
Output2$dPdn<-as.numeric(Output2$dPdn)
Output2$dP_dP<-as.numeric(Output2$dP_dP)
Output2$intercept<-as.numeric(Output2$intercept)
Output2$time<-as.numeric(Output2$time)
Output2$Treatment<-as.factor(Output2$Treatment)

Output2 <- Output2 %>%
  filter_all(all_vars(!is.nan(.)))
Output2 %>% ggplot(aes(x=time, y=dPdn, col=factor(Cage_ID),group=factor(Cage_ID)))+
  geom_line()+
  theme_bw()+
  labs(y = expression("Jacobian element," *partialdiff*P / partialdiff*N))+
  geom_hline(yintercept = 0,linetype="dashed")+
  facet_wrap(.~Treatment)

Output2 %>% ggplot(aes(x=time, y=dP_dP, col=factor(Cage_ID),group=factor(Cage_ID)))+
  geom_line()+
  theme_bw()+
  labs(y = expression("Jacobian element," *partialdiff*P / partialdiff*P))+
  geom_hline(yintercept = 0,linetype="dashed")+
  facet_wrap(.~Treatment)


dndn_mean<-Output %>% select(Treatment, Cage_ID, dn_dn,time) #%>% 

dndp_mean<-Output %>% select(Treatment, Cage_ID, dn_dP,time)# %>% 


dpdn_mean<-Output2 %>% select(Treatment, Cage_ID, dPdn,time) #%>% 

dpdp_mean<-Output2 %>% select(Treatment, Cage_ID, dP_dP,time) #%>% 



combined_data <- dndn_mean %>% 
  full_join(dndp_mean,by=c("Cage_ID","Treatment","time")) %>% 
  full_join(dpdn_mean,by=c("Cage_ID","Treatment","time")) %>% 
  full_join(dpdp_mean,by=c("Cage_ID","Treatment","time")) 

# Create a list-column with matrices for each Treatment and Cage_ID
matrix_data <- combined_data %>%
  group_by(Treatment, Cage_ID,time) %>%
  dplyr::summarise(
    result_matrix = list(matrix(c(dn_dn, dn_dP, dPdn, dP_dP), nrow = 2, byrow = TRUE)),
    .groups = "drop"
  )

# Calculate dominant_eigenvalue and avg_robustness for each matrix
# Load necessary libraries
library(dplyr)
library(purrr)
library(dplyr)
matrix_data <- matrix_data %>%
  group_by(Cage_ID, Treatment, time) %>%  # Group by Time if it exists in the data
  mutate(
    dominant_eigenvalue = map_dbl(result_matrix, function(mat) {
      # Check if matrix is square, finite, and has no missing values
      if (is.matrix(mat) && nrow(mat) == ncol(mat) && all(is.finite(mat))) {
        Re(eigen(mat, only.values = TRUE)$values[1])
      } else {
        NA  # Return NA if matrix is not suitable for eigen calculation
      }
    }),
    avg_robustness = map_dbl(result_matrix, function(mat) {
      # Check if matrix is square, finite, and has no missing values
      if (is.matrix(mat) && nrow(mat) == ncol(mat) && all(is.finite(mat))) {
        e_values <- eigen(mat, only.values = TRUE)$values #eigen analysis
        exp(mean(log(abs(e_values))))
      } else {
        NA  # Return NA if matrix is not suitable for eigen calculation
      }
    })
  ) %>%
  ungroup()

matrix_data

# Plot avg_robustness over time, colored by Treatment
ggplot(matrix_data, aes(x = time, y = avg_robustness, color = Treatment, group = Cage_ID)) +
  geom_line(aes(linetype = Treatment)) +  # Lines for each Cage_ID, differentiated by Treatment
  geom_point(size = 1) +                  # Facet by Treatment for easier comparison
  labs(
    title = "",
    x = "Sampling day",
    y = "Average Robustness"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Plot dominant eigen over time, colored by Treatment
ggplot(matrix_data, aes(x = time, y = dominant_eigenvalue, color = Treatment, group = Cage_ID)) +
  geom_line(aes(linetype = Treatment)) +  # Lines for each Cage_ID, differentiated by Treatment
  geom_point(size = 1) +   
  geom_hline(yintercept = 0,linetype="dashed",color="black")+# Points for each Cage_ID over time
  facet_wrap(~ Treatment) +     
  # Facet by Treatment for easier comparison
  labs(
    title = "Dominant eigenvalue ",
    x = "Sampling day",
    y = expression(Re(lambda))
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


se_robustness<-matrix_data %>%
  group_by(time, Treatment) %>% dplyr::summarise(se_robustness=sd(avg_robustness,na.rm=T)/sqrt(n()))

average_data <- matrix_data %>%
  group_by(time, Treatment) %>% 
  dplyr::summarise(
    avg_robustness = mean(avg_robustness, na.rm = TRUE),
    avg_eigenvalue=mean(dominant_eigenvalue,na.rm=TRUE),
    se_eigen = sd(dominant_eigenvalue, na.rm = TRUE) / sqrt(n())) %>%
  ungroup()

average_data$se_robustness<-se_robustness$se_robustness




# Plot the average avg_robustness over time with SE bands for each Treatment
(p1<-ggplot(average_data, aes(x = time, y = avg_robustness, color = Treatment, group = Treatment)) +
    geom_line(size = 1.2) +                                 # Line for mean avg_robustness
    geom_point(size = 2) +                                  # Points for mean values
    geom_ribbon(aes(ymin = avg_robustness - (1.96 *se_robustness), 
                    ymax = avg_robustness + (1.96 *se_robustness), fill = Treatment), 
                alpha = 0.2, linetype = 0) +                # SE bands
    labs(
      x = "Sampling day",
      y = "Average Robustness"
    ) +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2", name = "Treatment") +
    scale_color_brewer(palette = "Dark2", name = "Treatment")+
    scale_color_manual(values = c("#B2182B","#FDDBC7","#2166AC"))+
    scale_fill_manual(values = c("#B2182B","#FDDBC7","#2166AC"))  +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ))


(p2<-ggplot(average_data, aes(x = time, y = avg_eigenvalue, color = Treatment, group = Treatment)) +
    geom_line(size = 1.2) +                                 # Line for mean avg_robustness
    geom_point(size = 2) +                                  # Points for mean values
    geom_ribbon(aes(ymin = avg_eigenvalue - (1.96 *se_eigen),
                    ymax = avg_eigenvalue + (1.96 *se_eigen), fill = Treatment), 
                alpha = 0.2, linetype = 0) +                # SE bands
    labs(
      x = "Sampling day",
      y = expression(Re(lambda))
    ) +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2", name = "Treatment") +
    scale_color_brewer(palette = "Dark2", name = "Treatment") +
    geom_hline(yintercept = 0,linetype="dashed")+
    scale_color_manual(values = c("#B2182B","#FDDBC7","#2166AC"))+
    scale_fill_manual(values = c("#B2182B","#FDDBC7","#2166AC")) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ))


ggpubr::ggarrange(p2,p1,labels = c("A","B"))

# Arrange the plots with a shared legend at the bottom
# png(file = "Figure 5.png",width = 8 ,height = 4,units = 'in', res = 600)
ggpubr::ggarrange(
  p2, p1,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "bottom"
)
# dev.off()

# Dominant eigenvalue plot (individual lines + mean)

# png(file = "Figure Supple S2.png",width = 8 ,height = 4,units = 'in', res = 600)
ggplot(matrix_data, aes(x = time, y = dominant_eigenvalue, group = Cage_ID, color = Treatment)) +
  geom_line(alpha = 0.3, size = 0.8) +  # Individual cage trajectories
  geom_smooth(aes(group = Treatment, fill = Treatment),
              method = "loess", se = TRUE, size = 1.2, alpha = 0.3) +  # LOESS with 95% CI
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Treatment) +
  scale_color_manual(values = c("#B2182B","#FDDBC7","#2166AC"))+
  scale_fill_manual(values = c("#B2182B","#FDDBC7","#2166AC")) +
  theme_classic() +
  labs(
    title = "Individual Cage Dynamics: Dominant Eigenvalue (LOESS ± 95% CI)",
    y = expression(Re(lambda)),
    x = "Sampling Day"
  )
# dev.off()

library(dplyr)
library(ggplot2)

# Average robustness plot (individual lines + mean)


# png(file = "Figure Supple S3.png",width = 8 ,height = 4,units = 'in', res = 600)
ggplot(matrix_data, aes(x = time, y = avg_robustness, group = Cage_ID, color = Treatment)) +
  geom_line(alpha = 0.3, size = 0.8) +  # Individual cage trajectories
  geom_smooth(aes(group = Treatment, fill = Treatment),
              method = "loess", se = TRUE, size = 1.2, alpha = 0.3) +  # LOESS with 95% CI
  facet_wrap(~ Treatment) +
  scale_color_manual(values = c("#B2182B","#FDDBC7","#2166AC"))+
  scale_fill_manual(values = c("#B2182B","#FDDBC7","#2166AC")) +
  theme_classic() +
  labs(
    title = "Individual Cage Dynamics: Robustness (LOESS ± 95% CI)",
    y = "Average Robustness",
    x = "Sampling Day"
  )
# dev.off()


# > sessionInfo()
# R version 4.3.3 (2024-02-29 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 26100)
# 