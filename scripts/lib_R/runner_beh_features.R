# ~~~~~~~~~~~~ library
#install.packages("hBayesDM", dependencies=TRUE)
#Sys.setenv(BUILD_ALL='true')  # Build all the models on installation
#Sys.setenv(MAKEFLAGS='-j 4')  # Use 4 cores for compilation (or the number you want)
#install.packages("hBayesDM")
#install.packages("gridExtra")  
#install.packages("psycho")
library(readr)
library(ggplot2)
library(tidyverse)
library(tidyquant)
library(ggdist)
library(ggthemes)
library(dplyr)



# ~~~~~~~~~~~ download the data
sessions = c('ses-1', 'ses-2', 'ses-3', 'ses-4')
fdir_get_across <-
  "/data/p_02825/cocoa/data/derivatives/stats/group/"
for (session in sessions) {
  print(session)
  
  fdir_orig <- gsub(" ", "", paste(fdir_get_across, session))
  fdir = paste(fdir_orig, '/beh/')
  fdir <- gsub(" ", "", fdir)

  if (dir.exists(fdir)) {
    list_of_files <- list.files(
      path = fdir,
      recursive = FALSE,
      pattern = "\\.csv$",
      full.names = TRUE
    )
    
    
    concatenated_df <- c()
    for (i in 1:length(list_of_files)) {
      variable_names <- c("tmp_var")
      df <- readr::read_csv(list_of_files[i], id = "file_name")
      assign(variable_names[1], df)
      print(df)
      
      # ~~~~~~~~~~~ pick up the correct, reaction time, participant name, letter, stimulus stype
      tmp_var_relevant = tmp_var[c("key_resp.corr",
                                   "key_resp.rt",
                                   "participant",
                                   "stim_images_order",
                                   "type")]
      tmp_var_relevant = tmp_var_relevant[complete.cases(tmp_var_relevant),] # remove NaN
      
      
      # combine all the subjects
      concatenated_df <- rbind(concatenated_df, tmp_var_relevant)
    }
  
  
  concatenated_df$stim_images_order <-
    gsub("images/alphabets/B.png",
         "B",
         concatenated_df$stim_images_order)
  concatenated_df$stim_images_order <-
    gsub("images/alphabets/F.png",
         "F",
         concatenated_df$stim_images_order)
  concatenated_df$stim_images_order <-
    gsub("images/alphabets/K.png",
         "K",
         concatenated_df$stim_images_order)
  concatenated_df$stim_images_order <-
    gsub("images/alphabets/M.png",
         "M",
         concatenated_df$stim_images_order)
  concatenated_df$stim_images_order <-
    gsub("images/alphabets/P.png",
         "P",
         concatenated_df$stim_images_order)
  concatenated_df$stim_images_order <-
    gsub("images/alphabets/T.png",
         "T",
         concatenated_df$stim_images_order)
  
  
  fdir_proc <- gsub(" ", "", paste(fdir, "proc_data/"))
  dir.create(fdir_proc)
  
  fname <- gsub(" ", "", paste(fdir_proc, "Beh_data_concatenated.csv"))
  write.csv(concatenated_df, fname)
  
  dr_across_session <-  data_frame(type = c('Total'),
                                   RT = c(concatenated_df$key_resp.rt))
  dr_across_session1 <-  data_frame(type = c('Congruent'),
                                    RT = c(concatenated_df$key_resp.rt[concatenated_df$type ==
                                                                         'congruent']))
  dr_across_session2 <-  data_frame(type = c('Incongruent'),
                                    RT = c(concatenated_df$key_resp.rt[concatenated_df$type ==
                                                                         'incongruent']))
  dr_across_session <-
    rbind(dr_across_session, dr_across_session1, dr_across_session2)
  
  fname <- gsub(" ", "", paste(fdir_proc, "Beh_RT_concatenated.csv"))
  write.csv(dr_across_session, fname)
  # transform stim_images to character
  colnames(concatenated_df)[colnames(concatenated_df) == "stim_images_order"] <-
    "character"
  
  
  
  
  
  # ~~~~~~~~~~~ plot an average accuracy for each participant or each letter ~~~~~~~~~~~~~~~~~~
  # summary
  # Define the data_summary function
  data_summary <- function(data, varname, groupnames) {
    require(plyr)
    summary_func <- function(x, col) {
      c(mean = mean(x[[col]], na.rm = TRUE),
        sd = sd(x[[col]], na.rm = TRUE))
    }
    data_sum <- ddply(data, groupnames, .fun = summary_func, varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  
  #summary result
  df1 <-
    data_summary(concatenated_df,
                 varname = "key_resp.corr",
                 groupnames = c("participant"))
  
  # standard error
  df1$SE <- df1$sd / sqrt(48)
  
  type_row <- data.frame(Type = c('Total'))
  df1 <- cbind(df1, type_row)
  #summary result
  tmp_df <-
    data.frame(key_resp.corr = concatenated_df$key_resp.corr[concatenated_df$type ==
                                                               'congruent'],
               participant = concatenated_df$participant[concatenated_df$type ==
                                                           'congruent'])
  df1_cong <-
    data_summary(tmp_df,
                 varname = "key_resp.corr",
                 groupnames = c("participant"))
  
  # standard error
  df1_cong$SE <- df1_cong$sd / sqrt(48)
  type_row <- data.frame(Type = c('Congruent'))
  df1_cong <- cbind(df1_cong, type_row)
  
  #summary result
  tmp_df <-
    data.frame(key_resp.corr = concatenated_df$key_resp.corr[concatenated_df$type ==
                                                               'incongruent'],
               participant = concatenated_df$participant[concatenated_df$type ==
                                                           'incongruent'])
  df1_incong <-
    data_summary(tmp_df,
                 varname = "key_resp.corr",
                 groupnames = c("participant"))
  
  # standard error
  df1_incong$SE <- df1_incong$sd / sqrt(48)
  
  type_row <- data.frame(Type = c('Incongruent'))
  df1_incong <- cbind(df1_incong, type_row)
  
  # Combine all types
  df1 <- rbind(df1, df1_cong, df1_incong)
  
  fname <-
    gsub(" ",
         "",
         paste(fdir_proc, "Beh_data_accuracy_concatenated.csv"))
  write.csv(df1, fname)
  
  
  # ~~~~~~~~~~~ (2.1) Group-average reaction time
  df_averageRT <-
    data.frame(
      x_axis = c('1'),
      Mean_RT = mean(concatenated_df$key_resp.rt),
      SD = sd(concatenated_df$key_resp.rt)
    )
  df_averageRT$SE <- df_averageRT$SD / sqrt(nrow(concatenated_df))
  
  Mean_RT <- df_averageRT$Mean_RT
  SE <- df_averageRT$SE
  SD <- df_averageRT$SD
  
  fname <-
    gsub(" ",
         "",
         paste(fdir_proc, "/Beh_data_avg_RT_concatenated.csv"))
  write.csv(df_averageRT, fname)
  
  
  #~~~~~~~~~~~ signal detection theory ~~~~~~~~~~~~~~~~~~
  library(psycho)
  
  SubName <- concatenated_df$participant
  SubLists <- unique(SubName)
  
  # signal detection theory B
  n_hit_B <- c()
  n_fa_B <- c()
  n_miss_B <- c()
  n_cr_B <- c()
  
  for (i in 1:length(list_of_files)) {
    N_hit <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'B', ])
    n_hit_B <- c(n_hit_B, N_hit)
    N_fa <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'B', ])
    n_fa_B <- c(n_fa_B, N_fa)
    N_miss <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'B', ])
    n_miss_B <- c(n_miss_B, N_miss)
    N_cr <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'B', ])
    n_cr_B <- c(n_cr_B, N_cr)
  }
  indices_B <- psycho::dprime(n_hit_B, n_fa_B, n_miss_B, n_cr_B)
  
  # signal detection theory F
  n_hit_F <- c()
  n_fa_F <- c()
  n_miss_F <- c()
  n_cr_F <- c()
  
  for (i in 1:length(list_of_files)) {
    N_hit <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'F', ])
    n_hit_F <- c(n_hit_F, N_hit)
    N_fa <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'F', ])
    n_fa_F <- c(n_fa_F, N_fa)
    N_miss <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'F', ])
    n_miss_F <- c(n_miss_F, N_miss)
    N_cr <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'F', ])
    n_cr_F <- c(n_cr_F, N_cr)
  }
  indices_F <- psycho::dprime(n_hit_F, n_fa_F, n_miss_F, n_cr_F)
  
  # signal detection theory K
  n_hit_K <- c()
  n_fa_K <- c()
  n_miss_K <- c()
  n_cr_K <- c()
  
  for (i in 1:length(list_of_files)) {
    N_hit <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'K', ])
    n_hit_K <- c(n_hit_K, N_hit)
    N_fa <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'K', ])
    n_fa_K <- c(n_fa_K, N_fa)
    N_miss <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'K', ])
    n_miss_K <- c(n_miss_K, N_miss)
    N_cr <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'K', ])
    n_cr_K <- c(n_cr_K, N_cr)
  }
  indices_K <- psycho::dprime(n_hit_K, n_fa_K, n_miss_K, n_cr_K)
  
  # signal detection theory M
  n_hit_M <- c()
  n_fa_M <- c()
  n_miss_M <- c()
  n_cr_M <- c()
  
  for (i in 1:length(list_of_files)) {
    N_hit <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'M', ])
    n_hit_M <- c(n_hit_M, N_hit)
    N_fa <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'M', ])
    n_fa_M <- c(n_fa_M, N_fa)
    N_miss <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'M', ])
    n_miss_M <- c(n_miss_M, N_miss)
    N_cr <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'M', ])
    n_cr_M <- c(n_cr_M, N_cr)
  }
  indices_M <- psycho::dprime(n_hit_M, n_fa_M, n_miss_M, n_cr_M)
  
  # signal detection theory P
  n_hit_P <- c()
  n_fa_P <- c()
  n_miss_P <- c()
  n_cr_P <- c()
  
  for (i in 1:length(list_of_files)) {
    N_hit <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'P', ])
    n_hit_P <- c(n_hit_P, N_hit)
    N_fa <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'P', ])
    n_fa_P <- c(n_fa_P, N_fa)
    N_miss <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'P', ])
    n_miss_P <- c(n_miss_P, N_miss)
    N_cr <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'P', ])
    n_cr_P <- c(n_cr_P, N_cr)
  }
  
  indices_P <- psycho::dprime(n_hit_P, n_fa_P, n_miss_P, n_cr_P)
  
  # signal detection theory T
  n_hit_T <- c()
  n_fa_T <- c()
  n_miss_T <- c()
  n_cr_T <- c()
  
  for (i in 1:length(list_of_files)) {
    N_hit <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'T', ])
    n_hit_T <- c(n_hit_T, N_hit)
    N_fa <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'T', ])
    n_fa_T <- c(n_fa_T, N_fa)
    N_miss <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent' &
                             concatenated_df$character == 'T', ])
    n_miss_T <- c(n_miss_T, N_miss)
    N_cr <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent' &
                             concatenated_df$character == 'T', ])
    n_cr_T <- c(n_cr_T, N_cr)
  }
  
  indices_T <- psycho::dprime(n_hit_T, n_fa_T, n_miss_T, n_cr_T)
  
  # signal detection theory TOTAL
  n_hit1 <- c()
  n_fa1 <- c()
  n_miss1 <- c()
  n_cr1 <- c()
  for (i in 1:length(list_of_files)) {
    N_hit1 <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent', ])
    n_hit1 <- c(n_hit1, N_hit1)
    N_fa1 <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent', ])
    n_fa1 <- c(n_fa1, N_fa1)
    N_miss1 <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 0 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'congruent', ])
    n_miss1 <- c(n_miss1, N_miss1)
    N_cr1 <-
      nrow(concatenated_df[concatenated_df$key_resp.corr == 1 &
                             concatenated_df$participant == SubLists[i] &
                             concatenated_df$type == 'incongruent', ])
    n_cr1 <- c(n_cr1, N_cr1)
  }
  
  indices_tot <- psycho::dprime(n_hit1, n_fa1, n_miss1, n_cr1)
  
  # visualization of d and C
  df_signaldetection <-
    data.frame(
      participant = rep(SubLists, times = 6),
      character = c(
        rep("B", times = length(list_of_files)),
        rep("F", times = length(list_of_files)),
        rep("K", times = length(list_of_files)),
        rep("M", times = length(list_of_files)),
        rep("P", times = length(list_of_files)),
        rep("T", times = length(list_of_files))
      )
    )
  df_signaldetection$d <-
    c(indices_B[["dprime"]],
      indices_F[["dprime"]],
      indices_K[["dprime"]],
      indices_M[["dprime"]],
      indices_P[["dprime"]],
      indices_T[["dprime"]])
  df_signaldetection$c <-
    c(indices_B[["c"]],
      indices_F[["c"]],
      indices_K[["c"]],
      indices_M[["c"]],
      indices_P[["c"]],
      indices_T[["c"]])
  
  #
  fname <-
    gsub(" ",
         "",
         paste(fdir_proc, "/Beh_data_sensitivity_cha_concatenated.csv"))
  write.csv(df_signaldetection, fname)
  
  df_signaldetection_tot <- data.frame(
    participant = c(SubLists),
    d = c(indices_tot[["dprime"]]),
    c = c(indices_tot[["c"]])
  )
  #
  fname <-
    gsub(" ",
         "",
         paste(fdir_proc, "/Beh_data_sensitivity_tot_concatenated.csv"))
  write.csv(df_signaldetection_tot, fname)
  
  }
}







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~ COLLECT DATA FROM SEESSIONS ~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sessions = c('ses-1', 'ses-2', 'ses-3', 'ses-4')
fdir_get_across <- "/data/p_02825/cocoa/data/derivatives/stats/group/"

RT_across_session  <- data.frame(RT= c(),
                                 Session=c(),
                                 Type=c())
Acc_across_session  <- data.frame(Acc= c(),
                                  SD= c(),
                                  SE= c(),
                                  Session=c(),
                                  Type=c())
Sen_across_session  <- data.frame(Sen= c(),
                                  Sen_bias= c(),
                                  Session=c(),
                                  Type=c())
RT_avg_across_session <- data.frame(RT_avg= c(),
                                    SE=c(),
                                    SD=c(),
                                    Session=c(),
                                    Type=c())
Acc_avg_across_session <- data.frame(Acc_avg= c(),
                                     SE=c(),
                                     SD=c(),
                                     Session=c(),
                                     Type=c())
Sen_avg_across_session <- data.frame(Sen_avg= c(),
                                     SD= c(),
                                     Sen_bias_avg=c(),
                                     Session=c())

for (session in sessions){
  print(session)
  tmp_dir = paste(fdir_get_across,session,"/beh/proc_data")
  tmp_dir <- gsub(" ", "", tmp_dir)
  if (dir.exists(tmp_dir)){
    
    
    # Load the  RT
    fname1 <- paste(tmp_dir, "/Beh_RT_concatenated.csv")
    fname1 <- gsub(" ", "", fname1)
    
    ses_row <- data.frame(Session = session)
    df_rt <- readr::read_csv(fname1,id = "file_name")
    df_rt <- cbind(df_rt, session)
    df_rt <- subset(df_rt, select = -file_name)
    names(df_rt)[names(df_rt) == "key_resp.rt"] <- "RT"
    
    RT_across_session  <- rbind(RT_across_session , df_rt)
    
    # Calculate the average RT
    df_tmp_RT_avg <- data.frame(Type = c('Total', 'Congruent', 'Incongruent'),
                                RT_avg = c(mean(RT_across_session$RT),
                                           mean(RT_across_session$RT[RT_across_session$type=='Congruent']),
                                           mean(RT_across_session$RT[RT_across_session$type=='Incongruent'])),
                                SD = c(sd(RT_across_session$RT),
                                       sd(RT_across_session$RT[RT_across_session$type=='Congruent']),
                                       sd(RT_across_session$RT[RT_across_session$type=='Incongruent'])),
                                SE = c(sd(RT_across_session$RT) / sqrt(nrow(RT_across_session)),
                                       sd(RT_across_session$RT[RT_across_session$type=='Congruent']) / sqrt(length(RT_across_session$RT[RT_across_session$type=='Congruent'])),
                                       sd(RT_across_session$RT[RT_across_session$type=='Incongruent']) / sqrt(length(RT_across_session$RT[RT_across_session$type=='Incongruent'])))
    )

    df_tmp_RT_avg <- cbind(df_tmp_RT_avg, session)
    
    RT_avg_across_session  <- rbind(RT_avg_across_session , df_tmp_RT_avg)
    
    # Load the accuracy
    fname4 <- paste(tmp_dir, "/Beh_data_accuracy_concatenated.csv")
    fname4 <- gsub(" ", "", fname4)
    
    ses_row <- data.frame(Session = session)
    df_acc <- readr::read_csv(fname4,id = "file_name")
    df_acc <- cbind(df_acc, session)
    df_acc <- subset(df_acc, select = -file_name)
    df_acc <- subset(df_acc, select = -participant)
    names(df_acc)[names(df_acc) == "key_resp.corr"] <- "acc"
    
    # Make it percent
    df_acc['acc']<-df_acc['acc']*100
    df_acc['sd']<-df_acc['sd']*100
    df_acc['SE']<-df_acc['SE']*100
    
    Acc_across_session  <- rbind(Acc_across_session , df_acc)
    
    # Calculate the average accuracy
    df_tmp_acc_avg <- data.frame(Type = c('Total', 'Congruent', 'Incongruent'),
                                 Acc_avg = c(mean(concatenated_df$key_resp.corr),
                                             mean(concatenated_df$key_resp.corr[concatenated_df$type=='congruent']),
                                             mean(concatenated_df$key_resp.corr[concatenated_df$type=='incongruent'])),
                                 SD = c(sd(concatenated_df$key_resp.corr),
                                        sd(concatenated_df$key_resp.corr[concatenated_df$type=='congruent']),
                                        sd(concatenated_df$key_resp.corr[concatenated_df$type=='incongruent'])),
                                 SE = c(sd(concatenated_df$key_resp.corr) / sqrt(length(concatenated_df$key_resp.corr)),
                                        sd(concatenated_df$key_resp.corr[concatenated_df$type=='congruent']) / sqrt(length(concatenated_df$key_resp.corr[concatenated_df$type=='congruent'])),
                                        sd(concatenated_df$key_resp.corr[concatenated_df$type=='incongruent']) 
                                        / sqrt(length(concatenated_df$key_resp.corr[concatenated_df$type=='incongruent'])))
                                 )
    
    df_tmp_acc_avg['Acc_avg'] <- df_tmp_acc_avg['Acc_avg'] *100
    df_tmp_acc_avg['SD'] <- df_tmp_acc_avg['SD'] *100
    df_tmp_acc_avg['SE'] <- df_tmp_acc_avg['SE'] *100
    
    df_tmp_acc_avg <- cbind(df_tmp_acc_avg, session)
    
    Acc_avg_across_session  <- rbind(Acc_avg_across_session , df_tmp_acc_avg)
    
    
    # Load the sensitivity
    fname5 <- paste(tmp_dir, "/Beh_data_sensitivity_tot_concatenated.csv")
    fname5 <- gsub(" ", "", fname5)
    
    
    ses_row <- data.frame(Session = session)
    df_sen <- readr::read_csv(fname5,id = "file_name")
    df_sen <- cbind(df_sen, session)
    df_sen <- subset(df_sen, select = -file_name)
    df_sen <- subset(df_sen, select = -participant)
    
    Sen_across_session  <- rbind(Sen_across_session , df_sen)
    
    # Calculate the average sensitivity
    df_tmp_sen_avg <- data.frame(Type = c('Total'),
                                 Sen_avg = c(mean(df_sen$d)),
                                 SD = c(sd(df_sen$d)),
                                 SE = c(sd(df_sen$d) / sqrt(length(df_sen$d))),
                                 Sen_bias_avg = c(sd(df_sen$c))
    )
    df_tmp_sen_avg <- cbind(df_tmp_sen_avg, session)
    
    Sen_avg_across_session  <- rbind(Sen_avg_across_session , df_tmp_sen_avg)
  }
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~           PLOTTING          ~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggridges)
fdir_plot <-  "/data/p_02825/cocoa/data/derivatives/stats/group/QC_Behavioral_across_Sessions/"
dir.create(paste(fdir_plot, "QC_full_retry"))


### RT
Ridge_RT <- ggplot(RT_across_session, aes(x = RT, y = session, fill=type, color=type)) +
  scale_fill_brewer(palette="Accent") +
  geom_density_ridges(alpha = .5, color = "white",
                      scale = 2.5, rel_min_height = .01) +
  labs(x = "Reaction time [s]", y = "Session", title="Reaction time across sessions") +
  guides(fill = "none") +
  theme_ridges() +
  guides(fill = guide_legend(title = "Stim. Type"))
ggsave(paste(fdir_plot, "/RidgePlot_AS_Reactiontime.png"))
print(Ridge_RT)

### Accuracy
Ridge_RT <- ggplot(Acc_across_session, aes(x = acc, y = session, fill=Type, color=Type)) +
  scale_fill_brewer(palette="Accent") +
  geom_density_ridges(alpha = .5, color = "white",
                      scale = 2.5, rel_min_height = .01) +
  labs(x = "Accuracy [%]", y = "Session", title="Accuracy across sessions") +
  guides(fill = "none") +
  theme_ridges() +
  geom_vline(xintercept = 50, linetype = "dotted", color = "black") +
  annotate(geom="text",x = 50, y = 0.5, label = "Chance level", size = 4, angle = 0) +
  xlim(0, 100) +
  guides(fill = guide_legend(title = "Stim. Type"))
ggsave(paste(fdir_plot, "/RidgePlot_AS_Accuracy.png"))
print(Ridge_RT)

### Sensitivity
Ridge_RT <- ggplot(Sen_across_session, aes(x = d, y = session)) +
  scale_fill_brewer(palette="Accent") +
  geom_density_ridges(alpha = .5, color = "white",
                      scale = 2.5, rel_min_height = .01) +
  labs(x = "Sensitivity", y = "Session", title="D Prime across sessions") +
  guides(fill = "none") +
  theme_ridges() +
  guides(fill = guide_legend(title = "Stim. Type"))
ggsave(paste(fdir_plot, "/RidgePlot_AS_Sensitivity.png"))
print(Ridge_RT)


# RT average
LinePLot <- ggplot(RT_avg_across_session, aes(x=session, y=RT_avg, group=Type, color=Type)) +
  geom_line(position=position_dodge(.9), size=1) +
  geom_point(position=position_dodge(.9), size=2)+
  geom_errorbar(aes(ymin=RT_avg-SD,
                    ymax=RT_avg+SD),width = 0,
                linetype = "dashed",
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=RT_avg-SE,
                    ymax=RT_avg+SE), width=.2,
                position=position_dodge(.9)) +
  scale_color_brewer(palette="Accent") +
  labs(y = "Reaction time [s]",x = "Session", title="Average reaction time across sessions")  + 
  ylim(0,3) +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.9, 0.9),
        plot.tag = element_text(size = 10))
ggsave(paste(fdir_plot, "/LinePlot_AS_RT_avg.png"))
print(LinePLot)

# Acc average
LinePLot <- ggplot(Acc_avg_across_session, aes(x=session, y=Acc_avg, group=Type, color=Type)) +
  geom_line(position=position_dodge(.9), size=1) +
  geom_point(position=position_dodge(.9), size=2) +
  geom_errorbar(aes(ymin=Acc_avg-SD,
                    ymax=Acc_avg+SD),width = 0,
                linetype = "dashed",
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=Acc_avg-SE,
                    ymax=Acc_avg+SE), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "black") +
  annotate(geom="text",x = 0.7, y = 52, label = "Chance level", size = 4, angle = 0) +
  scale_color_brewer(palette="Accent") +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.9, 0.9),
        plot.tag = element_text(size = 10) ) +
  labs(y = "Accuracy [%]",x = "Session", title="Average Accuracy across sessions") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, by = 10),expand = c(0, 0) )
ggsave(paste(fdir_plot, "/LinePlot_AS_Acc_avg.png"))
print(LinePLot)

# Sen average
LinePLot <- ggplot(Sen_avg_across_session, aes(x=session, y=Sen_avg, group=Sen_avg, color=Sen_avg)) +
  geom_line(position=position_dodge(.9), size=1) +
  geom_point(position=position_dodge(.9), size=2) +
  geom_errorbar(aes(ymin=Sen_avg-SD,
                    ymax=Sen_avg+SD),width = 0,
                linetype = "dashed",
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=Sen_avg-SE,
                    ymax=Sen_avg+SE), width=.2,
                position=position_dodge(.9)) +
  labs(y = "Sensitivity [dprime]",x = "Session", title="Average D Prime across sessions") +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.95, 0.9),
        plot.tag = element_text(size = 10) ) +
  ylim(0,3.5)
ggsave(paste(fdir_plot, "/LinePlot_AS_Sen_avg.png"))
print(LinePLot)



