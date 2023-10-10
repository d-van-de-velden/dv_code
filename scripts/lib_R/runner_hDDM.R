

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
library(hBayesDM)


# ~~~~~~~~~~~ Get and process the data for all sessions

minimum_reaction_time <- 0.3 # in seconds
iterations <- 4000           # Iterations for HDDM

sessions = c('ses-1', 'ses-2', 'ses-3', 'ses-4')
fdir_get_across <-
  "/data/p_02825/cocoa/data/derivatives/stats/group/"
for (session in sessions) {
  print(session)
  
  
  fdir = paste(fdir_get_across, session, '/beh/')
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
    
    # transform stim_images to character
    colnames(concatenated_df)[colnames(concatenated_df) == "stim_images_order"] <-
      "character"
    
    
    
    
    
    #~~~~~~~~~~ Drift diffusion model ~~~~~~~~~~~~~~~~~~
    
    
    
    # create data frame for DDM
    DDMdf_tot <- concatenated_df
    colnames(DDMdf_tot)[colnames(DDMdf_tot) == "key_resp.corr"] <-
      "choice" # change the name of columns
    colnames(DDMdf_tot)[colnames(DDMdf_tot) == "participant"] <-
      "subjID"
    colnames(DDMdf_tot)[colnames(DDMdf_tot) == "key_resp.rt"] <- "RT"
    
    # remove outlier RT < 300ms
    outlier_idx <- which(DDMdf_tot$RT > minimum_reaction_time)
    outlier_idx <- outlier_idx[1:length(DDMdf_tot$RT)]
    DDMdf_tot$RT <- DDMdf_tot$RT[outlier_idx]
    
    # specify choice number
    DDMdf_tot$choice <- gsub("1", "2", DDMdf_tot$choice) # 2 = correct
    DDMdf_tot$choice <- gsub("0", "1", DDMdf_tot$choice) # 1 = incorrect
    
    DDMdf_cong   <- DDMdf_tot
    DDMdf_incong <- DDMdf_tot
    
    DDMdf_cong   <- DDMdf_cong[DDMdf_cong$type == 'congruent', ]
    DDMdf_incong <- DDMdf_incong[DDMdf_incong$type == 'incongruent', ]
    
    
    # ~~~~~~~~~ Run the model with a given data.frame (DDMdf)   TOTAL
    dir.create(gsub(" ", "", paste(fdir, "hDDM_csv/") ))
    
    output <- choiceRT_ddm(
      data = DDMdf_tot,
      niter = iterations,
      nwarmup = iterations * 0.25,
      nchain = 4,
      ncore = 4
    )
    # Get each parameter
    df_eachparameter <- output$parVals
    TOT_save_dv <- data.frame(
      Alpha = c(df_eachparameter$alpha),
      Beta = c(df_eachparameter$beta),
      Delta = c(df_eachparameter$delta),
      Tau = c(df_eachparameter$tau),
      Type = 'Total'
    )
    fname <-
      gsub(" ", "", paste(fdir, "hDDM_csv/hDDM_results_total.csv"))
    write.csv(TOT_save_dv, fname)
    
    
    output_cong <- choiceRT_ddm(
      data = DDMdf_cong,
      niter = iterations,
      nwarmup = iterations * 0.25,
      nchain = 4,
      ncore = 4
    )
    # Get each parameter
    df_eachparameter_cong <- output_cong$parVals
    CONG_save_dv <- data.frame(
      Alpha = c(df_eachparameter_cong$alpha),
      Beta = c(df_eachparameter_cong$beta),
      Delta = c(df_eachparameter_cong$delta),
      Tau = c(df_eachparameter_cong$tau),
      Type = 'Congruent'
    )
    fname <-
      gsub(" ",
           "",
           paste(fdir, "hDDM_csv/hDDM_results_congruent.csv"))
    write.csv(CONG_save_dv, fname)
    
    output_incong <- choiceRT_ddm(
      data = DDMdf_incong,
      niter = iterations,
      nwarmup = iterations * 0.25,
      nchain = 4,
      ncore = 4
    )
    # Get each parameter
    df_eachparameter_incong <- output_incong$parVals
    INCONG_save_dv <-
      data.frame(
        Alpha = c(df_eachparameter_incong$alpha),
        Beta = c(df_eachparameter_incong$beta),
        Delta = c(df_eachparameter_incong$delta),
        Tau = c(df_eachparameter_incong$tau),
        Type = 'Incongruent'
      )
    fname <-
      gsub(" ",
           "",
           paste(fdir, "hDDM_csv/hDDM_results_incongruent.csv"))
    write.csv(INCONG_save_dv, fname)
    
    
    
    # ~~~~~~~~~~ visualizaton of DDM
    # (5.1) Visually check convergence of the sampling chains (should look like 'hairy caterpillars')
    sanity_4000 <- plot(output, type = "trace")
    
    # check the mean of each subject's parameters
    output$allIndPars
    
    # plot each parameter with multiple bar plots in one figure
    library(gridExtra)
    
    # plot the each parameter
    plot1 <-
      ggplot(output$allIndPars, aes(x = subjID, y = alpha, fill = subjID)) +
      geom_bar(stat = "identity",
               width = 0.7,
               position = position_dodge()) +
      labs(title = "Alpha per participant", x = "Participants", y = "Alpha [a.U.]") +
      theme_classic()
    plot1 + theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
    plot1 <- plot1 + guides(fill = "none")
    ggsave(paste(fdir, "QC/hBayesDDModeling_Alpha_per_participant.png"))
    
    plot2 <-
      ggplot(output$allIndPars, aes(x = subjID, y = beta, fill = subjID)) +
      geom_bar(stat = "identity",
               width = 0.7,
               position = position_dodge()) +
      labs(title = "Beta per participant", x = "Participants", y = "Beta [a.U.]") +
      theme_classic()
    plot2 + theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
    plot2 <- plot2 + guides(fill = "none")
    ggsave(paste(fdir, "QC/hBayesDDModeling_Beta_per_participant.png"))
    
    plot3 <-
      ggplot(output$allIndPars, aes(x = subjID, y = delta, fill = subjID)) +
      geom_bar(stat = "identity",
               width = 0.7,
               position = position_dodge()) +
      labs(title = "Delta per participant", x = "Participants", y = "Delta [a.U.]") +
      theme_classic()
    plot3 + theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
    plot3 <- plot3 + guides(fill = "none")
    ggsave(paste(fdir, "QC/hBayesDDModeling_Delta_per_participant.png"))
    
    plot4 <-
      ggplot(output$allIndPars, aes(x = subjID, y = tau, fill = subjID)) +
      geom_bar(stat = "identity",
               width = 0.7,
               position = position_dodge()) +
      labs(title = "Tau per participant", x = "Participants", y = "Tau [s]") +
      theme_classic()
    plot4 + theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
    plot4 <- plot4 + guides(fill = "none")
    ggsave(paste(fdir, "QC/hBayesDDModeling_Tau_per_participant.png"))
    
    # concatenate all the plots
    png(
      file = paste(fdir, "QC/hBayesDDModeling_per_participant.png"),
      width = 1800,
      height = 1800
    )
    combined_plot <-
      grid.arrange(plot1,
                   plot2,
                   plot3,
                   plot4,
                   nrow = 2,
                   ncol = 2)
    # Display the concatenated plot
    
    dev.off()
    print(combined_plot)
    
    
    SubName <- concatenated_df$participant
    SubLists <- unique(SubName)
    
    
    # alpha data frame
    df_alpha <- data.frame(as.vector(df_eachparameter$alpha), SubLists)
    colnames(df_alpha) <- c("alpha", "participant")
    
    # beta data frame
    df_beta <- data.frame(as.vector(df_eachparameter$beta), SubLists)
    colnames(df_beta) <- c("beta", "participant")
    
    # delta data frame
    df_delta <- data.frame(as.vector(df_eachparameter$delta), SubLists)
    colnames(df_delta) <- c("delta", "participant")
    
    # tau data frame
    df_tau <- data.frame(as.vector(df_eachparameter$tau), SubLists)
    colnames(df_tau) <- c("tau", "participant")
    
    
    
    
    # ~~~~~~~~~ plot alpha
    alpha_4000 <- df_alpha %>%
      ggplot(aes(x = participant, y = alpha, fill = participant)) +
      
      # add half-violin from {ggdist} package
      stat_halfeye(
        # adjust bandwidth
        adjust = 0.5,
        # move to the right
        justification = -0.2,
        # remove the slub interval
        .width = 0,
        point_colour = NA
      ) +
      geom_boxplot(width = 0.2,
                   # removing outliers
                   outlier.color = NA,
                   alpha = 0.5) +
      labs(title = "4000\nResponse caution (alpha)", x = "Participant", y = "Alpha [a.U.]") +
      scale_y_continuous(
        limits = c(0, 4.5),
        breaks = seq(0, 6, by = 1),
        expand = c(0, 0)
      )
    alpha_4000 <- alpha_4000 + guides(fill = "none")
    alpha_4000 + theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
    ggsave(paste(fdir, "QC/hBayesDDModeling_Response caution_alpha.png"))
    
    
    # ~~~~~~~~ (5.2) plot group-average of alpha
    df_averageAlpha <- data.frame(
      x_axis = c('1'),
      Mean = mean(df_alpha$alpha),
      SD = sd(df_alpha$alpha)
    )
    df_averageAlpha$SE <- df_averageAlpha$SD / sqrt(nrow(df_alpha))
    
    plot_avgAlpha <- ggplot(df_averageAlpha, aes(x = x_axis, y = Mean)) +
      geom_bar(
        stat = "identity",
        width = 0.5,
        position = position_dodge(),
        fill = "transparent",
        color = "black"
      ) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                    width = 0,
                    position = position_dodge(.9)) +
      geom_errorbar(
        aes(ymin = Mean - SD, ymax = Mean + SD),
        width = 0,
        linetype = "dashed",
        position = position_dodge(.9)
      ) +
      labs(title = "Group-average response caution (alpha)", x = "Sessions", y = "Alpha [a.U.]") +
      theme_classic() +
      labs(tag = "──  SE\n- - -  SD") +
      theme(plot.tag.position = c(0.85, 0.9),
            plot.tag = element_text(size = 10)) +
      scale_y_continuous(
        limits = c(0, 3),
        breaks = seq(0, 3, by = 0.5),
        expand = c(0, 0)
      )
    ggsave(paste(
      fdir,
      "QC/hBayesDDModeling_Group_average_response_caution_alpha.png"
    ))
    print(plot_avgAlpha)
    
    
    
    # ~~~~~~~~~~ plot beta
    beta_4000 <- df_beta %>%
      ggplot(aes(x = participant, y = beta, fill = participant)) +
      
      # add half-violin from {ggdist} package
      stat_halfeye(
        # adjust bandwidth
        adjust = 0.5,
        # move to the right
        justification = -0.2,
        # remove the slub interval
        .width = 0,
        point_colour = NA
      ) +
      geom_boxplot(width = 0.2,
                   # removing outliers
                   outlier.color = NA,
                   alpha = 0.5) +
      labs(title = " \nResponse bias (beta)", x = "Participant", y = "Beta [a.U.]") +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = c(0, 0)
      ) +
      annotate(
        geom = "text",
        x = 1,
        y = 0.03,
        label = "incorrect",
        size = 5,
        angle = 0
      ) +
      annotate(
        geom = "text",
        x = 1,
        y = 0.97,
        label = "correct",
        size = 5,
        angle = 0
      )
    beta_4000 <- beta_4000 + guides(fill = "none")
    beta_4000 + theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
    ggsave(paste(fdir, "QC/hBayesDDModeling_Response bias_beta.png"))
    
    # ~~~~~~~~~~ (5.2) plot average of beta
    df_averageBeta <- data.frame(
      x_axis = c('1'),
      Mean = mean(df_beta$beta),
      SD = sd(df_beta$beta)
    )
    df_averageBeta$SE <- df_averageBeta$SD / sqrt(nrow(df_beta))
    
    plot_avgBeta <- ggplot(df_averageBeta, aes(x = x_axis, y = Mean)) +
      geom_bar(
        stat = "identity",
        width = 0.5,
        position = position_dodge(),
        fill = "transparent",
        color = "black"
      ) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                    width = 0,
                    position = position_dodge(.9)) +
      geom_errorbar(
        aes(ymin = Mean - SD, ymax = Mean + SD),
        width = 0,
        linetype = "dashed",
        position = position_dodge(.9)
      ) +
      labs(title = "Group-average response bias (beta)", x = "Sessions", y = "Beta [a.U.]") +
      theme_classic() +
      labs(tag = "──  SE\n- - -  SD") +
      theme(plot.tag.position = c(0.85, 0.9),
            
            plot.tag = element_text(size = 10)) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = c(0, 0)
      )
    ggsave(paste(
      fdir,
      "QC/hBayesDDModeling_Group_average_Response bias_beta.png"
    ))
    print(plot_avgBeta)
    
    
    
    # ~~~~~~~~~~ plot delta
    delta_4000 <- df_delta %>%
      ggplot(aes(x = participant, y = delta, fill = participant)) +
      
      # add half-violin from {ggdist} package
      stat_halfeye(
        # adjust bandwidth
        adjust = 0.5,
        # move to the right
        justification = -0.2,
        # remove the slub interval
        .width = 0,
        point_colour = NA
      ) +
      geom_boxplot(width = 0.2,
                   # removing outliers
                   outlier.color = NA,
                   alpha = 0.5) +
      labs(title = " \nSpeed of evidence accumulation (delta)", x = "Participant", y = "Delta [a.U.]") +
      scale_y_continuous(
        limits = c(0, 3),
        breaks = seq(0, 6, by = 1),
        expand = c(0, 0)
      )
    delta_4000 <- delta_4000 + guides(fill = "none")
    delta_4000 + theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
    ggsave(paste(
      fdir,
      "QC/hBayesDDModeling_Speed_of_evidence_accumulation_delta.png"
    ))
    
    # ~~~~~~~~~~~~ (5.2) plot average delta
    df_averageDelta <- data.frame(
      x_axis = c('1'),
      Mean = mean(df_delta$delta),
      SD = sd(df_delta$delta)
    )
    
    df_averageDelta$SE <- df_averageDelta$SD / sqrt(nrow(df_delta))
    
    plot_avgDelta <- ggplot(df_averageDelta, aes(x = x_axis, y = Mean)) +
      geom_bar(
        stat = "identity",
        width = 0.5,
        position = position_dodge(),
        fill = "transparent",
        color = "black"
      ) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                    width = 0,
                    position = position_dodge(.9)) +
      geom_errorbar(
        aes(ymin = Mean - SD, ymax = Mean + SD),
        width = 0,
        linetype = "dashed",
        position = position_dodge(.9)
      ) +
      labs(title = "Group-average speed of evidence accumulation (delta)", x =
             "Sessions", y = "Delta [a.U.]") +
      theme_classic() +
      labs(tag = "──  SE\n- - -  SD") +
      theme(plot.tag.position = c(0.85, 0.9),
            plot.tag = element_text(size = 10)) +
      scale_y_continuous(
        limits = c(0, 4),
        breaks = seq(0, 4, by = 0.5),
        expand = c(0, 0)
      )
    ggsave(
      paste(
        fdir,
        "QC/hBayesDDModeling_Group_average_Speed_of_evidence_accumulation_delta.png"
      )
    )
    print(plot_avgDelta)
    
    
    # ~~~~~~~~~~~~ plot tau
    tau_4000 <- df_tau %>%
      ggplot(aes(x = participant, y = tau, fill = participant)) +
      
      # add half-violin from {ggdist} package
      stat_halfeye(
        # adjust bandwidth
        adjust = 0.5,
        # move to the right
        justification = -0.2,
        # remove the slub interval
        .width = 0,
        point_colour = NA
      ) +
      geom_boxplot(width = 0.2,
                   outlier.color = NA,
                   alpha = 0.5) +
      labs(title = " \nSensory encoding and motor execution", x = "Participant", y = "Tau [s]") +
      scale_y_continuous(
        limits = c(0, 1.75),
        breaks = seq(0, 1.75, by = 0.1),
        expand = c(0, 0)
      )
    tau_4000 <- tau_4000 + guides(fill = "none")
    tau_4000 + theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
    ggsave(
      paste(
        fdir,
        "QC/hBayesDDModeling_Sensory_encoding_and_motor_execution_Tau.png"
      )
    )
    
    
    # ~~~~~~~~~~~ (5.2) plot average tau
    df_averageTau <- data.frame(
      x_axis = c('1'),
      Mean = mean(df_tau$tau),
      SD = sd(df_tau$tau)
    )
    df_averageTau$SE <- df_averageTau$SD / sqrt(nrow(df_tau))
    a <- df_averageTau$Mean - df_averageTau$SE
    b <- df_averageTau$Mean + df_averageTau$SE
    c <- df_averageTau$Mean - df_averageTau$SD
    d <- df_averageTau$Mean + df_averageTau$SD
    
    plot_avgTau <-
      ggplot(df_averageTau, aes(x = x_axis, y = df_averageTau$Mean)) +
      geom_bar(
        stat = "identity",
        width = 0.5,
        position = position_dodge(),
        fill = "transparent",
        color = "black"
      ) +
      geom_errorbar(aes(ymin = a, ymax = b),
                    width = 0,
                    position = position_dodge(.9)) +
      geom_errorbar(
        aes(ymin = c, ymax = d),
        width = 0,
        linetype = "dashed",
        position = position_dodge(.9)
      ) +
      labs(title = "Group-average sensory encoding and motor execution (tau)", x =
             "Sessions", y = "Tau [s]") +
      theme_classic() + labs(tag = "──  SE\n- - -  SD") +
      theme(plot.tag.position = c(0.85, 0.9),
            plot.tag = element_text(size = 10)) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = c(0, 0)
      )
    ggsave(
      paste(
        fdir,
        "QC/hBayesDDModeling_Group_average_sensory_encoding_and_motor_execution_Tau.png"
      )
    )
    print(plot_avgTau)
    
    
    
    # ~~~~~~~~ Plot together
    # concatenate all the plots
    png(
      file = paste(fdir, "QC/hBayesDDModeling_Parameters_per_Subject.png"),
      width = 1800,
      height = 1800
    )
    
    
    # concatenate all the parameter for each participant
    combined_plot <-
      grid.arrange(alpha_4000,
                   beta_4000,
                   delta_4000,
                   tau_4000,
                   nrow = 1,
                   ncol = 4)
    dev.off()
    print(combined_plot)
    
    # concatenate all the average plots
    png(
      file = paste(fdir, "QC/hBayesDDModeling_Group_average.png"),
      width = 1800,
      height = 1200
    )
    combined_plot_Avg <-
      grid.arrange(
        plot_avgAlpha,
        plot_avgBeta,
        plot_avgDelta,
        plot_avgTau,
        nrow = 1,
        ncol = 4
      )
    dev.off()
    print(combined_plot_Avg)
    
    
  }
}
  



# ~~~~~~~~ RIDGE plot group-average of alpha


dr_across_session <- data.frame(Alpha = c(),
                                Beta = c(),
                                Delta= c(),
                                Tau= c(),
                                Session=c(),
                                Type=c())
for (session in sessions){
  print(session)
  tmp_dir = paste(fdir_get_across,session,"/beh/hDDM_csv")
  tmp_dir <- gsub(" ", "", tmp_dir)
  if (dir.exists(tmp_dir)){
  
    fname1 <- paste(tmp_dir, "/hDDM_results_total.csv")
    fname1 <- gsub(" ", "", fname1)
    fname2 <- paste(tmp_dir, "/hDDM_results_congruent.csv")
    fname2 <- gsub(" ", "", fname2)
    fname3 <- paste(tmp_dir, "/hDDM_results_incongruent.csv")
    fname3 <- gsub(" ", "", fname3)
    
    ses_row <- data.frame(Session = session)
    df1 <- readr::read_csv(fname1,id = "file_name")
    df1 <- cbind(df1, session)
    df2 <- readr::read_csv(fname2, id = "file_name")
    df2 <- cbind(df2, session)
    df3 <- readr::read_csv(fname3, id = "file_name")
    df3 <- cbind(df3, session)
    
    dr_across_session <- rbind(dr_across_session, df1, df2, df3 )
  }
}

fdir <- '/data/p_02825/cocoa/data/derivatives/stats/group/QC_Behavioral_across_Sessions/ QC_full'
library(ggridges)


Alpha_plot <- ggplot(dr_across_session, aes(x = Alpha, y = session, fill=Type, color=Type)) +
  scale_fill_brewer(palette="Accent") +
  geom_density_ridges(alpha = .8, color = "white",
                      scale = 2.5, rel_min_height = .01) +
  labs(x = "Alpha [a.U.]", y = "Session") +
  guides(fill = "none") +
  theme_ridges() +
  guides(fill = guide_legend(title = "Stim. Type"))
ggsave(paste(fdir, "QC/hBayesDDModeling_Across_session_Alpha.png"))
print(Alpha_plot)


Beta_plot <- ggplot(dr_across_session, aes(x = Beta, y = session, fill=Type)) +
  scale_fill_brewer(palette="Accent") +
  geom_density_ridges(alpha = .8, color = "white",
                      scale = 2.5, rel_min_height = .01) +
  labs(x = "Beta [a.U.]", y = "Session") +
  guides(fill = "none") +
  theme_ridges() +
  guides(fill = guide_legend(title = "Stim. Type")) +
  scale_x_continuous( limits = c(0.2, .8), breaks = seq(0.2, .8, by = 0.1), expand = c(0, 0))
ggsave(paste(fdir, "QC/hBayesDDModeling_Across_session_Betaa.png"))
print(Beta_plot)


Delta_plot <- ggplot(dr_across_session, aes(x = Delta, y = session, fill=Type)) +
  scale_fill_brewer(palette="Accent") +
  geom_density_ridges(alpha = .8, color = "white",
                      scale = 2.5, rel_min_height = .01) +
  labs(x = "Delta [a.U.]", y = "Session") +
  guides(fill = "none") +
  theme_ridges() +
  guides(fill = guide_legend(title = "Stim. Type"))
ggsave(paste(fdir, "QC/hBayesDDModeling_Across_session_Delta.png"))
print(Delta_plot)



Tau_plot <- ggplot(dr_across_session, aes(x = Tau, y = session, fill=Type)) +
  scale_fill_brewer(palette="Accent") +
  geom_density_ridges(alpha = .8, color = "white",
                      scale = 2.5, rel_min_height = .01) +
  labs(x = "Tau [s]", y = "Session") +
  guides(fill = "none") +
  theme_ridges() +
  guides(fill = guide_legend(title = "Stim. Type")) +
  scale_x_continuous( limits = c(.0, 2.0), breaks = seq(.0, 2.0, by = 0.2), expand = c(0, 0))
ggsave(paste(fdir, "QC/hBayesDDModeling_Across_session_Tau.png"))
print(Tau_plot)

