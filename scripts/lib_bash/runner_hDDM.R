# Maxplanck_BDA

# cwd <- '/Users/muku/Desktop/Maxplanck_BDA/CSV_data/'

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
fdir <- getwd()
list_of_files <- list.files(path = fdir,
                            recursive = TRUE,
                            pattern = "\\.csv$",
                            full.names = TRUE)


concatenated_df <- c()
for (i in 1:length(list_of_files)){
  variable_names <- c("tmp_var")
  df <- readr::read_csv(list_of_files[i], id = "file_name")
  assign(variable_names[1],df)
  print(df)
  
  # ~~~~~~~~~~~ pick up the correct, reaction time, participant name, letter, stimulus stype
  tmp_var_relevant = tmp_var[c("key_resp.corr", "key_resp.rt", "participant","stim_images_order","type")]
  tmp_var_relevant = tmp_var_relevant[complete.cases(tmp_var_relevant), ] # remove NaN
  
  
  # combine all the subjects
  concatenated_df <- rbind(concatenated_df, tmp_var_relevant)
}


concatenated_df$stim_images_order <- gsub("images/alphabets/B.png", "B", concatenated_df$stim_images_order)
concatenated_df$stim_images_order <- gsub("images/alphabets/F.png", "F", concatenated_df$stim_images_order)
concatenated_df$stim_images_order <- gsub("images/alphabets/K.png", "K", concatenated_df$stim_images_order)
concatenated_df$stim_images_order <- gsub("images/alphabets/M.png", "M", concatenated_df$stim_images_order)
concatenated_df$stim_images_order <- gsub("images/alphabets/P.png", "P", concatenated_df$stim_images_order)
concatenated_df$stim_images_order <- gsub("images/alphabets/T.png", "T", concatenated_df$stim_images_order)

# transform stim_images to character
colnames(concatenated_df)[colnames(concatenated_df) == "stim_images_order"] <- "character"





# ~~~~~~~~~~~ plot an average accuracy for each participant or each letter ~~~~~~~~~~~~~~~~~~
# summary
# Define the data_summary function
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#summary result
df1 <- data_summary(concatenated_df, varname="key_resp.corr", groupnames=c("participant"))

# standard error
N_trial <- c(nrow(concatenated_df[concatenated_df$participant=='anonymous',]),
             nrow(concatenated_df[concatenated_df$participant=='anonymous2',]),
             nrow(concatenated_df[concatenated_df$participant=='EDX3',]),
             nrow(concatenated_df[concatenated_df$participant=='EVY3',]))

df1$SE <- df1$sd/sqrt(N_trial)

# chance level
chance_level <- 50


# plot average for each participant

p1<- ggplot(df1, aes(x=participant, y=key_resp.corr*100, fill = participant)) + 
  geom_bar(stat="identity",
           width = 0.7, position=position_dodge()) +
  geom_errorbar(aes(ymin=key_resp.corr*100-SE*100, ymax=key_resp.corr*100+SE*100), width=0,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=key_resp.corr*100-sd*100, ymax=100), width=0, linetype = "dashed",
                position=position_dodge(.9)) + 
  labs(title="Accuracy per participant", x="Participants", y = "Response Accuracy (%)") +
  theme_classic() +
  geom_hline(yintercept = chance_level, linetype = "dotted", color = "black") +
  ylim(0, 100) +
  annotate(geom="text",x = 4.5, y = 48, label = "Chance level", size = 4, angle = 90) +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.85, 0.9),
        plot.tag = element_text(size = 10) )
p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste(fdir, "QC/Accuracy_per_participant.png"))
print(p1)


# ~~~~~~~~~~~~ (1.2) plot average for each letter
# summarys
df2 <- data_summary(concatenated_df, varname="key_resp.corr", 
                    groupnames=c("character"))

# standard error
N_trial <- c(nrow(concatenated_df[concatenated_df$character=='B',]),
             nrow(concatenated_df[concatenated_df$character=='F',]),
             nrow(concatenated_df[concatenated_df$character=='K',]),
             nrow(concatenated_df[concatenated_df$character=='M',]),
             nrow(concatenated_df[concatenated_df$character=='P',]),
             nrow(concatenated_df[concatenated_df$character=='T',]))

df2$SE <- df2$sd/sqrt(N_trial)
df2$letter <- df2$character

# chance level
chance_level <- 50

# plot

p2 <- ggplot(df2, aes(x=character, y=key_resp.corr*100, fill = letter)) + 
    geom_bar(stat="identity",
             width = 0.5, position=position_dodge()) +
    geom_errorbar(aes(ymin=key_resp.corr*100-SE*100, ymax=key_resp.corr*100+SE*100), width=0,
                  position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=key_resp.corr*100-sd*100, ymax=100), width=0, linetype = "dashed",
                  position=position_dodge(.9)) + 
    labs(title="Group-average accuracy for each letter", x="Letter", y = "Response Accuracy (%)") +
    theme_classic() +
    geom_hline(yintercept = chance_level, linetype = "dotted", color = "black") +
    ylim(0, 100) +
    annotate(geom="text",x = 6.5, y = 48, label = "Chance level", size = 4, angle = 90) +
    labs(tag = "──  SE\n- - -  SD") +
    theme(plot.tag.position = c(0.95, 0.9),
          plot.tag = element_text(size = 10) ) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100, by = 10),expand = c(0,0))
ggsave(paste(fdir, "QC/Group_average_accuracy_for_each_letter.png"))
print(p2)


# ~~~~~~~~~~~~ (1.1) Group average accuracy, one bar plot
df3 <- data.frame(x_axis = c('1'),key_resp.corr = mean(concatenated_df$key_resp.corr), sd = sd(concatenated_df$key_resp.corr))
df3$SE <- df3$sd/sqrt(nrow(concatenated_df))

# chance level
chance_level <- 50

# plot
p3 <- ggplot(df3, aes(x = x_axis, y=key_resp.corr*100)) + 
  geom_bar(stat="identity",
           width = 0.5, position=position_dodge(), fill="transparent", color="black") +
  geom_errorbar(aes(ymin=key_resp.corr*100-SE*100, ymax=key_resp.corr*100+SE*100), width=0,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=key_resp.corr*100-sd*100, ymax=100), width=0, linetype = "dashed",
                position=position_dodge(.9)) + 
  labs(title="Group-average accuracy", x="Sessions", y = "Response Accuracy (%)") +
  theme_classic() +
  geom_hline(yintercept = chance_level, linetype = "dotted", color = "black") +
  ylim(0, 100) +
  annotate(geom="text",x = 1.4, y = 48, label = "Chance level", size = 4, angle = 90) +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.85, 0.9),
        plot.tag = element_text(size = 10) ) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, by = 10),expand = c(0,0))
ggsave(paste(fdir, "QC/Group_average_accuracy_one_bar_plot.png"))
print(p3)






# ~~~~~~~~~~~~ reaction time ~~~~~~~~~~~~~~~~~~

concatenated_df %>% 
  ggplot(aes(x = participant, y = key_resp.rt*1000, fill = participant)) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.5,
    # move to the right
    justification = -0.2,
    # remove the slub interval
    .width = 0,
    point_colour = NA) +
  geom_boxplot(
    width = 0.2,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  stat_dots(
    # ploting on left side
    side = "left",
    # adjusting position
    justification = 1.2,
    # adjust grouping (binning) of observations
    binwidth = 40
  ) +
  labs(
    title="Reaction time per participant", x="Participant", y = "Reaction time (ms)") 
ggsave(paste(fdir, "QC/Reaction_time_per_participant.png"))


# ~~~~~~~~~~~~~ (2.2) Reaction time for each character
reactiontime_df <- concatenated_df
reactiontime_df$letter <- reactiontime_df$character

#reactiontime_df %>% 

RT_plot_letter <- ggplot(reactiontime_df,aes(x = character, y = key_resp.rt*1000, fill = letter)) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.3,
    # move to the right
    justification = -0.2,
    # remove the slub interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.2,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  labs(title="Group-average reaction time for each letter", x="Letter", y = "Reaction time (ms)") 
ggsave(paste(fdir, "QC/Group_average_reaction_time_for_each_letter.png"))

# calculate the values of the box plot
# install function
ggplot2_boxplot <- function(x){
  quartiles <- as.numeric(quantile(x,
                                   probs = c(0.25, 0.5, 0.75)))
  
  names(quartiles) <- c("25th percentile",
                        "50th percentile\n(median)",
                        "75th percentile")
  
  IQR <- diff(quartiles[c(1,3)])
  upper_whisker <- max(x[x < (quartiles[3] + 1.5 * IQR)])
  lower_whisker <- min(x[x > (quartiles[1] - 1.5 * IQR)])
  
  upper_dots <- x[x > (quartiles[3] + 1.5*IQR)]
  lower_dots <- x[x < (quartiles[1] - 1.5*IQR)]
  
  return(list("quartiles" = quartiles,
              "25th percentile" = as.numeric(quartiles[1]),
              "50th percentile\n(median)" = as.numeric(quartiles[2]),
              "75th percentile" = as.numeric(quartiles[3]),
              "IQR" = IQR,
              "upper_whisker" = upper_whisker,
              "lower_whisker" = lower_whisker,
              "upper_dots" = upper_dots,
              "lower_dots" = lower_dots))
}
ggsave(paste(fdir, "QC/saving_plot6.png"))

# calculate Min, Q1, Median, Q3, Max
ggplot_output_B <- ggplot2_boxplot(reactiontime_df$key_resp.rt[reactiontime_df$letter=='B'])
ggsave(paste(fdir, "QC/saving_plot_B.png"))
ggplot_output_F <- ggplot2_boxplot(reactiontime_df$key_resp.rt[reactiontime_df$letter=='F'])
ggsave(paste(fdir, "QC/saving_plot_F.png"))
ggplot_output_K <- ggplot2_boxplot(reactiontime_df$key_resp.rt[reactiontime_df$letter=='K'])
ggsave(paste(fdir, "QC/saving_plot_K.png"))
ggplot_output_M <- ggplot2_boxplot(reactiontime_df$key_resp.rt[reactiontime_df$letter=='M'])
ggsave(paste(fdir, "QC/saving_plot_M.png"))
ggplot_output_P <- ggplot2_boxplot(reactiontime_df$key_resp.rt[reactiontime_df$letter=='P'])
ggsave(paste(fdir, "QC/saving_plot_P.png"))
ggplot_output_T <- ggplot2_boxplot(reactiontime_df$key_resp.rt[reactiontime_df$letter=='T'])
ggsave(paste(fdir, "QC/saving_plot_T.png"))




# ~~~~~~~~~~~ (2.1) Group-average reaction time
df_averageRT <- data.frame(x_axis = c('1'),Mean_RT = mean(concatenated_df$key_resp.rt), SD = sd(concatenated_df$key_resp.rt))
df_averageRT$SE <- df_averageRT$SD/sqrt(nrow(concatenated_df))

Mean_RT<-df_averageRT$Mean_RT
SE<-df_averageRT$SE
SD<-df_averageRT$SD

# plot
plot_avgRT <- ggplot(df_averageRT, aes(x=x_axis, y=Mean_RT*1000)) + 
  geom_bar(stat="identity", width = 0.5, position=position_dodge(),
           fill="transparent", color="black") +
  geom_errorbar(aes(ymin=Mean_RT*1000-SE*1000, ymax=Mean_RT*1000+SE*1000), width=0,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=Mean_RT*1000-SD*1000, ymax=Mean_RT*1000+SD*1000), width=0, 
                linetype = "dashed", position=position_dodge(.9)) + 
  labs(title="Group-average reaction time", x="Sessions", y = "Reaction time (ms)") +
  theme_classic() +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.85, 0.9),
        plot.tag = element_text(size = 10) ) +
  scale_y_continuous(limits = c(0,2500), breaks = seq(0,2500, by = 500),expand = c(0,0))
ggsave(paste(fdir, "QC/Group_average_reaction_time.png"))
print(plot_avgRT)




# ~~~~~~~~~~~~~ (3.1) Integrated speed-accuracy plot for each letter ~~~~~~~~~~~~~~~~~~
# accuracy over a period of time

rt_range <- 3.5
n_bins <- 10
break_seq <- seq(0, rt_range, rt_range/n_bins)

timeslice_range <- concatenated_df %>%
  mutate(RT_bin = cut(key_resp.rt, breaks = break_seq)) %>%
  group_by(RT_bin, character) %>%
  mutate(RT_bin_avg = mean(key_resp.rt, na.rm = T))

count_range <- timeslice_range %>%
  group_by(RT_bin, character) %>%
  summarise(subjcount = n_distinct(participant),totalcount=n())

timeslice_range <- timeslice_range %>%
  group_by(RT_bin_avg, character, participant) %>% 
  summarise(ss_acc = mean(key_resp.corr, na.rm=T)) %>% 
  group_by(RT_bin_avg, character) %>%
  summarise(mean = mean(ss_acc),n = n())
timeslice_range$n <- count_range$totalcount
timeslice_range$letter <- timeslice_range$character

ggplot(aes(x=RT_bin_avg*1000, y=mean, weight = n, 
           color = letter, group = letter), 
       data = timeslice_range) + 
  geom_point(aes(size = n)) +
  geom_smooth(method = "lm", formula = y ~ poly(x,2), se = FALSE) +
  geom_hline(yintercept = 0.5, lty = "dashed") +
  xlab("Average RT (ms)") +
  ylab("Proportion Correct") +
  labs(title="Speed-accuracy distribution plot for each letter") +
  ylim(0,1)
ggsave(paste(fdir, "QC/Speed_accuracy_distribution_plot_for_each_letter.png"))


# ~~~~~~~~~~~~~ integrated speed accuracy plot for each participant
# accuracy over a period of time

rt_range <- 3.5
n_bins <- 10
break_seq <- seq(0, rt_range, rt_range/n_bins)

timeslice_range <- concatenated_df %>%
  mutate(RT_bin = cut(key_resp.rt, breaks = break_seq)) %>%
  group_by(RT_bin, participant) %>%
  mutate(RT_bin_avg = mean(key_resp.rt, na.rm = T))

count_range <- timeslice_range %>%
  group_by(RT_bin, participant) 

timeslice_range <- timeslice_range %>%
  group_by(RT_bin_avg, character, participant) %>% 
  summarise(ss_acc = mean(key_resp.corr, na.rm=T)) %>% 
  group_by(RT_bin_avg, participant) %>%
  summarise(mean = mean(ss_acc),
            n = n())

# count the total number of bin
n <- c()
for (i in 1:length(timeslice_range$RT_bin_avg)){
  N <- sum(count_range$RT_bin_avg == timeslice_range$RT_bin_avg[i])
  n <- c(n,N)
}
timeslice_range$n <- n

ggplot(aes(x=RT_bin_avg*1000, y=mean, weight = n, 
           color = participant, group = participant), 
       data = timeslice_range) + 
  geom_point(aes(size = n)) +
  geom_smooth(method = "lm", formula = y ~ poly(x,2), se = FALSE) +
  geom_hline(yintercept = 0.5, lty = "dashed") +
  xlab("Average RT (ms)") +
  ylab("Proportion Correct") +
  ylim(0,1)
ggsave(paste(fdir, "QC/Average_RT.png"))




#~~~~~~~~~~~ signal detection theory ~~~~~~~~~~~~~~~~~~
library(psycho)

SubName <- concatenated_df$participant
SubLists <- unique(SubName) #[c(1,49,97,145)]

# signal detection theory B
n_hit_B <- c()
n_fa_B <- c()
n_miss_B <- c()
n_cr_B <- c()

for (i in 1:length(list_of_files)){
  N_hit <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent' & concatenated_df$character=='B',])
  n_hit_B <- c(n_hit_B,N_hit)
  N_fa <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='B',])
  n_fa_B <- c(n_fa_B,N_fa)
  N_miss <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent'& concatenated_df$character=='B',])
  n_miss_B <- c(n_miss_B,N_miss)
  N_cr <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='B',])
  n_cr_B <- c(n_cr_B,N_cr)
}
indices_B <- psycho::dprime(n_hit_B, n_fa_B, n_miss_B, n_cr_B)

# signal detection theory F
n_hit_F <- c()
n_fa_F <- c()
n_miss_F <- c()
n_cr_F <- c()

for (i in 1:length(list_of_files)){
  N_hit <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent' & concatenated_df$character=='F',])
  n_hit_F <- c(n_hit_F,N_hit)
  N_fa <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='F',])
  n_fa_F <- c(n_fa_F,N_fa)
  N_miss <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent'& concatenated_df$character=='F',])
  n_miss_F <- c(n_miss_F,N_miss)
  N_cr <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='F',])
  n_cr_F <- c(n_cr_F,N_cr)
}
indices_F <- psycho::dprime(n_hit_F, n_fa_F, n_miss_F, n_cr_F)

# signal detection theory K
n_hit_K <- c()
n_fa_K <- c()
n_miss_K <- c()
n_cr_K <- c()

for (i in 1:length(list_of_files)){
  N_hit <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent' & concatenated_df$character=='K',])
  n_hit_K <- c(n_hit_K,N_hit)
  N_fa <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='K',])
  n_fa_K <- c(n_fa_K,N_fa)
  N_miss <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent'& concatenated_df$character=='K',])
  n_miss_K <- c(n_miss_K,N_miss)
  N_cr <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='K',])
  n_cr_K <- c(n_cr_K,N_cr)
}
indices_K <- psycho::dprime(n_hit_K, n_fa_K, n_miss_K, n_cr_K)

# signal detection theory M
n_hit_M <- c()
n_fa_M <- c()
n_miss_M <- c()
n_cr_M <- c()

for (i in 1:length(list_of_files)){
  N_hit <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent' & concatenated_df$character=='M',])
  n_hit_M <- c(n_hit_M,N_hit)
  N_fa <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='M',])
  n_fa_M <- c(n_fa_M,N_fa)
  N_miss <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent'& concatenated_df$character=='M',])
  n_miss_M <- c(n_miss_M,N_miss)
  N_cr <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='M',])
  n_cr_M <- c(n_cr_M,N_cr)
}
indices_M <- psycho::dprime(n_hit_M, n_fa_M, n_miss_M, n_cr_M)

# signal detection theory P
n_hit_P <- c()
n_fa_P <- c()
n_miss_P <- c()
n_cr_P <- c()

for (i in 1:length(list_of_files)){
  N_hit <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent' & concatenated_df$character=='P',])
  n_hit_P <- c(n_hit_P,N_hit)
  N_fa <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='P',])
  n_fa_P <- c(n_fa_P,N_fa)
  N_miss <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent'& concatenated_df$character=='P',])
  n_miss_P <- c(n_miss_P,N_miss)
  N_cr <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='P',])
  n_cr_P <- c(n_cr_P,N_cr)
}

indices_P <- psycho::dprime(n_hit_P, n_fa_P, n_miss_P, n_cr_P)

# signal detection theory T
n_hit_T <- c()
n_fa_T <- c()
n_miss_T <- c()
n_cr_T <- c()

for (i in 1:length(list_of_files)){
  N_hit <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent' & concatenated_df$character=='T',])
  n_hit_T <- c(n_hit_T,N_hit)
  N_fa <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='T',])
  n_fa_T <- c(n_fa_T,N_fa)
  N_miss <- nrow(concatenated_df[concatenated_df$key_resp.corr==0 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='congruent'& concatenated_df$character=='T',])
  n_miss_T <- c(n_miss_T,N_miss)
  N_cr <- nrow(concatenated_df[concatenated_df$key_resp.corr==1 & concatenated_df$participant==SubLists[i] & concatenated_df$type=='incongruent'& concatenated_df$character=='T',])
  n_cr_T <- c(n_cr_T,N_cr)
}

indices_T <- psycho::dprime(n_hit_T, n_fa_T, n_miss_T, n_cr_T)



# visualization of d and C
df_signaldetection <- data.frame(participant = rep(SubLists,times=6), 
                                 character = c(rep("B",times=length(list_of_files)),
                                               rep("F",times=length(list_of_files)),
                                               rep("K",times=length(list_of_files)),
                                               rep("M",times=length(list_of_files)),
                                               rep("P",times=length(list_of_files)),
                                               rep("T",times=length(list_of_files))))
df_signaldetection$d <- c(indices_B[["dprime"]],indices_F[["dprime"]],indices_K[["dprime"]],indices_M[["dprime"]],indices_P[["dprime"]],indices_T[["dprime"]])
df_signaldetection$c <- c(indices_B[["c"]],indices_F[["c"]],indices_K[["c"]],indices_M[["c"]],indices_P[["c"]],indices_T[["c"]])

# plot d'
ggplot(df_signaldetection, aes(x = character, y = d, color = participant, group = participant)) + 
  # geom_point(position = position_jitter(width = 0.2, height = 0), size=2) + 
  geom_point(size=2) +
  geom_line() +
  facet_wrap(~ participant)+
  xlab('character') +
  ylab("d'") +
  labs(title="Sensitivity (DPrime) NSubjects / Letters")
ggsave(paste(fdir, "QC/DPrime_per_subject_per_letter.png"))

# bar plot d'
ggplot(df_signaldetection, aes(x = character, y = d, color = participant, group = participant, fill=participant)) +   
  geom_bar(stat = "identity", position = "dodge",width = 0.7) +
  theme_classic()  +
  labs(title="Sensitivity (DPrime) NSubjects / Letters")
ggsave(paste(fdir, "QC/DPrime_per_subject_per_letter_BAR_PLOT.png"))

# plot c
ggplot(df_signaldetection, aes(x = character, y = c, color = participant, group = participant)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size=2) + 
  # geom_point(size=2) +
  geom_line() +
  # facet_wrap(~ participant)+
  xlab('character') +
  ylab("C") +
  labs(title="Sensitivity Bias - NSubjects / Letters")
ggsave(paste(fdir, "QC/DPrime_Bias_per_subject_per_letter.png"))

# bar plot c
ggplot(df_signaldetection, aes(x = character, y = c, color = participant, group = participant, fill=participant)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.7) +
  theme_classic()  +
  labs(title="Sensitivity Bias - NSubjects / Letters")
ggsave(paste(fdir, "QC/DPrime_Bias_per_subject_per_letter_BAR_PLOT.png"))



#~~~~~~~~~ (4.1) plot group average of d for each letter
df_average_d <- data.frame(letter =c('B','F','K','M','P','T'),
                           Mean = c(mean(df_signaldetection$d[1:4]),mean(df_signaldetection$d[5:8]),mean(df_signaldetection$d[9:12]),mean(df_signaldetection$d[13:16]),mean(df_signaldetection$d[17:20]),mean(df_signaldetection$d[21:24])),
                           SD = c(sd(df_signaldetection$d[1:4]),sd(df_signaldetection$d[5:8]),sd(df_signaldetection$d[9:12]),sd(df_signaldetection$d[13:16]),sd(df_signaldetection$d[17:20]),sd(df_signaldetection$d[21:24])))

df_average_d$SE <- df_average_d$SD/sqrt(nrow(df_signaldetection[df_signaldetection$character=='B',]))


plot_avg_d <- ggplot(df_average_d, aes(x = letter, y=Mean, fill=letter)) + 
  geom_bar(stat="identity",
           width = 0.5, position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0, linetype = "dashed",
                position=position_dodge(.9)) + 
  labs(title="Group-average sensitivity for each letter", x="Letter", y = "d'") +
  theme_classic() +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.92, 0.9),
        plot.tag = element_text(size = 10) ) +
  scale_y_continuous(limits = c(0,3.5), breaks = seq(0,3.5, by = 0.5),expand = c(0,0))
ggsave(paste(fdir, "QC/Group_average_sensitivity_for_each_letter.png"))
print(plot_avg_d)



#~~~~~~~~~ (4.2) plot group average of c for each letter
df_average_c <- data.frame(letter =c('B','F','K','M','P','T'),
                           Mean = c(mean(df_signaldetection$c[1:4]),mean(df_signaldetection$c[5:8]),mean(df_signaldetection$c[9:12]),mean(df_signaldetection$c[13:16]),mean(df_signaldetection$c[17:20]),mean(df_signaldetection$c[21:24])),
                           SD = c(sd(df_signaldetection$c[1:4]),sd(df_signaldetection$c[5:8]),sd(df_signaldetection$c[9:12]),sd(df_signaldetection$c[13:16]),sd(df_signaldetection$c[17:20]),sd(df_signaldetection$c[21:24])))

df_average_c$SE <- df_average_c$SD/sqrt(nrow(df_signaldetection[df_signaldetection$character=='B',]))


plot_avg_c <- ggplot(df_average_c, aes(x = letter, y=Mean, fill=letter)) + 
  geom_bar(stat="identity",
           width = 0.5, position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0, linetype = "dashed",
                position=position_dodge(.9)) + 
  labs(title="Group-average bias for each letter", x="Letter", y = "C") +
  theme_classic() +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.92, 0.9),
        plot.tag = element_text(size = 10) ) +
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5, by = 0.1),expand = c(0,0))+
  annotate(geom="text",x = 1.4, y = -0.47, label = "Liberal (more 'yes' responses)", size = 3, angle = 0) +
  annotate(geom="text",x = 1.5, y = 0.47, label = "Conservative (more 'no' responses)", size = 3, angle = 0) +
  geom_hline(yintercept = 0.0, lty = "solid") 
ggsave(paste(fdir, "QC/Group_average_bias_for_each_letter.png"))
print(plot_avg_c)








#~~~~~~~~~~ Drift diffusion model ~~~~~~~~~~~~~~~~~~

library(hBayesDM)

# create data frame for DDM
DDMdf <- concatenated_df
colnames(DDMdf)[colnames(DDMdf) == "key_resp.corr"] <- "choice" # change the name of columns
colnames(DDMdf)[colnames(DDMdf) == "participant"] <- "subjID"
colnames(DDMdf)[colnames(DDMdf) == "key_resp.rt"] <- "RT"

# remove outlier RT < 300ms
outlier_idx <- which(DDMdf$RT > 0.3)
outlier_idx <- outlier_idx[1:length(DDMdf$RT)]
DDMdf$RT <- DDMdf$RT[outlier_idx]

# specify choice number
DDMdf$choice <- gsub("1", "2", DDMdf$choice) # 2 = correct
DDMdf$choice <- gsub("0", "1", DDMdf$choice) # 1 = incorrect



# ~~~~~~~~~ Run the model with a given data.frame (DDMdf)
iterations <- 1500
output <- choiceRT_ddm(data = DDMdf, 
                       niter = iterations, nwarmup = iterations*0.2, 
                       nchain = 4, ncore = 4)




# ~~~~~~~~~~ visualizaton of DDM
# (5.1) Visually check convergence of the sampling chains (should look like 'hairy caterpillars')
sanity_4000 <- plot(output, type = "trace")

# check the mean of each subject's parameters
output$allIndPars

# plot each parameter with multiple bar plots in one figure
library(gridExtra)

# plot the each parameter
plot1 <- ggplot(output$allIndPars, aes(x=subjID, y=alpha, fill = subjID)) + 
  geom_bar(stat="identity",
           width = 0.7, position=position_dodge()) +
  labs(title="Alpha per participant", x="Participants", y = "Alpha") +
  theme_classic()
plot1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot1 <- plot1 + guides(fill = "none")
ggsave(paste(fdir, "QC/hBayesDDModeling_Alpha_per_participant.png"))

plot2 <- ggplot(output$allIndPars, aes(x=subjID, y=beta, fill=subjID)) + 
  geom_bar(stat="identity",
           width = 0.7, position=position_dodge()) +
  labs(title="Beta per participant", x="Participants", y = "Beta") +
  theme_classic() 
plot2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot2 <- plot2 + guides(fill = "none")
ggsave(paste(fdir, "QC/hBayesDDModeling_Beta_per_participant.png"))

plot3 <- ggplot(output$allIndPars, aes(x=subjID, y=delta, fill=subjID)) + 
  geom_bar(stat="identity",
           width = 0.7, position=position_dodge()) +
  labs(title="Delta per participant", x="Participants", y = "Delta") +
  theme_classic() 
plot3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot3 <- plot3 + guides(fill = "none")
ggsave(paste(fdir, "QC/hBayesDDModeling_Delta_per_participant.png"))

plot4 <- ggplot(output$allIndPars, aes(x=subjID, y=tau, fill=subjID)) + 
  geom_bar(stat="identity",
           width = 0.7, position=position_dodge()) +
  labs(title="Tau per participant", x="Participants", y = "Tau") +
  theme_classic() 
plot4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot4 <- plot4 + guides(fill = "none")
ggsave(paste(fdir, "QC/hBayesDDModeling_Tau_per_participant.png"))

# concatanate all the plots
png(file=paste(fdir, "QC/hBayesDDModeling_per_participant.png"),
    width=1800, height=1800)
combined_plot <- grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
# Display the concatenated plot

dev.off()
print(combined_plot)


# ~~~~~~~ plot each parameter of each individual
df_eachparameter <- output$parVals
#df_SUBname <- data.frame(participant = c(rep('anonymous',times=12000),
#                                         rep('anonymous2',times=12000),
#                                         rep('EDX3',times=12000),
#                                         rep('EVY3',times=12000)))

# alpha data frame
df_alpha <- data.frame(as.vector(df_eachparameter$alpha),SubLists)
colnames(df_alpha) <- c("alpha", "participant")

# beta data frame
df_beta <- data.frame(as.vector(df_eachparameter$beta),SubLists)
colnames(df_beta) <- c("beta", "participant")

# delta data frame
df_delta <- data.frame(as.vector(df_eachparameter$delta),SubLists)
colnames(df_delta) <- c("delta", "participant")

# tau data frame
df_tau <- data.frame(as.vector(df_eachparameter$tau),SubLists)
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
    point_colour = NA) +
  geom_boxplot(
    width = 0.2,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5) +
  labs(title="4000\nResponse caution (alpha)", x="Participant", y = "Alpha [a.U.]") +
  scale_y_continuous(limits = c(0,4.5), breaks = seq(0, 6, by = 1), expand = c(0,0))
alpha_4000 <- alpha_4000 + guides(fill = "none")
alpha_4000 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste(fdir, "QC/hBayesDDModeling_Response caution_alpha.png"))


# ~~~~~~~~ (5.2) plot group-average of alpha
df_averageAlpha <- data.frame(x_axis = c('1'),
                              Mean = mean(df_alpha$alpha), 
                              SD = sd(df_alpha$alpha))
df_averageAlpha$SE <- df_averageAlpha$SD/sqrt(nrow(df_alpha))

plot_avgAlpha <- ggplot(df_averageAlpha, aes(x = x_axis, y=Mean)) + 
  geom_bar(stat="identity",
           width = 0.5, position=position_dodge(), fill="transparent", color="black") +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0, linetype = "dashed",
                position=position_dodge(.9)) + 
  labs(title="Group-average response caution (alpha)", x="Sessions", y = "Alpha [a.U.]") +
  theme_classic() +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.85, 0.9),
        plot.tag = element_text(size = 10) ) +
  scale_y_continuous(limits = c(0,3), breaks = seq(0,3, by = 0.5),expand = c(0,0))
ggsave(paste(fdir, "QC/hBayesDDModeling_Group_average_response_caution_alpha.png"))
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
  geom_boxplot(
    width = 0.2,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  labs(
    title=" \nResponse bias (beta)", x="Participant", y = "Beta [a.U.]"
  ) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2), expand = c(0,0)
  ) +
  annotate(geom="text",x = 1, y = 0.03, label = "incorrect", size = 5, angle = 0) +
  annotate(geom="text",x = 1, y = 0.97, label = "correct", size = 5, angle = 0)
beta_4000 <- beta_4000 + guides(fill = "none")
beta_4000 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste(fdir, "QC/hBayesDDModeling_Response bias_beta.png"))

# ~~~~~~~~~~ (5.2) plot average of beta
df_averageBeta <- data.frame(x_axis = c('1'),
                             Mean = mean(df_beta$beta),
                             SD = sd(df_beta$beta))
df_averageBeta$SE <- df_averageBeta$SD/sqrt(nrow(df_beta))

plot_avgBeta <- ggplot(df_averageBeta, aes(x = x_axis, y=Mean)) + 
  geom_bar(stat="identity",
           width = 0.5, position=position_dodge(), fill="transparent", color="black") +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0, linetype = "dashed",
                position=position_dodge(.9)) + 
  labs(title="Group-average response bias (beta)", x="Sessions", y = "Beta [a.U.]") +
  theme_classic() +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.85, 0.9),
        
        plot.tag = element_text(size = 10) ) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.2),expand = c(0,0))
ggsave(paste(fdir, "QC/hBayesDDModeling_Group_average_Response bias_beta.png"))
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
    point_colour = NA) +
  geom_boxplot(
    width = 0.2,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5) +
  labs(title=" \nSpeed of evidence accumulation (delta)", x="Participant", y = "Delta [a.U.]") +
  scale_y_continuous(limits = c(0,3), breaks = seq(0, 6, by = 1), expand = c(0,0))
delta_4000 <- delta_4000 + guides(fill = "none")
delta_4000 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste(fdir, "QC/hBayesDDModeling_Speed_of_evidence_accumulation_delta.png"))

# ~~~~~~~~~~~~ (5.2) plot average delta
df_averageDelta <- data.frame(x_axis = c('1'),
                              Mean = mean(df_delta$delta),
                              SD = sd(df_delta$delta))

df_averageDelta$SE <- df_averageDelta$SD/sqrt(nrow(df_delta))

plot_avgDelta <- ggplot(df_averageDelta, aes(x = x_axis, y=Mean)) +
  geom_bar(stat="identity",
           width = 0.5, position=position_dodge(), fill="transparent", color="black") +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0, linetype = "dashed",
                position=position_dodge(.9)) + 
  labs(title="Group-average speed of evidence accumulation (delta)", x="Sessions", y = "Delta [a.U.]") +
  theme_classic() +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.85, 0.9),
        plot.tag = element_text(size = 10) ) +
  scale_y_continuous(limits = c(0,4), breaks = seq(0,4, by = 0.5),expand = c(0,0))
ggsave(paste(fdir, "QC/hBayesDDModeling_Group_average_Speed_of_evidence_accumulation_delta.png"))
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
    point_colour = NA) +
  geom_boxplot(width = 0.2, outlier.color = NA, alpha = 0.5 ) +
  labs(title=" \nSensory encoding and motor execution", x="Participant", y = "Tau [s]" ) +
  scale_y_continuous(limits = c(0,1.75), breaks = seq(0, 1.75, by = 0.1), expand = c(0,0))
tau_4000 <- tau_4000 + guides(fill = "none")
tau_4000 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste(fdir, "QC/hBayesDDModeling_Sensory_encoding_and_motor_execution_Tau.png"))


# ~~~~~~~~~~~ (5.2) plot average tau
df_averageTau <- data.frame(x_axis = c('1'),
                            Mean = mean(df_tau$tau),
                            SD = sd(df_tau$tau))
df_averageTau$SE <- df_averageTau$SD/sqrt(nrow(df_tau))
a<-df_averageTau$Mean-df_averageTau$SE
b<-df_averageTau$Mean+df_averageTau$SE
c<-df_averageTau$Mean-df_averageTau$SD
d<-df_averageTau$Mean+df_averageTau$SD

plot_avgTau <- ggplot(df_averageTau, aes(x = x_axis, y=df_averageTau$Mean)) + 
  geom_bar(stat="identity", width = 0.5, position=position_dodge(), fill="transparent", color="black") +
  geom_errorbar(aes(ymin=a, ymax=b),
                width=0, position=position_dodge(.9)) + 
  geom_errorbar(aes(ymin=c, ymax=d),
                width=0, linetype = "dashed", position=position_dodge(.9)) + 
  labs(title="Group-average sensory encoding and motor execution (tau)", x="Sessions", y = "Tau [s]") +
  theme_classic() + labs(tag = "──  SE\n- - -  SD") + 
  theme(plot.tag.position = c(0.85, 0.9), plot.tag = element_text(size = 10) ) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.2),expand = c(0,0))
ggsave(paste(fdir, "QC/hBayesDDModeling_Group_average_sensory_encoding_and_motor_execution_Tau.png"))
print(plot_avgTau)



# ~~~~~~~~ Plot together
# concatanate all the plots
png(file=paste(fdir, "QC/hBayesDDModeling_Parameters_per_Subject.png"),
    width=1800, height=1800)


# concatanate all the parameter for each participant
combined_plot <- grid.arrange(alpha_4000, beta_4000, delta_4000, tau_4000, nrow = 1, ncol = 4)
dev.off()
print(combined_plot)

# concatanate all the average plots
png(file=paste(fdir, "QC/hBayesDDModeling_Group_average.png"),
    width=1800, height=1200)
combined_plot_Avg <- grid.arrange(plot_avgAlpha, plot_avgBeta, plot_avgDelta, plot_avgTau, nrow = 1, ncol = 4)
dev.off()
print(combined_plot_Avg)




## ~~~~~~~~~~~~~~~~~ accuracy (congruent/incongruent) ~~~~~~~~~~~~~~~~~

SubName <- concatenated_df$participant
SubName<- unique(SubName)
SUBplot_garage <- list()
for (i in 1:length(SubLists)){
  print(i)
  SUB <- concatenated_df[concatenated_df$participant==SubLists[i],]
  SUB_mean <- data_frame(type=c('Congruent','Incongruent', 'Total'),
                         mean = c(mean(SUB$key_resp.corr[SUB$type=='congruent']),
                                  mean(SUB$key_resp.corr[SUB$type=='incongruent']),
                                  mean(SUB$key_resp.corr)))
  variable_name <- paste("SUB", i, sep = "_") 
  
  SUB_plot <- ggplot(SUB_mean, aes(x=type, y=mean*100, fill = type)) + 
    geom_bar(stat="identity", width = 0.7, position=position_dodge()) + 
    labs(title=SubLists[i], x="Stimulus", y = "Response Accuracy (%)") +
    theme_classic() + geom_hline(yintercept = chance_level, linetype = "dotted", color = "black") +
    annotate(geom="text",x = 2.5, y = 48, label = "Chance level", size = 4, angle = 90) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100, by = 25),expand = c(0,0))
  
  SUB_plot <- SUB_plot + guides(fill = "none") 
  SUBplot_garage[[variable_name]] <- SUB_plot
  
  ggsave(paste(fdir, paste("QC/perSubject/Accuracy_", SubLists[i],".png")))
}



# ~~~~~~~~~~~~~~ (1.3) plot average-accuracy for each stimulus
Stimu_Avg <- concatenated_df
Stimu_Avg_df <- data_frame(type=c('Congruent','Incongruent','Total'),
                           Mean = c(mean(Stimu_Avg$key_resp.corr[Stimu_Avg$type=='congruent']),
                                    mean(Stimu_Avg$key_resp.corr[Stimu_Avg$type=='incongruent']),
                                    mean(Stimu_Avg$key_resp.corr)),
                           SD = c(sd(Stimu_Avg$key_resp.corr[Stimu_Avg$type=='congruent']),
                                  sd(Stimu_Avg$key_resp.corr[Stimu_Avg$type=='incongruent']),
                                  sd(Stimu_Avg$key_resp.corr)))

Stimu_Avg_df$SE <- Stimu_Avg_df$SD/c(sqrt(nrow(Stimu_Avg[Stimu_Avg$type=='congruent',])))

plot_avgStimu <- ggplot(Stimu_Avg_df, aes(x =type, y=Mean*100, fill = type)) +
  geom_bar(stat="identity", width = 0.5, position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean*100-SE*100, ymax=Mean*100+SE*100), width=0,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=Mean*100-SD*100, ymax=100), width=0, linetype = "dashed",
                position=position_dodge(.9)) + 
  labs(title="Group-average accuracy", x="Stimulus", y = "Response Accuracy (%)") +
  theme_classic() +
  geom_hline(yintercept = chance_level, linetype = "dotted", color = "black") +
  annotate(geom="text",x = 2.5, y = 48, label = "Chance level", size = 4, angle = 90) +
  labs(tag = "──  SE\n- - -  SD") +
  theme(plot.tag.position = c(0.85, 0.9), plot.tag = element_text(size = 10) ) +
  scale_y_continuous(limits = c(0,100.), breaks = seq(0,100, by = 10),expand = c(0,0))
ggsave(paste(fdir, paste("QC/Group_average_accuracy_per_Stim.png")))
print(plot_avgStimu)
