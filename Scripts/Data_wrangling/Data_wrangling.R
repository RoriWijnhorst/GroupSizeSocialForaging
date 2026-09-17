library(dplyr)
library(tidyr)

# load group trial data and merge with trial information
df_group <- read.csv("Data/Raw_data/PS_RFID_group_calibration_data_2023.csv")
str(df_group)
Trial_ID_merged <- read.csv("Data/Raw_data/Trial_ID_merged.csv")
Trial_ID_merged<-Trial_ID_merged %>% 
  data.frame %>%
  unite(TrialID,c("Trial.number", "Feeder"), sep = "_", remove=F)

# Add scores from both feeder boards together
df_group <- df_group[,c("Tag.ID", "TrialID", "Dir_scr_visits_n", "Prim_prd_visits_n",  "Total_feeds_RFID")]
colnames(df_group) <-c("Tag.ID","TrialID","Scr", "Prd", "Feeds")

# Add scores from both feeder boards together
df_group<- df_group %>% data.frame %>% separate(TrialID, into = c("Trial.number", "Feeder"), sep = "_")
df_group<- df_group %>% data.frame %>% separate(Trial.number, into = c("G", "TrialID"), sep = 1) 

df_group <- df_group %>% group_by(G, Tag.ID, TrialID) %>% summarise(across(c(Scr, Prd, Feeds),sum),
                                                                    .groups = 'drop') 


dftrial <- data.frame(unique(Trial_ID_merged$Trial.number))
colnames(dftrial) <- "X"
dftrial<- dftrial %>% 
  data.frame %>%
  separate(X, into = c("G", "TrialID"), sep = 1)
dftrial$TrialID <- as.integer(dftrial$TrialID)
dftrial <- dftrial[order(dftrial$TrialID),]

dftrial$Seq <- rep(rep(1:10, each=2),6)


## Group trial first (0: no, 1: yes)
dftrial$Grp <- c(rep(rep(0:1, each=20),2), rep(rep(0:1, each=10),2))
dftrial$Group <- c(rep(32:31, 10), rep(34:33, 10), rep(41:42, 10), rep(44:43, 10), rep(51:50, 10), rep(52:53, each=10))
df_group<- merge(df_group,dftrial[,c("TrialID","Seq", "Grp", "Group")], all.x=T)


# get triadic trial data 
df_triad <- read.csv("Data/Raw_data/Rstan_PS_rstr_data_22_23_25.04.2024.csv")
df_triad <- df_triad[,c("Tag.ID", "RingNR","Year","Group_index","TrialID_index","Trial_day_bin","TrialNR_Indiv", "Grp_tr_b4", "TrialID","Focal_ID_index", "Opponent1_ID_index", "Opponent2_ID_index", "Dir_scr_visits_n", "Prim_prd_visits_n", "Total_feeds_RFID")]
colnames(df_triad)<-c("Tag.ID", "RingNR","Year","Group","Trial","Day","Seq","Grp", "TrialID", "IDi", "IDj", "IDk","Scr", "Prd", "Feeds")
df_triad <- filter(df_triad, df_triad$Tag.ID %in% df_group$Tag.ID, Year!=2022) # exclude triadic trial from 2022

# add group RingNR and Tag.ID to group trial data
df_info <- unique(df_triad[,c("Tag.ID", "RingNR","Group","IDi")])
df_group<- merge(df_group, df_info, by=c("Tag.ID", "Group"), all.x=T)

## Calculate proportion scrounging (optional, but done safely)
df_triad$Events  <- df_triad$Scr + df_triad$Prd
df_group$Events <- df_group$Scr + df_group$Prd

df_triad  <- df_triad[df_triad$Events  > 0, ]
df_group <- df_group[df_group$Events > 0, ]

df_triad$Prop  <- df_triad$Scr  / df_triad$Events
df_group$Prop <- df_group$Scr / df_group$Events

## For group assays
# Remove trials with only 1 individual
df_group <- df_group[ave(df_group$Tag.ID, df_group$TrialID, FUN = length) > 1, ]

# Partner totals: Scr_opp and Events_opp
df_group$Scr_opp <- ave(df_group$Scr, df_group$TrialID, FUN = sum) - df_group$Scr
df_group$Prd_opp <- ave(df_group$Prd, df_group$TrialID, FUN = sum) - df_group$Prd
df_group$Events_opp <- ave(df_group$Events, df_group$TrialID, FUN = sum) - df_group$Events

# Partner exact proportion
df_group$Prop_opp <- df_group$Scr_opp / df_group$Events_opp

## For triadic assays
# Remove trials with only 1 individual
df_triad <- df_triad[ave(df_triad$Tag.ID, df_triad$TrialID, FUN = length) > 1, ]

# Partner totals: Scr_opp and Events_opp
df_triad$Scr_opp <- ave(df_triad$Scr,    df_triad$TrialID, FUN = sum) - df_triad$Scr
df_triad$Prd_opp <- ave(df_triad$Prd, df_triad$TrialID, FUN = sum) - df_triad$Prd
df_triad$Events_opp <- ave(df_triad$Events, df_triad$TrialID, FUN = sum) - df_triad$Events

# (Optional) partner mean proportion
df_triad$Prop_opp <- df_triad$Scr_opp / df_triad$Events_opp

## New trial index: Triadic
df_t <- data.frame(unique(df_triad[,"TrialID"]))
colnames(df_t) <- "TrialID"
df_t$Trial_new <- seq_along(df_t$TrialID)
df_triad  <- merge(df_triad,  df_t)

## New trial index: Group
df_t <- data.frame(unique(df_group[,"TrialID"]))
colnames(df_t) <- "TrialID"
df_t$Trial_new <- seq_along(df_t$TrialID)
df_group <- merge(df_group, df_t)

## New ID index
df_t <- data.frame(unique(df_triad[,"IDi"]))
colnames(df_t) <- "IDi"
df_t$IDi_new <- seq_along(df_t$IDi)
df_triad  <- merge(df_triad,  df_t)
df_group <- merge(df_group, df_t)

## New Group index
df_t <- data.frame(unique(df_group[,"Group"]))
colnames(df_t) <- "Group"
df_t$Group_new <- seq_along(df_t$Group)
df_triad  <- merge(df_triad,  df_t)
df_group <- merge(df_group, df_t)


df_triad <- subset(df_triad, !is.na(Feeds))
df_group <- subset(df_group, !is.na(Feeds))

## Data list for Stan: use partner counts instead of Prop_opp
dfl <- list(
  N_z1 = nrow(df_triad),
  N_z2 = nrow(df_group),
  Ng   = length(unique(df_triad$Group)),
  Nt_t = length(unique(df_triad$TrialID)),
  Nt_g = length(unique(df_group$TrialID)),
  
  Day   = df_triad$Day*1-0.5,
  Seq_t = as.vector(scale(df_triad$Seq,scale=F),),
  Grp_t = df_triad$Grp*1-0.5,
  Seq_g = as.vector(scale(df_group$Seq,scale=F),),
  Grp_g = df_group$Grp*1-0.5,
  
  Trial_t = df_triad$Trial_new,
  Group_t = df_triad$Group_new,
  Trial_g = df_group$Trial_new,
  Group_g = df_group$Group_new,
  
  ind_z1 = df_triad$IDi_new,
  ind_z2 = df_group$IDi_new,
  J      = length(unique(df_triad$RingNR)),
  
  # focal responses
  prop1 = df_triad$Prop,
  prop2 = df_group$Prop, 
  Scr1 = df_triad$Scr,
  Scr2 = df_group$Scr,
  Prd1 = df_triad$Prd,
  Prd2 = df_group$Prd,
  n_events1 = df_triad$Events,
  n_events2 = df_group$Events,
  Events1 = as.vector(scale(df_triad$Events,scale=F),),
  Events2 = as.vector(scale(df_group$Events,scale=F),),
  
  # opponents
  Scr_opp1   = df_triad$Scr_opp,
  Scr_opp2   = df_group$Scr_opp,
  Prd_opp1   = df_triad$Prd_opp,
  Prd_opp2   = df_group$Prd_opp,
  For_opp1   = df_triad$Events_opp,   # total partner events
  For_opp2   = df_group$Events_opp,    # total partner events
  Events_opp1 = as.vector(scale(df_triad$Events_opp,scale=F),),
  Events_opp2 = as.vector(scale(df_group$Events_opp,scale=F),),
  
  # Feed consumed
  Rel_fit1 = (df_triad$Feeds)/mean(c(df_triad$Feeds,df_group$Feeds)),
  Rel_fit2 = (df_group$Feeds)/mean(c(df_triad$Feeds,df_group$Feeds))
)

saveRDS(dfl, "Data/Wrangled_data/STAN_input_data.rds")
