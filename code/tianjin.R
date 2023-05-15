library(survminer)
library(survival)
library(tidyverse)
library(lubridate)
library(icenReg)
library(igraph)
library(visNetwork)
library(mvtnorm)
library(viridis)
options(digits=3)
set.seed(1235) 

setwd("/Users/padmashri/Desktop/IDD/Project_COVID")

# Load the data
data <- read_csv("Tianjin135cases_revised.csv")

# Convert the date variable to a Date object
data$date_confirmation <- as.Date(data$confirm_date, format = "%d/%m/%Y")
data$date_onset_symptoms <- as.Date(data$symptom_onset, format = "%d/%m/%Y")

# Calculate daily counts
daily_counts <- data %>%
  count(date_confirmation) %>%
  rename(date = date_confirmation, count = n)

tj_plotinc <- daily_counts

N <- 15000000 # population at risk
C <- sum(daily_counts$count)
total_i <- 0
for(i in (N-C+1):(N-1)){
  total_i <- 1/i + total_i
}

R01.tj <- ((N-1)/C)*total_i

#daily_counts2 <- daily_counts[1:26,] 

# Calculate the exponential growth rate
growth_rate <- lm(log(count) ~ decimal_date(date), data = daily_counts)$coefficients[2]

tdata=read.csv("Tianjin135cases_revised.csv",na.strings = "", stringsAsFactors = F)
tdata$symptom_onset=as.Date(tdata$symptom_onset, format = "%d/%m/%Y")
tdata$start_source=as.Date(tdata$start_source, format = "%d/%m/%Y")
tdata$end_source=as.Date(tdata$end_source,format = "%d/%m/%Y" )
tdata$confirm_date=as.Date(tdata$confirm_date,format = "%d/%m/%Y" )

# remove cases that do NOT have a date of symptom onset to determine the ICC. 
# There are 10 cases without symptom onset, but all have a confirmation date. 
# We will re-run the estimates with imputed data later

# Make a copy of the original dataset for use later on in imputed serial interval estimates
tdata_org <- tdata

# Remove any cases that are missing data for date of symptom onset
tdata <- tdata[which(!is.na(tdata$symptom_onset)),]  #removes 10 observations

# fix a data entry error in Infection_source column, where TJ is listed as JN. 
tdata$Infection_source <- str_replace(tdata$Infection_source, pattern = "JN", replacement = "TJ")

# Turn 'Infection_source' into lower case and get trim any whitespace so don't have issues with case sensitivity, etc
tdata$Infection_source <- str_to_lower(tdata$Infection_source)
tdata$Infection_source <- str_trim(tdata$Infection_source)
#table(tdata$Infection_source)
sum(is.na(tdata$Infection_source)) #1 NA

#Create a duplicated Infection_source column for the purposes of string manipulation so we can appropriately group sources of infection
tdata$Infection_source_dup <- tdata$Infection_source

#Change the string "person" to "individual" or else it will get picked up as a relative (matches "son") in the code below
tdata$Infection_source_dup <- str_replace_all(tdata$Infection_source_dup, pattern = "person", replacement = "individual")

#Case TJ27 is infected by a coworker from Wuhan, so we want the source label to be coworker, not Wuhan
#so need to remove "Wuhan" from the reason or it will get labeled as Wuhan below
tdata$Infection_source_dup <- str_replace(tdata$Infection_source_dup, 
                                          pattern = "coworker of a individual from wuhan",
                                          replacement = "coworker")

#Note that the order the data are selected in is VERY important to which case goes into which source_group category
#For those that meet multiple criteria (e.g. wuhan; tj1), the str_match which is highest in the case_when call (i.e. "wuhan|hubei") will have priority over those matching later 
#so that the 'source' column contain "wuhan; tj1" would be labelled as infection from a "wuhan" rather than from a "known relationship" origin 
#There are only a small number of cases for which this matters (n = 12)

#We will emphasize the wuhan and travel cases over known relationships
#This seems logical, given that the epicenter of the outbreak was Wuhan
tdata <- mutate(tdata, source_group = case_when(!is.na(str_match(Infection_source_dup, "wuhan|hubei")) ~ "Wuhan and Hubei", #Priority 1
                                                !is.na(str_match(Infection_source_dup, "mall|store|shopper|shopping")) ~ "Mall", #Priority 1
                                                !is.na(str_match(Infection_source_dup, "family|relative|wife|mother|son|sister|daughter|brother|husband|duaghtor|whife|hunsband")) ~ "Relative", #Priority 2
                                                !is.na(str_match(Infection_source_dup, "coworker|business|workplace|colleague|colleage")) ~ "Coworker", #Priority 2
                                                !is.na(str_match(Infection_source_dup, "tj|patient")) ~ "Other relationship", #Priority 2
                                                !is.na(str_match(Infection_source_dup, "train|travel|trip|hebei|dalian")) ~ "Other travel", #Priority 3
                                                !is.na(str_match(Infection_source_dup, "unknown|unclear")) ~ "Unknown", #Priority 5
                                                is.na(Infection_source_dup) ~ "Unknown", #Priority 5
                                                T ~ "other")) #there should be none of these, so this is just a sanity check!  

#What is distribution of probably source of infection (grouped)?
table(tdata$source_group) 

mynodes <- tdata$case_id

#Transform everything to lower case to make sure there aren't any issues with matching due to case inconsistencies
mynodes <- str_to_lower(mynodes) 
tdata$case_id <- str_to_lower(tdata$case_id)

edges = data.frame(from=mynodes[9],to=mynodes[21],stringsAsFactors = F) # i read this one manually 

for (id in 1:nrow(tdata)) {
  tonode=tdata$case_id[id]
  fromnodes=str_extract_all(tdata$Infection_source[id], "tj\\d+", simplify = T) #in lower case due to above early/late split on infection source
  if (length(fromnodes)>0) {
    for (k in 1:length(fromnodes)) {
      edges=rbind(edges, c(fromnodes[k], tonode))
    }
  }
}
head(edges)
edges=edges[-1,] #Remove the initial relationship we gave so it isn't duplicated
edges=edges[-which(is.na(edges[,1])),] # NAs arose from a few empty entries for Infection_source 

# Make a smaller dataset of original tpdata that contains only the CaseID and date of symptom onset
tdata_sympt <- tdata[,c(1,4)]

# Add the date of symptom onset -for the caseID of the 'from' case - to the case pairs dataset (edges)
#Do some renaming so the join is based on the caseID in the from column and that name of date column reflects this
names(tdata_sympt) <- str_replace(names(tdata_sympt), "case_id", "from")
undir_tdates <- left_join(edges, tdata_sympt, by = "from")
names(undir_tdates) <- str_replace(names(undir_tdates), "symptom_onset", "from_sympt_date")

# Repeat, but add the date of symptom onset for the caseID of the 'to' case
names(tdata_sympt) <- str_replace(names(tdata_sympt), "from", "to")
undir_tdates <- left_join(undir_tdates, tdata_sympt, by = "to")
names(undir_tdates) <- str_replace(names(undir_tdates), "symptom_onset", "to_sympt_date")

# Now add some extra columns which give us the raw serial interval (i.e. number of days between symptom onset in infector-infectee pairs)
#As well as the absolute value of the serial interval (as some cases in the "from" and "to" columns should be switched around!)
#And finally a 'direction' column in case we need to sort out which directions the arrows should be going in for a network graph and where we have missing dates
undir_tdates <- mutate(undir_tdates, earliest_sympt_onset = pmin(to_sympt_date, from_sympt_date, na.rm = T), 
                       raw_serial_interval = to_sympt_date - from_sympt_date,   
                       abs_serial_interval = abs(raw_serial_interval))

# Split dataset into positive (or 0) serial interval vs. negative vs. NA 
#A negative serial interval means our "to" and "from" cases are mixed up
pos <- filter(undir_tdates, raw_serial_interval >= 0)
neg <- filter(undir_tdates, raw_serial_interval < 0)
onlyone <- filter(undir_tdates, is.na(raw_serial_interval)) #3 NAs where date of symptom onset is not known; keep for making columns for now

names(neg)
names(neg)[1] <- "to"
names(neg)[2] <- "from"
names(neg)[3] <- "to_sympt_date"
names(neg)[4] <- "from_sympt_date"

# Now bind the rows of the seperated datasets back together based on column names
#Must use dplyr::bind_rows to bind based on column name rather than position
undir_tdates <- bind_rows(pos, neg, onlyone)

# For plotting - Add a column with padded to and from caseID numbers so they print in numerical order
#Add a zero on the left of the number so all numbers have three digits  
undir_tdates$pto <- str_replace(undir_tdates$to, pattern = "tj", replacement = "")
undir_tdates$pto <- str_pad(undir_tdates$pto, width = 3, side = "left", pad = "0")

undir_tdates$pfrom <- str_replace(undir_tdates$from, pattern = "tj", replacement = "")
undir_tdates$pfrom <- str_pad(undir_tdates$pfrom, width = 3, side = "left", pad = "0")

# For plotting - Make a new column with case pair ID
undir_tdates <- mutate(undir_tdates, pairID = factor(paste("tj", pfrom, "-", "tj", pto, sep = "")))

rm(pos, neg, onlyone)

# Make data frame of edges, where the cases as the 'earliest' date of symptom onset are labeled as the "from" cases
#Filter down to cases that do NOT have NAs for date of symptom onset (some added back in by being related cases)
tedges <- filter(undir_tdates, !is.na(raw_serial_interval))
tedges <- tedges[,c(1,2)]
tedges$arrows <- "to" 

# Select down to only the unique cases-pairs 
#Note that if you keep the negative serial intervals as negatives, there are no duplicates; 
#But instead there are two instances where you get a loop (tj114-tj115 and tj96-tj97)
#ex) tj114 is supposed to infect tj115 and tj115 is supposed to infect tj114
tedges <- distinct(tedges)

nodes = data.frame(id=tdata$case_id, 
                   label=tdata$case_id,
                   group=tdata$source_group)

#### contact graph ####
visNetwork(nodes, tedges) %>% visLegend()

##### ICC #####
#The interval case to case (ICC) data are the times between the (presumed) 
#index case for a small cluster and the other cases in the cluster. 
#The Vink et al approach allows these intervals to be one of 4 types, and 
#estimates the serial interval and the probability of each type. 
#To extract ICC intervals, we let the clusters be the components of the graph, 
#and we let the presumed index case be the first to develop symptoms. 
#For each cluster, we subtract the index cases' symptom time from the symtom times of 
#the rest of the cluster (or just the first few; it turns out that the estimate is not 
#sensitive to this). This results in a list of time intervals between symptom onset in 
#presumed index cases and symptom onset in other cases in the same cluster (graph component). 

# construct graph
tgraph = graph_from_edgelist(as.matrix(edges[,1:2]), directed = FALSE)
ccs=components(tgraph)

tdata$component=vapply(tdata$case_id, function(x)
{ if (x %in% names(ccs$membership)) { return(ccs$membership[match(x, names(ccs$membership))])
} else { 
  return(NA)}}, FUN.VALUE = 3)

# extract ICC interval
getICCs <- function(thisdata, ccs, K, orderby= "onset" ) {
  iccs=1
  for (n in 1:max(ccs$membership)) {
    mycases  = which(thisdata$component==n)
    if(orderby == "onset")
    {  myonsets = sort(thisdata$symptom_onset[mycases])[1:min(K, length(mycases))]}
    if(orderby == "exposure") {
      myonsets =thisdata$symptom_onset[mycases][order(thisdata$end_source[mycases])][1:min(K,length(mycases))]
    }
    iccs =c(iccs, myonsets[-1]-myonsets[1])
  }
  return(iccs[-1]) 
}

# remove cases that don't have a symptom onset date
icc3 = getICCs(tdata,ccs,3)
icc4 = getICCs(tdata,ccs,4)
icc5 = getICCs(tdata,ccs,5)
icc6 = getICCs(tdata,ccs,6)
icc_expose = getICCs(tdata, ccs, 4, orderby ="exposure")

#### serial interval estimates ####

# Vinkwallinga 2014

source("TianjinSI_VinkWallinga_CC.R")
myest3 = serial_mix_est(data=icc3, N=100, startmu=10, startsig =4)
myest4 = serial_mix_est(data=icc4, N=100, startmu=10, startsig =4)
myest5 = serial_mix_est(data=icc5, N=100, startmu=10, startsig =4)
myest6 = serial_mix_est(data=icc6, N=100, startmu=10, startsig =4)
myest_exp= serial_mix_est(data=icc_expose, N=100, startmu=10, startsig =4)

mm=rbind(myest3, myest4, myest5,myest6, myest_exp)
colnames(mm)=c("mu","sig")
mm=as.data.frame(mm)
mm$NumCasesPerCluster=c( 3, 4, 5, 6, 4) 
mm$ordering = c("Onset","Onset","Onset","Onset","LastExposure")
print(mm[,c(4,3,1,2)]) 

#The mean SI (using first 4 cases per cluster) is `r myest4[1]`.
# The standard deviation of the serial intervals is `r myest4[2]`.

###  Cleveland dotplot of raw serial intervals per possible case pair 

# We want to exclude any case-pairs that have NA for length of serial interval
#i.e. one of the pair does not have a date of symptom onset
undir_tdates_org <- undir_tdates  #Just in case....
undir_tdates <- filter(undir_tdates, !is.na(raw_serial_interval))

#Pivot the to/from dates of symptom onset column to a long format, so that can make a legend based on this variable
undir_dotplot <- pivot_longer(undir_tdates, 
                              cols = contains("sympt_date"),
                              names_to = "pair_member",
                              values_to = "onset_date")

#Let's rename the values so it makes more sense in the legend
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "from_sympt_date", replacement = "Presumed infector")
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "to_sympt_date", replacement = "Presumed infectee")

#Make the Cleaveland dotplot
p <- ggplot(undir_dotplot, aes(y = reorder(pairID, earliest_sympt_onset))) +
  geom_segment(aes(x = earliest_sympt_onset, xend = earliest_sympt_onset + abs_serial_interval, yend = pairID), 
               color = "#404788FF") +
  geom_point(aes(x = onset_date, color = pair_member, fill = pair_member, shape = pair_member)) +
  scale_x_date(date_breaks = "1 day") +
  scale_color_manual(name = "Pair member for \ndate of symptom onset", values = c("#D44842FF", "#FAC127FF")) +
  scale_fill_manual(name = "Pair member for \ndate of symptom onset", values = c("#D44842FF", "#FAC127FF")) +
  scale_shape_manual(name = "Pair member for \ndate of symptom onset", values = c(23, 21)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey80", linetype = "dashed"),
        panel.background = element_rect(fill = "white")) +
  labs(title = "Serial intervals of possible case pairs from Tianjin, over time",
       x = "Date of symptom onset",
       y = "Case pairs")
p

#### Direct serial interval estimates: all possible case-pairs ####
# To determine an estimate of the raw serial intervals, we will find the mean, median, and standard deviation of all possible case-pairs. Note that this is the raw data where cases without date of symptom onset have been removed (i.e. NOT imputed!).

### Arrange the dataset by the date of earliest symptom onset for each case pair
undir_tdates <- arrange(undir_tdates, earliest_sympt_onset)


fit.ti <- generation.time("gamma",as.numeric(undir_tdates$raw_serial_interval))

# Extract mean and coefficient of variation
mu <- fit.ti$mean
nu <- fit.ti$sd/fit.ti$mean


# Calculate reproductive number using exponential growth rate - method
r <- growth_rate 
R.tj <- (1 + r*mu*(nu^2))^(1/(nu^2)) ## a non value - due to negative exponential growth rate!! 

##### third method #####

daily_counts <- daily_counts %>% mutate(day=as.numeric(difftime(daily_counts$date,daily_counts$date[1],
                                                                units="days")))
df <- daily_counts[,2:3]

# calculate R0 using White and Pagano #
R0.ML.tj <- est.R0.ML(df$count,t=df$day,begin=1,end=32,GT=fit.si)

# calculate R0 using sequential bayesian approach: 
R0.SB.tj <- est.R0.SB(df$count,GT=fit.si)



