library(tidyverse)
library(lubridate)
library(incidence)
library(MASS)
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
library(ggplot2)

library(EpiEstim)
library(R0)
library(reshape2)

setwd("~/Desktop/IDD/Project_COVID")
# Load the data
data <- read_csv("~/Desktop/IDD/Project_COVID/singapore_cleaned.csv")

# Convert the date variable to a Date object
data$date_confirmation <- as.Date(data$date_confirmation, format = "%Y.%m.%d")
data$date_onset_symptoms <- as.Date(data$date_onset_symptoms, format = "%Y-%m-%d")

#### first method ####

# Calculate daily counts
daily_counts <- data %>%
  count(date_confirmation) %>%
  rename(date = date_confirmation, count = n)

sg_plotinc <- daily_counts

N <- 5000000 # population at risk
C <- sum(daily_counts$count) # total cases

total_i <- 0
for(i in (N-C+1):(N-1)){
  total_i <- 1/i + total_i
}

R01.sg <- ((N-1)/C)*total_i

#### second method ####

# calculate the exponential growth rate using a regression model
growth_rate <- lm(log(count) ~ decimal_date(date), data = daily_counts)$coefficients[2]

# calculating serial interval 
# all serial interval calculations from line list data are based on Tindale et al., 2020

# load data with symptom onset 
data_sg <- read_csv("~/Desktop/IDD/Project_COVID/COVID-19_Singapore_data_revised.csv", col_types = list(presumed_infected_date = col_datetime()))
#table(spdata$`Related cases`)
data_sg$`Related cases`[which(data_sg$`Related cases` == "\n")] <- NA

data_sg <- rename(data_sg, related_cases = starts_with("Related"),
                 cluster_links = "Cluster links")
data_sg <- filter(data_sg, !is.na(date_onset_symptoms)) 

# we assume that in each pair, the person who exhibited symptoms first is considered the one who transmitted the infection.
nodes.sg <- data_sg$CaseID

# separating the contacts each case had into columns
data_sg <- separate(data_sg,col = related_cases,into = paste("contactID", 1:7, sep = "_"),
                   fill = "right")
data_sg <- mutate(data_sg, 
                 contactID_1 = as.numeric(contactID_1),
                 contactID_2 = as.numeric(contactID_2),
                 contactID_3 = as.numeric(contactID_3),
                 contactID_4 = as.numeric(contactID_4),
                 contactID_5 = as.numeric(contactID_5),
                 contactID_6 = as.numeric(contactID_6),
                 contactID_7 = as.numeric(contactID_7))

edges.sg <- data_sg[,1:8] # selects the contact columns in the dataset
edges.sg <- filter(edges.sg, !is.na(edges.sg$contactID_1)) 

SG_edges = data.frame(from=3,to=4) # create dummy edge in data frame

# fill the data frame with 'edges' as contacts for each case, and 'node' as each case
for (n in 1:nrow(edges.sg)) {
  for (k in 2:ncol(edges.sg)) { 
    if (!is.na(edges.sg[n,k])) {
      SG_edges=rbind(SG_edges, c(edges.sg[[n,k]],edges.sg[[n,1]])) 
    }  
  }
}

SG_edges=SG_edges[-1,] # removing the initial dummy edge

# remove duplicates and ensures edges aren't repeated
sg_undir=data.frame(from = pmin(SG_edges[,1],SG_edges[,2]),  
                 to = pmax(SG_edges[,1], SG_edges[,2]))
sg_undir = unique(sg_undir)
sg_undir = sg_undir[-which(sg_undir[,1]==sg_undir[,2]),]

data_sg_sympt <- data_sg[,c(1,22)] # selecting CaseID and date_onset_symptoms

# cleans up the data frame: rename columns, add info on symptom onset data for both 'from' and 'to' nodes
names(data_sg_sympt) <- str_replace(names(data_sg_sympt), "CaseID", "from")
sg_undir_dates <- left_join(sg_undir, data_sg_sympt, by = "from")
names(sg_undir_dates) <- str_replace(names(sg_undir_dates), "date_onset_symptoms", "from_sympt_date")
names(data_sg_sympt) <- str_replace(names(data_sg_sympt), "from", "to")
sg_undir_dates <- left_join(sg_undir_dates, data_sg_sympt, by = "to")
names(sg_undir_dates) <- str_replace(names(sg_undir_dates), "date_onset_symptoms", "to_sympt_date")
sg_undir_dates <- mutate(sg_undir_dates, earliest_sympt_onset = pmin(to_sympt_date, from_sympt_date, na.rm = T), 
                      raw_serial_interval = to_sympt_date - from_sympt_date,  
                      abs_serial_interval = abs(raw_serial_interval)) # serial interval calculation

# groups positive and negative serial interval values
pos <- filter(sg_undir_dates, raw_serial_interval >= 0)
neg <- filter(sg_undir_dates, raw_serial_interval < 0)
onlyone <- filter(sg_undir_dates, is.na(raw_serial_interval))

# renaming columns
names(neg)[1] <- "to"
names(neg)[2] <- "from"
names(neg)[3] <- "to_sympt_date"
names(neg)[4] <- "from_sympt_date"

sg_undir_dates <- bind_rows(pos, neg, onlyone)
sg_undir_dates$pto <- str_pad(sg_undir_dates$to, width = 2, side = "left", pad = "0")
sg_undir_dates$pfrom <- str_pad(sg_undir_dates$from, width = 2, side = "left", pad = "0")
sg_undir_dates <- mutate(sg_undir_dates, pairID = factor(paste("case", pfrom, "-", "case", pto, sep = "")))
rm(pos, neg, onlyone)

# calculating incubation periods
data_sg$minIncTimes <- data_sg$date_onset_symptoms - data_sg$end_source
data_sg$maxIncTimes <- data_sg$date_onset_symptoms - data_sg$start_source
data_sg$maxIncTimes = pmax(3, data_sg$maxIncTimes)
data_sg$minIncTimes = pmax(1, data_sg$minIncTimes)

# cleaning and storing the final dataset

data_sgclean <- data.frame(infector = sg_undir_dates$from, infectee = sg_undir_dates$to, serial.interval = sg_undir_dates$abs_serial_interval, inc.infector.min = data_sg$minIncTimes[match(sg_undir_dates$from,data_sg$CaseID)], inc.infector.max = data_sg$maxIncTimes[match(sg_undir_dates$from,data_sg$CaseID)], inc.infectee.min =  data_sg$minIncTimes[match(sg_undir_dates$to,data_sg$CaseID)], inc.infectee.max = data_sg$maxIncTimes[match(sg_undir_dates$to,data_sg$CaseID)])
data_sgclean = data_sgclean[!is.na(data_sgclean$serial.interval),] # filter out NAs
data_sgclean$serial.interval = as.numeric(data_sgclean$serial.interval)

#Pivot the to/from dates of symptom onset column to a long format, so that can make a legend based on this variable
undir_dotplot <- pivot_longer(sg_undir_dates, 
                              cols = contains("sympt_date"),
                              names_to = "pair_member",
                              values_to = "onset_date")

#Let's rename the values so it makes more sense in the legend
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "from_sympt_date", replacement = "Presumed infector")
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "to_sympt_date", replacement = "Presumed infectee")

#Make the Cleveland dotplot
s <- ggplot(undir_dotplot, aes(y = reorder(pairID, earliest_sympt_onset))) +
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
  labs(title = "Serial intervals of possible case pairs from Singapore, over time",
       x = "Date of symptom onset",
       y = "Case pairs")
s

# fitting a gamma distribution to the serial interval data
fit.si <- generation.time("gamma",data_sgclean$serial.interval)

# store mean and coefficient of variation
# method of calculation follows from Nishiura et al., 2009
mu <- fit.si$mean
nu <- fit.si$sd/fit.si$mean

# calculate reproductive number using exponential growth rate - method
r <- growth_rate
R.exp.sg <- (1 + r*mu*(nu^2))^(1/(nu^2))
R.exp.sg

##### third method #####

# defining the MLE method function

# set an arbitrary range of R0 values to optimize from
range = c(0.01,50)

# check.incid() formats incidence data for further estimation
epid <- check.incid(df$count, t=df$day, time.step=1)
# storing number of days epidemic has lasted 
begin <- 0
end <- max(epid$t)

# converting dates to number of days since first case
daily_counts <- daily_counts %>% mutate(day=as.numeric(difftime(daily_counts$date,daily_counts$date[1],
                                                                units="days")))
df <- daily_counts[,2:3]

# estimating R0 using the maximum likelihood function as developed by White and Pagano, 2009 and implemented by Boelle and Obadia, 2011.

# inside the function below: 
## fit.epid() computes the Poisson log-likelihood of the epidemic and predicts the curve
# res.R <- optimize(f=fit.epid, log(range),GT=GT,epid=epid,import=NULL,maximum=TRUE) 
# pred <- fit.epid(res.R$maximum,epid,GT,import=import,pred=TRUE)
# returns the estimate of the reproduction number using MLE and its confidence interval

# calculate R0 using White and Pagano 
R0.ML.sg <- est.R0.ML(df$count,t=df$day,begin=begin,end=end,GT=fit.si)
R0.ML.sg

# calculate R0 using sequential bayesian approach: 
R0.SB.sg <- est.R0.SB(df$count,GT=fit.si)

#### sensitivity analysis ####

tmp <- sensitivity.analysis(incid=df$count,GT.type="gamma", GT.mean=seq(0.9,3.9,0.5), 
                            GT.sd.range=4.49, begin=1, end=20, est.method="EG",sa.type="GT")
plot(x=tmp[,"GT.Mean"], xlab="mean GT (days)", y=tmp[,"R"], ylim=c(0,2), ylab="R0 (95% CI)", 
     type="p", pch=19, col="black", main="Sensitivity of R0 to mean GT")
arrows(x0=as.numeric(tmp[,"GT.Mean"]), y0=as.numeric(tmp[,"CI.lower"]), 
       y1=as.numeric(tmp[,"CI.upper"]), angle=90, code=3, col="black", length=0.05)

tmp <- sensitivity.analysis(incid=df$count,GT.type="gamma", GT.mean=fit.si$mean, 
                            GT.sd.range=seq(2,6,0.5), begin=1, end=20, est.method="EG",sa.type="GT")
plot(x=tmp[,"GT.Mean"], xlab="mean GT (days)", y=tmp[,"R"], ylim=c(0,2), ylab="R0 (95% CI)", 
     type="p", pch=19, col="black", main="Sensitivity of R0 to mean GT")
arrows(x0=as.numeric(tmp[,"GT.Mean"]), y0=as.numeric(tmp[,"CI.lower"]), 
       y1=as.numeric(tmp[,"CI.upper"]), angle=90, code=3, col="black", length=0.05)