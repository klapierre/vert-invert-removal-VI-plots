library(ggplot2)
library(grid)
library(codyn)
library(lme4)
library(plyr)
library(dplyr)
library(tidyr)


setwd('C:\\Users\\Kim\\Dropbox\\konza projects\\VI plots\\data\\analysis')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

############################################################################
############################################################################
#read in plan
plan <- read.csv('Vert_Invert_plan.csv')%>%
  select(plot, NPK, exclose, insecticide)%>%
  mutate(trt=paste(NPK, exclose, insecticide, sep='_'))


biomass2009 <- read.csv('Vert_Invert_2009_biomass.csv')%>%
  mutate(year=2009)%>%
  select(year, plot, subplot, form, mass)%>%
  group_by(year, plot, subplot)%>%
  summarise(total=sum(mass))
biomass2010 <- read.csv('Vert_Invert_2010_biomass.csv')%>%
  mutate(year=2010)%>%
  select(year, plot, subplot, total)
biomass2011 <- read.csv('Vert_Invert_2011_biomass.csv')%>%
  mutate(year=2011)%>%
  select(year, plot, subplot, total)
biomass2012 <- read.csv('Vert_Invert_2012_biomass.csv')%>%
  mutate(year=2012)%>%
  select(year, plot, subplot, total)
biomass2013 <- read.csv('Vert_Invert_2013_biomass.csv')%>%
  mutate(year=2013)%>%
  select(year, plot, subplot, total)%>%
  filter(plot!='NA')
biomass2014 <- read.csv('Vert_Invert_2014_biomass.csv')%>%
  mutate(year=2014, total=(gram+forb+woody))%>%
  select(year, plot, subplot, total)

biomass <- rbind(biomass2009, biomass2010, biomass2011, biomass2012, biomass2013, biomass2014)%>%
  select(year, plot, subplot, total)%>%
  group_by(year, plot)%>%
  summarise(biomass=mean(total))%>%
  left_join(plan)



############################################################################
############################################################################

#mixed model
summary(lmer(biomass ~ NPK*insecticide*exclose*year + (1|plot), data=biomass))
summary(lmer(biomass ~ NPK*insecticide*year + (1|plot), data=biomass))

#insecticide effect alone
ggplot(data=barGraphStats(data=biomass, variable='biomass', byFactorNames=c('insecticide')), aes(x=insecticide, y=mean)) +
  geom_bar(stat='identity', colour='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  scale_x_discrete(breaks=c('x', 'insecticide'), labels=c('-inverts', 'control')) +
  ylab(expression(paste('Aboveground Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26))
#export at 500x600



#insecticide x nutrient effect
test <- biomass%>%
  group_by(NPK, insecticide)%>%
  summarise(mean=mean(biomass), N=length(biomass), sd=sd(biomass))%>%
  mutate(se=sd/N^.5, NPK_1=ifelse(NPK=='x', 0, 1))
ggplot(data=test, aes(x=NPK_1, y=mean, color=insecticide)) +
  geom_point(position=position_dodge(0.1), size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  geom_line(position=position_dodge(0.1), size=2) +
  scale_color_manual(values=c('#73409D', '#008000'), breaks=c('x', 'insecticide'), labels=c(' control', ' -invert')) +
  scale_x_continuous(breaks=c(0,1), labels=c('control', 'NPK')) +
  ylab(expression(paste('Aboveground Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1))
#export at 500x600














