library(grid)
library(codyn)
library(lme4)
library(MASS)
library(vegan)
library(tidyverse)


setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\VI plots\\data\\analysis')


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
  select(plot, NPK, exclose, insecticide)

#read in data
cover2009 <- read.csv('Vert_Invert_2009_cover.csv')
cover2010 <- read.csv('Vert_Invert_2010_cover.csv')
cover2011 <- read.csv('Vert_Invert_2011_cover.csv')
cover2012 <- read.csv('Vert_Invert_2012_cover.csv')
cover2013 <- read.csv('Vert_Invert_2013_cover.csv')
cover2014 <- read.csv('Vert_Invert_2014_cover.csv')
cover2015 <- read.csv('Vert_Invert_2015_cover.csv')
cover2016 <- read.csv('Vert_Invert_2016_cover.csv')
cover2017 <- read.csv('Vert_Invert_2017_cover.csv')
cover2018 <- read.csv('Vert_Invert_2018_cover.csv')
cover2019 <- read.csv('Vert_Invert_2019_cover.csv')

cover <- rbind(cover2009, cover2010, cover2011, cover2012, cover2013, cover2014, cover2015, cover2016, cover2017, cover2018, cover2019)%>%
  select(plot, sppnum, cover, season, year)%>%
  group_by(year, plot, sppnum)%>%
  summarise(cover=max(cover))%>%
  left_join(plan)%>%
  filter(sppnum!=0, sppnum!=999)

coverSum <- cover%>%
  group_by(year, plot)%>%
  summarise(cover_sum=sum(cover))

coverRel <- cover%>%
  left_join(coverSum)%>%
  mutate(rel_cover=cover/cover_sum)%>%
  select(-cover_sum, -cover)

coverRelWide <- coverRel%>%
  spread(key=sppnum, value=rel_cover, fill=0)

############################################################################
############################################################################

###MDS
mds<- read.csv('VI_mds_2017.csv')

ggplot(data=subset(mds,year==2013), aes(x=X1, y=X2, shape=interaction(NPK, insecticide), colour=interaction(NPK, insecticide))) +
  geom_point(size=8) +
  xlab('MDS 1') +
  ylab('MDS 2') +
  scale_colour_manual(values=c("#FF9900", "#FF9900", "#009900", "#009900"), labels=c('+NPK-inverts    ', '-inverts    ', '+NPK    ', 'control    '), guide=guide_legend(reverse=T)) +
  scale_shape_manual(values=c(16, 17, 16, 17), labels=c('+NPK-inverts    ', '-inverts    ', '+NPK    ', 'control    '), guide=guide_legend(reverse=T)) +
  theme(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=30), axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=30), legend.title=element_blank(), legend.text=element_text(size=30),  legend.position='bottom', legend.margin=unit(2, "lines"))
#export at 1400x1500




###turnover
turnover <- turnover(df = coverRel, time.var = 'year', species.var = 'sppnum', abundance.var = 'rel_cover', replicate.var = 'plot')

#change plot to numeric to merge with plan
turnover$plot <- as.integer(turnover$plot)

#merge with plan
turnover <- turnover%>%
  left_join(plan)

#mixed model
summary(lmer(total ~ NPK*insecticide*exclose + (1|plot), data=turnover))

#plot turnover
#all interactions
ggplot(data=barGraphStats(data=turnover, variable='total', byFactorNames=c('year', 'NPK', 'insecticide', 'exclose')), aes(x=year, y=mean, colour=interaction(NPK, insecticide, exclose))) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_point() +
  geom_line()

#only NPK, significant interaction
ggplot(data=barGraphStats(data=turnover, variable='total', byFactorNames=c('year', 'NPK')), aes(x=year, y=mean, colour=NPK)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_point() +
  geom_line()


#parsing into apperance and disappearance
#calculate appearance
appearance <- turnover(df = coverRel, time.var = 'year', species.var = 'sppnum', abundance.var = 'rel_cover', replicate.var = 'plot', metric='appearance')

#change plot to numeric to merge with plan
appearance$plot <- as.integer(appearance$plot)

#merge with plan
appearance <- appearance%>%
  left_join(plan)

#mixed model
summary(lmer(appearance ~ NPK*insecticide*exclose + (1|plot), data=appearance))

#calculate disappearance
disappearance <- turnover(df = coverRel, time.var = 'year', species.var = 'sppnum', abundance.var = 'rel_cover', replicate.var = 'plot', metric='disappearance')

#change plot to numeric to merge with plan
disappearance$plot <- as.integer(disappearance$plot)

#merge with plan
disappearance <- disappearance%>%
  left_join(plan)

#mixed model
summary(lmer(disappearance ~ NPK*insecticide*exclose + (1|plot), data=disappearance))

#only NPK, significant interaction
ggplot(data=barGraphStats(data=disappearance, variable='disappearance', byFactorNames=c('year', 'NPK')), aes(x=year, y=mean, colour=NPK)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_point() +
  geom_line()



###rank shifts
rankShift <- rank_shift(df=coverRel, time.var='year', species.var='sppnum', abundance.var='rel_cover', replicate.var='plot')

#Create a column with the final year from the returned time.var_pair
rankShift$year <- as.numeric(substr(rankShift$year_pair, 6, 9))

#change plot to numeric to merge with plan
rankShift$plot <- as.integer(rankShift$plot)

#merge with plan
rankShift <- rankShift%>%
  left_join(plan)

#mixed model
summary(lmer(MRS ~ NPK*insecticide*exclose + (1|plot), data=rankShift))
summary(lmer(MRS ~ NPK*insecticide + (1|plot), data=rankShift))





###rate of change
rateChange <- rate_change_interval(coverRel, time.var='year', species.var='sppnum', abundance.var='rel_cover', replicate.var='plot')

#change plot to numeric to merge with plan
rateChange$plot <- as.integer(rateChange$plot)

#merge with plan
rateChange <- rateChange%>%
  left_join(plan)

#mixed model
summary(lmer(distance ~ NPK*insecticide*exclose + (1|interval), data=rateChange))
summary(lmer(distance ~ NPK*insecticide + (1|interval), data=subset(rateChange, exclose=='x')))

#only NPK, significant interaction
ggplot(data=barGraphStats(data=subset(rateChange, exclose=='x'), variable='distance', byFactorNames=c('interval', 'NPK', 'insecticide', 'exclose')), aes(x=interval, y=mean, colour=interaction(NPK, insecticide, exclose), shape=interaction(NPK, insecticide, exclose))) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.05) +
  geom_point(size=6) +
  geom_smooth(method=lm, se=F) +
  xlab('Time Interval') +
  ylab('Plant Community Change') +
  scale_colour_manual(values=c("#FF9900", "#FF9900", "#009900", "#009900"), labels=c(' +NPK-inverts', ' -inverts', ' +NPK', ' control'), guide=guide_legend(reverse=T)) +
  scale_shape_manual(values=c(16, 17, 16, 17), labels=c(' +NPK-inverts', ' -inverts', ' +NPK', ' control'), guide=guide_legend(reverse=T)) +
  theme(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=30), axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=30), legend.title=element_blank(), legend.text=element_text(size=30), legend.justification=c(0,1), legend.position=c(0,1))
#export at 1400x1400




###temporal community stability
###rate of change
stability <- community_stability(coverRel, time.var='year', abundance.var='rel_cover', replicate.var='plot')

#change plot to numeric to merge with plan
stability$plot <- as.integer(stability$plot)

#merge with plan
stability <- stability%>%
  left_join(plan)

#mixed model
summary(glm(stability ~ NPK*insecticide*exclose, data=stability))




###variance ratio
varianceRatio <- variance_ratio(df=coverRel, time.var='year', species.var='sppnum', abundance.var='rel_cover', replicate.var='plot', bootnumber=10, average.replicates = FALSE)

#change plot to numeric to merge with plan
varianceRatio$plot <- as.integer(varianceRatio$plot)

#merge with plan
varianceRatio <- varianceRatio%>%
  left_join(plan)

#mixed model
summary(glm(VR ~ NPK*insecticide*exclose, data=varianceRatio))




###rank clocks
#with dominant grass (schiz) and panicum
rankSpp <- coverRel%>%
  # #remove vert removal plots
  # filter(exclose!='caged')%>%
  #select only the relevant rare species
  filter(sppnum==111|sppnum==46|sppnum==55|sppnum==15|sppnum==121|sppnum==3)%>%
  mutate(spp=ifelse(sppnum==3, 'SCSC', ifelse(sppnum==15, 'PAVI', ifelse(sppnum=='46', 'AMPS', ifelse(sppnum==55, 'ASVI', ifelse(sppnum==111, 'PHPU', 'SAAZ'))))))

ggplot(rankSpp, aes(x=year, y=rel_cover, color=spp)) + 
  geom_line(size = 2) +
  coord_polar() +
  facet_wrap(~plot, ncol=6)


#without dominant grass (schiz) and panicum
rareSpp <- coverRel%>%
  # #remove vert removal plots
  # filter(exclose!='caged')%>%
  #select only the relevant rare species
  filter(sppnum==111|sppnum==46|sppnum==55|sppnum==121)%>%
  mutate(spp=ifelse(sppnum==3, 'SCSC', ifelse(sppnum==15, 'PAVI', ifelse(sppnum=='46', 'AMPS', ifelse(sppnum==55, 'ASVI', ifelse(sppnum==111, 'PHPU', 'SAAZ'))))))

ggplot(rareSpp, aes(x=year, y=rel_cover, color=spp)) + 
  geom_line(size = 2) +
  coord_polar() +
  facet_wrap(~plot, ncol=6)


#only plots that are intersting; without dominant grass (schiz) and panicum
rareSpp <- coverRel%>%
  # #remove vert removal plots
  # filter(exclose!='caged')%>%
  #select only the relevant rare species
  filter(sppnum==111|sppnum==46|sppnum==55|sppnum==121|sppnum==15)%>%
  filter(plot==1|plot==2|plot==10|plot==12|plot==17|plot==20|plot==8|plot==4|plot==3)%>%
  mutate(spp=ifelse(sppnum==3, 'S. scoparium  ', ifelse(sppnum==15, 'P. virgatum  ', ifelse(sppnum=='46', 'A. psilostachya  ', ifelse(sppnum==55, 'A. verticillata  ', ifelse(sppnum==111, 'P. pumilis  ', 'S. azura  '))))))%>%
  mutate(type=ifelse(plot==3, 'control', ifelse(plot==4, '+NPK', ifelse(plot==8, '-inverts', ifelse(plot==1, '+NPK-inverts #1*', ifelse(plot==2, '+NPK-inverts #2', ifelse(plot==10, '+NPK-inverts #3*', ifelse(plot==12, '+NPK-inverts #4', ifelse(plot==17, '+NPK-inverts #5', '+NPK-inverts #6*')))))))))

rareSpp$order <- factor(rareSpp$type, levels = c('control', '-inverts', '+NPK', '+NPK-inverts #1*', '+NPK-inverts #2', '+NPK-inverts #3*', '+NPK-inverts #4', '+NPK-inverts #5', '+NPK-inverts #6*'))
type <- rareSpp$type
ggplot(rareSpp, aes(x=year, y=rel_cover, color=spp)) + 
  geom_line(size = 2) +
  coord_polar() +
  facet_wrap(~order, ncol=3) +
  xlab('Year') +
  ylab('Relative Abundance (%)') +
  theme(strip.text = element_text(size=24), 
        legend.position='bottom', legend.margin=unit(2, "lines"),
        axis.title.x=element_text(size=24, vjust=-0.35, margin=margin(t=15)),
        axis.text.x=element_text(size=20),
        axis.title.y=element_text(size=24, angle=90, vjust=0.5, margin=margin(r=15)),
        axis.text.y=element_text(size=20))
#export at 1600x1500



