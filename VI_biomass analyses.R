library(grid)
library(PerformanceAnalytics)
library(nlme)
library(emmeans)
library(tidyverse)

#NOTE: treatments run from 2009-2018, treatments ceased starting in 2019 to follow recovery
#NOTE: biomass data has been through QA/QC checks through 2020 and outliers are confirmed to be true values

#set path
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\VI plots\\data\\analysis')


#set options
options(contrasts=c('contr.sum','contr.poly'))

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

#homemade functions
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
  spread(key=form, value=mass, fill=0)%>%
  mutate(gram2=(andro+schiz+sorg+gram), forb2=(forb+aster), pdead=0, total=(gram2+forb2+woody))%>%
  select(year, plot, subplot, gram2, forb2, woody, pdead, total)%>%
  rename(gram=gram2, forb=forb2)
biomass2010 <- read.csv('Vert_Invert_2010_biomass.csv')%>%
  mutate(year=2010)%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)
biomass2011 <- read.csv('Vert_Invert_2011_biomass.csv')%>%
  mutate(year=2011)%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)
biomass2012 <- read.csv('Vert_Invert_2012_biomass.csv')%>%
  mutate(year=2012)%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)
biomass2013 <- read.csv('Vert_Invert_2013_biomass.csv')%>%
  mutate(year=2013)%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)%>%
  filter(plot!='NA')
biomass2014 <- read.csv('Vert_Invert_2014_biomass.csv')%>%
  mutate(year=2014, total=(gram+forb+woody))%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)
biomass2015 <- read.csv('Vert_Invert_2015_biomass.csv')%>%
  mutate(year=2015)%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)%>%
  filter(plot!='NA')
biomass2016 <- read.csv('Vert_Invert_2016_biomass.csv')%>%
  mutate(year=2016, total=(gram+forb+woody))%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)%>%
  filter(plot!='NA')
biomass2017 <- read.csv('Vert_Invert_2017_biomass.csv')%>%
  mutate(year=2017, total=(gram+forb+woody))%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)%>%
  filter(plot!='NA')
biomass2018 <- read.csv('Vert_Invert_2018_biomass.csv')%>%
  mutate(year=2018, total=(gram+forb+woody))%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)%>%
  filter(plot!='NA')
biomass2019 <- read.csv('Vert_Invert_2019_biomass.csv')%>%
  filter(plot!='NA')%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  mutate(year=2019, total=(gram+forb+woody))%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)
biomass2020 <- read.csv('Vert_Invert_2020_biomass.csv')%>%
  mutate(drop=ifelse(plot==24&subplot==3, 1, 0))%>% #remove missing sample
  filter(drop!=1)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  mutate(year=2020, total=(gram+forb+woody))%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)

biomass <- rbind(biomass2009, biomass2010, biomass2011, biomass2012, biomass2013, biomass2014, biomass2015, biomass2016, biomass2017, biomass2018, biomass2019, biomass2020)%>%
  select(year, plot, subplot, gram, forb, woody, pdead, total)%>%
  left_join(plan)

#checking for outliers
dataVis <- biomass%>%
  select(gram, forb, woody, pdead, total) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)

chart.Correlation(subset(biomass, trt=='x_x_x')%>%
  select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='x_x_insecticide')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='x_caged_x')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='x_caged_insecticide')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='NPK_x_x')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='NPK_x_insecticide')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='NPK_caged_x')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)
chart.Correlation(subset(biomass, trt=='NPK_caged_insecticide')%>%
                    select(gram, forb, woody, pdead, total), histogram=T, pch=19)

#generating mean across clip strips and scaling to g/m2
biomassMean <- biomass%>%
  group_by(year, plot)%>%
  summarise(total=10*mean(total), gram=10*mean(gram), forb=10*mean(forb), woody=10*mean(woody), pdead=10*mean(pdead))%>%
  left_join(plan)





############################################################################
############################################################################

#mixed models

#total biomass
summary(biomassTrtYrs <- lme(total~NPK*insecticide*exclose,
                             data=subset(biomassMean, year<2019),
                             random=~1|plot,
                             correlation=corCompSymm(form=~year|plot), 
                             control=lmeControl(returnObject=T)))
anova.lme(biomassTrtYrs, type='sequential') 
emmeans(biomassTrtYrs, pairwise~NPK*insecticide*exclose, adjust="tukey")

#grass biomass
summary(gramTrtYrs <- lme(gram~NPK*insecticide*exclose,
                             data=subset(biomassMean, year<2019),
                             random=~1|plot,
                             correlation=corCompSymm(form=~year|plot), 
                             control=lmeControl(returnObject=T)))
anova.lme(gramTrtYrs, type='sequential') 
emmeans(gramTrtYrs, pairwise~NPK*insecticide*exclose, adjust="tukey")

#forb biomass
summary(forbTrtYrs <- lme(forb~NPK*insecticide*exclose,
                          data=subset(biomassMean, year<2019),
                          random=~1|plot,
                          correlation=corCompSymm(form=~year|plot), 
                          control=lmeControl(returnObject=T)))
anova.lme(forbTrtYrs, type='sequential') 
emmeans(forbTrtYrs, pairwise~NPK*insecticide*exclose, adjust="tukey")

#forb:grass ratio
summary(grassforbTrtYrs <- lme((forb/gram)~NPK*insecticide*exclose,
                               data=subset(biomassMean, year<2019),
                               random=~1|plot,
                               correlation=corCompSymm(form=~year|plot), 
                               control=lmeControl(returnObject=T)))
anova.lme(grassforbTrtYrs, type='sequential') 
emmeans(grassforbTrtYrs, pairwise~NPK*insecticide*exclose, adjust="tukey")


#bar graphs
ggplot(data=barGraphStats(data=subset(biomassMean, year<2019), variable="gram", byFactorNames=c("trt")), aes(x=trt, y=mean)) +
  geom_point(position=position_dodge(0.1), size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  geom_line(position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1))
#export at 500x600

ggplot(data=barGraphStats(data=subset(biomassMean, year<2019), variable="forb", byFactorNames=c("trt")), aes(x=trt, y=mean)) +
  geom_point(position=position_dodge(0.1), size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  geom_line(position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, angle=90, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1))
#export at 500x600


#time series
ggplot(data=barGraphStats(data=subset(biomassMean, year<2019), variable="gram", byFactorNames=c("year", "trt")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  geom_line(position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1))
#export at 500x600

ggplot(data=barGraphStats(data=subset(biomassMean, year<2019), variable="forb", byFactorNames=c("year", "trt")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1), size=2) +
  geom_line(position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1))
#export at 500x600


##### recovery #####
#time series
ggplot(data=barGraphStats(data=subset(biomassMean, year>2017), variable="gram", byFactorNames=c("year", "trt")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, size=2) +
  geom_line(position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Graminoid Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1))
#export at 500x600

ggplot(data=barGraphStats(data=subset(biomassMean, year>2017), variable="forb", byFactorNames=c("year", "trt")), aes(x=year, y=mean, color=trt)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, size=2) +
  geom_line(position=position_dodge(0.1), size=2) +
  ylab(expression(paste('Forb Biomass (g m'^'-2',')'))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30, angle=90), axis.title.y=element_text(size=30, vjust=1, margin=margin(r=15)), axis.text.y=element_text(size=26), legend.position=c(0, 1), legend.justification=c(0,1))
#export at 500x600