
library(tidyverse)
library(readr)

custom_order = c('Jlux_08945BO', 'Jlux_08946BO', 'Jlux_08947BO', 'Jlux_17668BF', 'Jlux_17669BF', 'Jlux_17670BF', #luxata
                 'Jmul_17901MT', 'Jmul_17902MT', 'Jmul_17903MT', 'Jmul_17909RS', 'Jmul_17912RS', 'Jmul_17915RS', 'Jspp_06607SU', 'Jspp_06616SU', 'Jspp_06618SU', #lineataN
                 'Jspp_06617SU', 'Jspp_06619SU', 'Jspp_06620SU', 'Jmul_03795AR', 'Jmul_03796ar', 'Jmul_03797AR', 'Jmul_08354SO', 'Jmul_08377SO', 'Jmul_08400SO', 'Jmul_06589SO', 'Jmul_06594SO', 'Jmul_06604SO', 'Jmul_08352SO','Jmul_03771SR', 'Jmul_03772SR', 'Jmul_03773SR', 'Jmul_03778CE', 'Jmul_03780CE', 'Jmul_03781CE', 'Jmul_03851CE', 'Jmul_03852CE', 'Jmul_03853CE',  'Jmul_08353CE',  #lineataW
                 'Jmul_08311LD', 'Jmul_08314LD', 'Jmul_08315LD', 'Jmul_08316Ld', 'Jmul_08317Ld', 'Jmul_08318Ld','Jlin_03760RC', 'Jlin_03761RC', 'Jlin_03762RC', 'Jlin_08228Rc', 'Jlin_08230Rc', 'Jlin_08231Rc','Jlin_03751LM', 'Jlin_03752LM', 'Jlin_03753LM',  'Jmul_06621PT','Jmul_06623PT', 'Jmul_06625PT',  #lineataE
                 'Jmul_08671LP', 'Jmul_08672LP', 'Jmul_08674LP', 'Jmul_09891LP', 'Jmul_09892LP', 'Jmul_09893LP','Jmul_08748LL', 'Jmul_08749LL', 'Jmul_08750LL','Jmul_08721LC', 'Jmul_08722LC', 'Jmul_08723LC', #darwini
                 'Jonc_03716T1', 'Jonc_03717T1', 'Jonc_03718T1', 'Jonc_03734T3', 'Jonc_03735T3', 'Jonc_03736T3', 'Jonc_03739T4', 'Jonc_03740T4', 'Jonc_03741T4','Jonc_03626OG', 'Jonc_03627RN', 'Jonc_03621RN','Jonc_03696RY', 'Jonc_03697RY', 'Jonc_03699RY') #onca


#For k=7

K7 <- read_table("data/plink.7.Q.csv", 
                  col_names = FALSE)
popmap<-read.table("popmap.txt")
K7<-cbind(K7,popmap)
long<-K7 %>% 
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5, X6, X7) ,
    names_to = 'Cluster',
    values_to = 'Probability')

long$V1 <- factor(long$V1, levels = custom_order)

K7<-ggplot(long,aes(x=V1,y=Probability,fill=Cluster))+geom_bar(stat='identity',position='stack', color="black")+theme_classic()+ 
  scale_fill_manual(
    values = c("#6CC7AB", "#6a59cdff","slateblue4","#B95251","#1c90b9ff","#4D4D4D","#EE9A00"))+
  labs(x=NULL)+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour='black'))

#For k=6

K6 <- read_table("data/plink.6.Q.csv", 
                 col_names = FALSE)
popmap<-read.table("popmap.txt")
K6<-cbind(K6,popmap)


long<-K6 %>% 
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5, X6) ,
    names_to = 'Cluster',
    values_to = 'Probability')

long$V1 <- factor(long$V1, levels = custom_order)

K6<-ggplot(long,aes(x=V1,y=Probability,fill=Cluster))+geom_bar(stat='identity',position='stack', color="black")+theme_classic()+ 
  scale_fill_manual(
    values = c("#1c90b9ff", "#6CC7AB",  "#B95251", "#4D4D4D","#EE9A00","#6a59cdff"))+
  labs(x=NULL)+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour='black'))

#For k=5

K5 <- read_table("data/plink.5.Q.csv", 
                 col_names = FALSE)
popmap<-read.table("popmap.txt")
K5<-cbind(K5,popmap)

long<-K5 %>% 
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5) ,
    names_to = 'Cluster',
    values_to = 'Probability')

long$V1 <- factor(long$V1, levels = custom_order)

K5<-ggplot(long,aes(x=V1,y=Probability,fill=Cluster))+geom_bar(stat='identity',position='stack', color="black")+theme_classic()+ 
  scale_fill_manual(
    values = c( "#B95251", "#6a59cdff","#EE9A00",  "#4D4D4D","#6CC7AB"))+
  labs(x=NULL)+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour='black'))
