### Analysis of Campanulaceae dN/dS data

#### Libraries ####
library(tidyverse)
library(ggtree)
library(treeio)
library(patchwork)
library(DHARMa)
library(car)

#### Read in data sets ####
atpsub <- read.codeml_mlc("BrysonPAML/InputFiles/atpfinsub_unroot.txt")
ndhsub <- read.codeml_mlc("BrysonPAML/InputFiles/ndhfinsub_unroot.txt")
petsub <- read.codeml_mlc("BrysonPAML/InputFiles/petfinsub_unroot.txt")
psasub <- read.codeml_mlc("BrysonPAML/InputFiles/psafinsub_unroot.txt")
psbsub <- read.codeml_mlc("BrysonPAML/InputFiles/psbfinsub_unroot.txt")
rplsub <- read.codeml_mlc("BrysonPAML/InputFiles/rplfinsub_unroot.txt")
rposub <- read.codeml_mlc("BrysonPAML/InputFiles/rpofinsub_unroot.txt")
rpssub <- read.codeml_mlc("BrysonPAML/InputFiles/rpsfinsub_unroot.txt")
clpPsub <- read.codeml_mlc("BrysonPAML/InputFiles/clpP1finsub_unroot.txt")
ycf2sub <- read.codeml_mlc("BrysonPAML/InputFiles/ycf2finsub_unroot.txt")
matKsub <- read.codeml_mlc("BrysonPAML/InputFiles/matKfinsub_unroot.txt")
ccsAsub <- read.codeml_mlc("BrysonPAML/InputFiles/ccsAfinsub_unroot.txt")
cemAsub <- read.codeml_mlc("BrysonPAML/InputFiles/cemAfinsub_unroot.txt")
pafsub <- read.codeml_mlc("BrysonPAML/InputFiles/paffinsub_unroot.txt")

### Pull out tables of dN/dS and other statistics, add gene family column, then concatenate
#Look at file to figure out how to pull out table of relevant info
str(atpsub)

# Pull out relevant table from atp datafile
atpsubtbl <- atpsub@data

# Add genefam column to atp table
atpsubtbl <- atpsubtbl %>%
  mutate(genefam="atp")

# Do the same with ndh
ndhsubtbl <- ndhsub@data
ndhsubtbl <- ndhsubtbl %>%
  mutate(genefam="ndh")

# Do the same with the rest
petsubtbl <- petsub@data
petsubtbl <- petsubtbl %>%
  mutate(genefam="pet")

psasubtbl <- psasub@data
psasubtbl <- psasubtbl %>%
  mutate(genefam="psa")

psbsubtbl <- psbsub@data
psbsubtbl <- psbsubtbl %>%
  mutate(genefam="psb")

rplsubtbl <- rplsub@data
rplsubtbl <- rplsubtbl %>%
  mutate(genefam="rpl")

rposubtbl <- rposub@data
rposubtbl <- rposubtbl %>%
  mutate(genefam="rpo")

rpssubtbl <- rpssub@data
rpssubtbl <- rpssubtbl %>%
  mutate(genefam="rps")

ycf2subtbl <- ycf2sub@data
ycf2subtbl <- ycf2subtbl %>%
  mutate(genefam="ycf2")

clpPsubtbl <- clpPsub@data
clpPsubtbl <- clpPsubtbl %>%
  mutate(genefam="clpP")

matKsubtbl <- matKsub@data
matKsubtbl <- matKsubtbl %>%
  mutate(genefam="matK")

ccsAsubtbl <- ccsAsub@data
ccsAsubtbl <- ccsAsubtbl %>%
  mutate(genefam="ccsA")

cemAsubtbl <- cemAsub@data
cemAsubtbl <- cemAsubtbl %>%
  mutate(genefam="cemA")

pafsubtbl <- pafsub@data
pafsubtbl <- pafsubtbl %>%
  mutate(genefam="paf")

# Concatenate individual tables
total <- rbind(atpsubtbl, ndhsubtbl, petsubtbl, psasubtbl, psbsubtbl, rplsubtbl, 
               rposubtbl, rpssubtbl, ycf2subtbl, clpPsubtbl, matKsubtbl, ccsAsubtbl,
               cemAsubtbl, pafsubtbl)

### Read out total table and fix dN/dS values that are actually uncaculable due to dN or dS of 0
write.csv(total, file="total_012226.csv", quote=FALSE, row.names=FALSE)

### Read back in fixed file
total <- read_csv("total_012226.csv")

### Organize gene family names
total <- total %>%
  mutate(genefam=factor(genefam, levels=c("atp", "ndh", "paf", "pet", "psa", "psb", "rbcl", "rpl", "rps", 
                                          "rpo", "clpP", "ycf2", "cemA", "ccsA", "matK")))

### Add functional category
total <- total %>%
  mutate(category=ifelse(genefam=="rpl" | genefam=="rps" | genefam=="rpo", "Gene Regulation", ifelse(
    genefam=="clpP" | genefam=="ycf2", "Proteostasis", ifelse(genefam=="matK" |
                                                                genefam=="cemA" | genefam=="ccsA", "Other", "Photosynthesis"
    ))))

#### Graph Campanula americana Only ####
CampOnly <- filter(total, Label=="Camp_IL6" | Label=="Camp_VA73")

CampSub <- select(CampOnly, c(6:7, 11:13))

CampSubLong <- pivot_longer(CampSub, 1:2, names_to = "type", values_to = "rate")

CampSubLong <- CampSubLong %>%
  mutate(Label=factor(Label, levels=c("Camp_VA73", "Camp_IL6")))

CampSubLong <- CampSubLong %>%
  mutate(genefam=factor(genefam, levels=c("atp", "ndh", "paf", "pet", "psa", "psb", "rbcl", "rpl", "rps", 
                                          "rpo", "clpP", "ycf2", "cemA", "ccsA", "matK")))

CampSubLong <- CampSubLong %>%
  mutate(category=factor(category, levels=c("Photosynthesis", "Gene Regulation", "Proteostasis", "Other")))

App <- ggplot(data=filter(CampSubLong, Label=="Camp_VA73"), aes(x=genefam, y=rate, fill=type)) + 
  scale_fill_manual(breaks=c("dN", "dS"), 
                    values=c("White", "Black")) +
  geom_bar(stat="identity", color="black", width=0.7) +
  ylim(0,0.01) +
  xlab("Gene Family") + ylab("Substitution Rate") + labs(fill="Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(.~category, scales="free_x", space="free") +
  theme(strip.text = element_blank()) +
  theme(panel.spacing = unit(1, "lines"))
App

West <- ggplot(data=filter(CampSubLong, Label=="Camp_IL6"), aes(x=genefam, y=rate, fill=type)) + 
  scale_fill_manual(breaks=c("dN", "dS"), 
                    values=c("White", "Black")) +
  geom_bar(stat="identity", color="black", width=0.7) +
  ylim(0,0.01) +
  xlab("Gene Family") + ylab("Substitution Rate") + labs(fill="Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(.~category, scales="free_x", space="free") +
  theme(strip.text = element_blank()) +
  theme(panel.spacing = unit(1, "lines"))
West

App + West + plot_layout(guides="collect", axis_titles = "collect", axes="collect", nrow=2)


#### Graph Campanualaceae ####
totaltiponly <- filter(total, Label>40)

totaltiponly <- totaltiponly %>%
  mutate(genefam=factor(genefam, levels=c("atp", "ndh", "paf", "pet", "psa", "psb", "rbcl", "rpl", "rps", 
                                          "rpo", "clpP", "ycf2", "cemA", "ccsA", "matK")))

### dN/dS by gene family
ggplot(data=filter(totaltiponly, Label!="Camp_VA73" & Label!="Camp_IL6"),
       aes(x=genefam, y=as.numeric(dN_vs_dS_B), color=category)) + 
  geom_jitter(width=0.2, height=0, alpha=0.5) + 
  geom_point(data=filter(totaltiponly, Label=="Camp_VA73" | Label=="Camp_IL6"), aes(x=genefam, y=as.numeric(dN_vs_dS_B)), 
             shape="diamond", size=3, color="black") +
  scale_color_manual(breaks=c("Photosynthesis", "Gene Regulation", "Proteostasis", "Other"), values=c("Green", "Red", "Blue", "Orange")) +
  geom_hline(yintercept=1, linetype="dashed") +
  xlab("Gene Family") +
  ylab("dN/dS") +
  scale_y_continuous(limits=c(0,10), breaks=seq(0,10,by=1)) +
  theme_classic()

### dN/dS by Species
totaltiponly <- totaltiponly %>%
  mutate(Label=factor(Label, c("Nicotiana", "Helianthus", "N36084", "N35624", "N63738", "N63739",
                               "N63740", "JMont", "N63742", "N26203", "Trachelium", "N24732", "M12303",
                               "N36221", "805474", "Tperoliata", "Camp_VA73", "Camp_IL6")),
         category=factor(category, levels=c("Photosynthesis", "Gene Regulation", "Proteostasis", "Other")))

LabelNew <- c("N.tabacuum", "H.annuus", "C.serratus", "P.grandiflorus", "C.bhutanica", "C.lobatus", "W. marginata",
              "J.montana", "C.pallida", "C.takesimana", "T.caeruleum", "H.asiatica", "A.racemosa",
              "A.divaricata", "A.japonicum", "T.perfoliata", "C.americana_App", "C.americana_West")

ggplot(data=filter(totaltiponly), aes(x=Label, y=dN_vs_dS_B, color=category, shape=category)) + 
  geom_jitter(width=0.15, height=0, alpha=0.6) +
  scale_color_manual(breaks=c("Photosynthesis", "Gene Regulation", "Proteostasis", "Other"), 
                     values=c("Green", "Red", "Blue", "Orange")) +
  scale_shape_manual(breaks=c("Photosynthesis", "Gene Regulation", "Proteostasis", "Other"), 
                     values=c(15, 16, 17, 18)) +
  stat_summary(mapping=aes(group=Label), fun="median", geom="crossbar", width=0.5, color="black", size=0.3) +
  xlab("Species") +
  ylab("dN/dS") +
  scale_x_discrete(labels=LabelNew)+
  geom_hline(yintercept=1, linetype="dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face="italic")) +
  scale_y_continuous(limits=c(0,10), breaks=seq(0,10,by=1))

#### Statistics ####

### Run model looking at impact of gene functional category and species (Label) on dN/dS
dNdS_mod <- glmmTMB(log(dN_vs_dS_B)~Label+category+(1|category:genefam), data=totaltiponly)

### Check that model look like a good fit
plot(simulateResiduals(dNdS_mod))

### Check significance
Anova(dNdS_mod, type=3)

### Look at individual contrasts
summary(dNdS_mod)
