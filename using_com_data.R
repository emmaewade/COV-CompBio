setwd("C:/School/compbio/project/Raw_Data/Combining-Different-Gene-Expresssion-Tables/")
com_analysis = read.table("OUTPUT.csv", sep = "\t", h =T)
colnames(com_analysis) = c("covid_genes", "covid_p", "covid_log", "sars_genes", "sars_p", "sars_log", "mers_genes", "mers_p",
                           "mers_log", "h1n1_genes", "h1n1_p", "h1n1_log", "ebola_genes", "ebola_p", "ebola_log")

neg_threshold = log2(1/1.5)
pos_threshold = log2(1.5)

covid = read.table("COVID-19.csv", sep = "\t", h=T)
covid$virus = "covid"
covid$DEG = NA
covid$DEG[which(covid$logFC <= neg_threshold)] = "up"
covid$DEG[which(covid$logFC >= pos_threshold) ] = "down"

sars = read.table("SARS-COV.csv", sep = "\t", h=T)
sars$virus = "sars"
sars$DEG = NA
sars$DEG[which(sars$logFC <= neg_threshold)] = "up"
sars$DEG[which(sars$logFC >= pos_threshold)] = "down"

mers = read.table("MERS-CoV.csv", sep = "\t", h=T)
mers$virus = "mers"
mers$DEG = NA
mers$DEG[which(mers$logFC >= pos_threshold)] = "up"
mers$DEG[which(mers$logFC <= neg_threshold)] = "down"

h1n1 = read.table("H1N1.csv", sep = "\t", h=T)
h1n1$virus = "h1n1"
h1n1$DEG = NA
h1n1$DEG[which(h1n1$logFC <= neg_threshold)] = "up"
h1n1$DEG[which(h1n1$logFC >= pos_threshold)] = "down"

ebola = read.table("Ebola.csv", sep = "\t", h=T)
ebola$virus = "ebola"
ebola$DEG = NA
ebola$DEG[which(ebola$logFC <= neg_threshold)] = "up"
ebola$DEG[which(ebola$logFC >= pos_threshold) ] = "down"

com = rbind(covid, sars, mers, h1n1, ebola)

deg_up = com[which(com$DEG == "up"), ]
deg_down = com[which(com$DEG == "down"), ]
deg_all = rbind(deg_up, deg_down)


##venn diagram##
X_up = list(covid = deg_up$Gene.symbol[which(deg_up$virus == "covid")], sars = deg_up$Gene.symbol[which(deg_up$virus == "sars")], mers = deg_up$Gene.symbol[which(deg_up$virus == "mers")], h1n1 = deg_up$Gene.symbol[which(deg_up$virus == "h1n1")], ebola = deg_up$Gene.symbol[which(deg_up$virus == "ebola")])
X_down = list(covid = deg_down$Gene.symbol[which(deg_down$virus == "covid")], sars = deg_down$Gene.symbol[which(deg_down$virus == "sars")], mers = deg_down$Gene.symbol[which(deg_down$virus == "mers")], h1n1 = deg_down$Gene.symbol[which(deg_down$virus == "h1n1")], ebola = deg_down$Gene.symbol[which(deg_down$virus == "ebola")])
X_all = list(covid = deg_all$Gene.symbol[which(deg_all$virus == "covid")], sars = deg_all$Gene.symbol[which(deg_all$virus == "sars")], mers = deg_all$Gene.symbol[which(deg_all$virus == "mers")], h1n1 = deg_all$Gene.symbol[which(deg_all$virus == "h1n1")], ebola = deg_all$Gene.symbol[which(deg_all$virus == "ebola")])

##higher in infected####
library(VennDiagram)
grid.newpage()
venn_object <- venn.diagram(X_up, filename = NULL, fill = c("pink", "green", "orange", "yellow", "cyan"))
grid.draw(venn_object)

##lower in infected####
library(VennDiagram)
grid.newpage()
venn_object <- venn.diagram(X_down, filename = NULL, fill = c("pink", "green", "orange", "yellow", "cyan"))
grid.draw(venn_object)

##all degs####
library(VennDiagram)
grid.newpage()
venn_object <- venn.diagram(X_all, filename = NULL, fill = c("pink", "green", "orange", "yellow", "cyan"))
grid.draw(venn_object)

##covid only###
#com$deglabel = NA
#com$deglabel[which(com$Gene.symbol == "SAA2" | com$Gene.symbol == "CCL20" | com$Gene.symbol == "IL8")] = com$Gene.symbol

cov_unique = com[which(com$Gene.symbol %in% setdiff(com$Gene.symbol[which(com$virus == "covid") ], com$Gene.symbol[which(com$virus != "covid") ])),]

library(ggplot2)
ggplot(cov_unique, aes(logFC, P.Value, color = DEG)) +
  geom_point()+
  geom_label(label = cov_unique$Gene.symbol, nudge_x=0.0000000001, nudge_y=0.0000000001, size = 2)+
  geom_vline(xintercept = c(neg_threshold, pos_threshold), linetype = 'dashed', color = "black", size = 0.5)+
  theme_bw()+
  theme(text = element_text (size = 15))#+
  #ylim(-.0000000001,.0000000005)+xlim(-2.2,-1.4)

##all degs####
ggplot(com, aes(logFC, P.Value, color = virus)) +
  geom_point()+
  #geom_label(label = com$Gene.symbol[com$virus == "covid" & com$P.Value < 0.001])+
  geom_vline(xintercept = c(neg_threshold, pos_threshold), linetype = 'dashed', color = "black", size = 0.5)+
  theme_bw()+
  theme(text = element_text (size = 15))


