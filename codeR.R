library(data.table)
library(sas7bdat)
library(sampling)
library(FactoMineR)
library(haven)
library(factoextra)
library(dplyr)
library(ade4)
library(descr)
library(crosstable)
library(rpart)
library(rpart.plot)


## Chargement des données


### Chargement du fichier caractéristiques : 
setwd("/Users/PC/Documents/INSEA/3A/2nd part/MOUSSANIF/projet")
caract <- read.sas7bdat("caract.sas7bdat")
caract<-as.data.table(caract)
caract

### Chargement du fichier marges : 

marge<-read.sas7bdat("marge_insee.sas7bdat")
marge
### Chargement du fichier réponses : 

reponse<-read.sas7bdat("Reponse.sas7bdat")
reponse<-as.data.table(reponse)
View(reponse)
## Modifications 

#Renommer la colonne _CODEPAN_ pour la matcher avec la colonne _hh_  :

reponse<-reponse[,":="(hh=CODEPAN)][,":="(CODEPAN=NULL)]
reponse
reponse<-reponse[!duplicated(reponse[,"hh"]),]
merge<-caract[reponse, on="hh", nomatch=0]
View(merge)

# Extrapolation de l’enquête

## Les critères d'extrapolation : *CSPCHEF, AGECHEF, REGION, NBPERS, TU*

A = merge[,c(3,4,2,5,6)]
View(A)
## Tableau disjonctif - éclatement des modalités :

Y<-data.frame(sapply(A,function(x) data.table(model.matrix(~x-1,data =A))[,]))
View(Y)

#Renommer les colonnes pour correspondre à _LIBELLE_ dans marge : 

names(Y)[1:4]<-c("<35", ">65",  "35-50" , "50-65")
names(Y)[5:9]<-c("centre EST","Centre Ouest","Nord", "R Parisienne", "Sud")
names(Y)[10:17]<-c("AGRI","CADM","CADS",  "COMM", "EMPLOY","INACTIF",  "OUVRIER","RETRAIT")
names(Y)[18:22]<-c("1", "2", "3", "4", "5 & PLUS")
names(Y)[23:27]<-c("< 20 000 individus", ">200 000 individus hors Agg. Parisienne","20 000-200 000 individus", "Agg. Parisienne" ,"Rurale"  )

View(Y)

#Suppression de la variable _ENF_ : 

marge<-marge[-c(1,2,3,4),]
marge$LIBELLE
View(marge)

#Réordonner les totaux dans _marge_ avec les libellés de _Y_, à travers une jointure : 

lib<-names(Y)
View(lib)
lib<-data.table("LIBELLE" = lib)
univers<-inner_join(lib, marge, by="LIBELLE")
View(univers)

#Taille de l'échantillon : 

n = nrow(merge)
n

#Taille de la population :

N = sum(univers$MAR[1:4])
N

#Totaux et Calcul des probabilités d'inclusion : 

totaux<-t(as.matrix(univers$MAR))
totaux

piks = rep(n/N, n)
head(piks)

#Calcul du poids : 

d=1/piks
total=as.integer(totaux)
g=calib(Y,
        d,
        method="logit",
        bounds=c(low=0.2,upp=4),
        total) 


min(g) #0.3340978
median(g) #0.8674619
max(g) #3.481854

plot(density(g))

## L'estimateur de _Horvitz-Thompson_ : 

colSums(Y*g/piks)

## L'estimateur calé de _Y_ : 

colSums(Y*g*d)

checkcalibration(Y, d, total=as.integer(totaux), g)

# Rassemblement des poids et des poids redressés : 

dk<-d
wk<-g/piks
hh<-merge$hh
XF<-data.frame(hh, dk, wk)
XF

# Traitement de l’enquête

#Reloading de _reponse_ : 

reponse<-read.sas7bdat("Reponse.sas7bdat")
reponse<-as.data.table(reponse)
reponse<-reponse[,":="(hh=CODEPAN)][,":="(CODEPAN=NULL)]

reponse

r<-reponse[,':='(varc=paste0("Q",sprintf("%02d",Q),sprintf("%02s",C)))][,reponse := ifelse(Q!=0,1,0)]
View(r)

# Merging _r_, _caract_ et _XF_ : 

r<-r[!duplicated(r[,c('hh','varc')]),]
merge<-r[caract,on="hh",nomatch=0,][XF,on="hh",nomatch=0]
merge

## Tri à plat

#Calcul de la somme des réponses de chaque choix sur chaque question sur l’ensemble de la population:


tri_a_plat <- merge[, .(
  varc=unique(varc),
  total= round(sum(wk*reponse),0)
), 
by=varc 
] 

tri_a_plat<-tri_a_plat[,2:3]

## Tri croisé : 

# Calcul de la somme des réponses de chaque choix sur chaque question sur l’ensemble de la populatio, en ventilant par les modalités de chaque variable :


### Variable : Catégorie socioprofessionnelle du chef du ménage _CSPCHEF_

CSPCHEF <- merge[, .(
  varc=unique(varc),
  CSPCHEF=unique(CSPCHEF),
  total_csp= sum(wk*reponse)), 
  by=.(varc,CSPCHEF) 
] 

CSPCHEF<-CSPCHEF[,c(1,2,5)]

index_cspchef<-tri_a_plat[CSPCHEF,on = c("varc"),nomatch=0][,':='(index=total_csp/total)]
index_cspchef

index_univers <- univers[, .(
  VAR=unique(VAR),
  LIBELLE,
  index_var= MAR/N
)
,
] 

index_univers

# Pondération de l'indice de chaque modalité


POND_index_cspchef<-index_cspchef[index_univers,on = c(CSPCHEF="LIBELLE"),nomatch=0][
  ,':='(indice_global=round((index/index_var)*100,digits=0))]

POND_index_cspchef

POND_index_cspchef<-POND_index_cspchef[order(POND_index_cspchef$varc, -POND_index_cspchef$indice_global),]

POND_index_cspchef

POND_index_cspchef<-POND_index_cspchef[,.(varc,CSPCHEF,indice_global)]
POND_index_cspchef<-POND_index_cspchef[!is.na(POND_index_cspchef$indice_global),]
POND_index_cspchef


### Variable : Age du chef du ménage _AGECHEF_

AGECHEF<- merge[, .(
  varc=unique(varc),
  AGECHEF=unique(AGECHEF),
  total_age= sum(wk*reponse)), 
  by=.(varc,AGECHEF) 
] 

AGECHEF<-AGECHEF[,c(1,2,5)]

index_age<-tri_a_plat[AGECHEF,on = c("varc"),nomatch=0][,':='(index=total_age/total)]
index_age

POND_index_age<-index_age[index_univers,on = c(AGECHEF="LIBELLE"),nomatch=0][
  ,':='(indice_global=round((index/index_var)*100,digits=0))]

POND_index_age<-POND_index_age[order(POND_index_age$varc, -POND_index_age$indice_global),]

POND_index_age<-POND_index_age[,.(varc,AGECHEF,indice_global)]
POND_index_age<-POND_index_age[!is.na(POND_index_age$indice_global),]
POND_index_age

### Variable : Région d'habitation _REGION_

REGION<- merge[, .(
  varc=unique(varc),
  REGION=unique(REGION),
  total_region= sum(wk*reponse)), 
  by=.(varc,REGION) 
] 


REGION<-REGION[,c(1,2,5)]

index_region<-tri_a_plat[REGION,on = c("varc"),nomatch=0][,':='(index=total_region/total)]
index_region

POND_index_region<-index_region[index_univers,on = c(REGION="LIBELLE"),nomatch=0][
  ,':='(indice_global=round((index/index_var)*100,digits=0))]

POND_index_region<-POND_index_region[order(POND_index_region$varc, -POND_index_region$indice_global),]

POND_index_region<-POND_index_region[,.(varc,REGION,indice_global)]
POND_index_region<-POND_index_region[!is.na(POND_index_region$indice_global),]


POND_index_region


### Variable : Nombre de personnes dans le ménage _NBPERS_

NBPERS<- merge[, .(
  varc=unique(varc),
  NBPERS=unique(NBPERS),
  total_nbpers= sum(wk*reponse)), 
  by=.(varc,NBPERS) 
] 


NBPERS<-NBPERS[,c(1,2,5)]

index_nbpers<-tri_a_plat[NBPERS,on = c("varc"),nomatch=0][,':='(index=total_nbpers/total)]
index_nbpers

POND_index_nbpers<-index_nbpers[index_univers,on = c(NBPERS="LIBELLE"),nomatch=0][
  ,':='(indice_global=round((index/index_var)*100,digits=0))]

POND_index_nbpers<-POND_index_nbpers[order(POND_index_nbpers$varc, -POND_index_nbpers$indice_global),]

POND_index_nbpers<-POND_index_nbpers[,.(varc,NBPERS,indice_global)]
POND_index_nbpers<-POND_index_nbpers[!is.na(POND_index_nbpers$indice_global),]


POND_index_nbpers


### Variable : Tranche Urbaine _TU_

TU<- merge[, .(
  varc=unique(varc),
  TU=unique(TU),
  total_tu= sum(wk*reponse)), 
  by=.(varc,TU) 
] 


TU<-TU[,c(1,2,5)]

index_tu<-tri_a_plat[TU,on = c("varc"),nomatch=0][,':='(index=total_tu/total)]
index_tu

POND_index_tu<-index_tu[index_univers,on = c(TU="LIBELLE"),nomatch=0][
  ,':='(indice_global=round((index/index_var)*100,digits=0))]

POND_index_tu<-POND_index_tu[order(POND_index_tu$varc, -POND_index_tu$indice_global),]

POND_index_tu<-POND_index_tu[,.(varc,TU,indice_global)]
POND_index_tu<-POND_index_tu[!is.na(POND_index_tu$indice_global),]


POND_index_tu

#Test de significativité
##H0 : Pas de d’association.
##H1 : Il existe une association.

#Variable CSPCHEF

tab_test_CSPCHEF=CrossTable(merge$varc,merge$CSPCHEF,chisq = TRUE)


#Variable AGECHEF

tab_test_AGECHEF=CrossTable(merge$varc,merge$AGECHEF,chisq = TRUE)
tab_test_AGECHEF


#Variable REGION

tab_test_REGION=CrossTable(merge$varc,merge$REGION ,chisq = TRUE)
tab_test_REGION


#Variable NBPERS

tab_test_NBPERS=CrossTable(merge$varc,merge$NBPERS,chisq = TRUE)
tab_test_NBPERS

#Variable TU

tab_test_TU=CrossTable(merge$varc,merge$TU,chisq = TRUE)
tab_test_TU
Y

### caractèrisation multivarié 

DM_C <- merge[,c(2:7)]

# Hyperparamètres
hp <- rpart.control(minsplit = 120, #L'effectif minimal pour séparer un noeud 
                    minbucket = 50, # L'effectif minimal dans chaque noeud terminal
                    maxdepth = 20, # Hauteur (profondeur) maximale de l'arbre
                    cp = 0 # Paramètre de pénalisation pour la complexité
)

# Construction de l'arbre
arbre_rpart <- rpart(Reponse~., # Y~X
                     DM_C, # Données
                     control = hp # Hyperparamètres
) 

# Afficher l'arbre
arbre_rpart

# Visualiser l'arbre
rpart.plot(arbre_rpart)

# Faire une prédiction avec l'arbre
hat_y <- predict(arbre_rpart, # Le modèle créé par rpart
                 DM_C, # Le jeu de données que l'on souhaite prédire
                 type = 'class' # Le type de prédiction que l'on souhaite obtenir (une classe ou une probabilité)
)

# Logistic regression

DM_L <- merge[,c(2:6,8)]
DM_L$Reponse_B<-as.factor(DM_L$Reponse_B)
model<-glm(Reponse_B~.,data=DM_L, family=binomial(link=logit))
summary(model)
ypredit<-model$fitted
o=order(ypredit)

# Créations des vecteurs correspondant aux 2 catégories :

hommes = c(50,70,110,60)

femmes = c(80,75,100,30)

# Création d'une matrice comparative :

tableau = matrix(c(hommes, femmes),2,4,byrow=T) # (2 : nombre de lignes et 4 nombres de colonnes (tranches salariales))

# Réalisation du test khi-deux - les résultats sont sauvegardés dans "khi_test"

khi_test = chisq.test(tableau)

khi_test # affiche le résultat du test

# Charger la librairie
library(party)
DM_P <- DM[,c(2:6,8)]
# Choisir les hyperparamètres
hp <- ctree_control(minsplit = 50, #L'effectif minimal pour séparer un noeud 
                    minbucket = 30, # L'effectif minimal dans chaque noeud terminal
                    maxdepth = 20, # Hauteur (profondeur) maximale de l'arbre
                    mincriterion = 0.3) # 1-p-valeur à partir de laquelle on souhaite cesser la croissance)

# Construire l'arbre
arbre_ctree <- ctree(Reponse_B~., DM_P,control = hp)

# Afficher l'arbre
arbre_ctree

# Visualiser l'arbre
plot(arbre_ctree)

# Faire une prédiction (donne une classe par défaut)

hat_y <- predict(arbre_ctree, # Le modèle
                 don) #