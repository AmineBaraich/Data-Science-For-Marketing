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
library(party)

caract <- as.data.table(read.sas7bdat("caract.sas7bdat"))

marge <- read.sas7bdat("marge_insee.sas7bdat")

reponse <- as.data.table(read.sas7bdat("Reponse.sas7bdat"))

reponse[, hh := CODEPAN][, CODEPAN := NULL]
reponse <- reponse[!duplicated(hh)]
merge <- caract[reponse, on="hh", nomatch=0]

A <- merge[, c(3, 4, 2, 5, 6), with=FALSE]

Y <- data.frame(sapply(A, function(x) data.table(model.matrix(~x-1, data=A))[,]))

names(Y)[1:4] <- c("<35", ">65", "35-50", "50-65")
names(Y)[5:9] <- c("centre EST", "Centre Ouest", "Nord", "R Parisienne", "Sud")
names(Y)[10:17] <- c("AGRI", "CADM", "CADS", "COMM", "EMPLOY", "INACTIF", "OUVRIER", "RETRAIT")
names(Y)[18:22] <- c("1", "2", "3", "4", "5 & PLUS")
names(Y)[23:27] <- c("< 20 000 individus", ">200 000 individus hors Agg. Parisienne", "20 000-200 000 individus", "Agg. Parisienne", "Rurale")

marge <- marge[-c(1, 2, 3, 4), ]

lib <- data.table(LIBELLE = names(Y))
univers <- inner_join(lib, marge, by="LIBELLE")

n <- nrow(merge)

N <- sum(univers$MAR[1:4])

totaux <- t(as.matrix(univers$MAR))

piks <- rep(n / N, n)

d <- 1 / piks
total <- as.integer(totaux)
g <- calib(Y, d, method="logit", bounds=c(low=0.2, upp=4), total)

min(g)
median(g)
max(g)

plot(density(g))

colSums(Y * g / piks)

colSums(Y * g * d)

checkcalibration(Y, d, total=as.integer(totaux), g)

dk <- d
wk <- g / piks
hh <- merge$hh
XF <- data.frame(hh, dk, wk)

reponse <- as.data.table(read.sas7bdat("Reponse.sas7bdat"))
reponse[, hh := CODEPAN][, CODEPAN := NULL]

r <- reponse[, varc := paste0("Q", sprintf("%02d", Q), sprintf("%02s", C))][, reponse := ifelse(Q != 0, 1, 0)]

r <- r[!duplicated(r[, c("hh", "varc"), with=FALSE])]
merge <- r[caract, on="hh", nomatch=0][XF, on="hh", nomatch=0]

tri_a_plat <- merge[, .(varc = unique(varc), total = round(sum(wk * reponse), 0)), by=varc][, 2:3, with=FALSE]

CSPCHEF <- merge[, .(varc = unique(varc), CSPCHEF = unique(CSPCHEF), total_csp = sum(wk * reponse)), by=.(varc, CSPCHEF)][, c(1, 2, 5), with=FALSE]

index_cspchef <- tri_a_plat[CSPCHEF, on="varc", nomatch=0][, index := total_csp / total]

index_univers <- data.table(univers)[, .(VAR = unique(VAR), LIBELLE, index_var = MAR / N)]

POND_index_cspchef <- index_cspchef[index_univers, on=c(CSPCHEF="LIBELLE"), nomatch=0][, indice_global := round((index / index_var) * 100, digits=0)]
POND_index_cspchef <- POND_index_cspchef[order(varc, -indice_global)]
POND_index_cspchef <- POND_index_cspchef[, .(varc, CSPCHEF, indice_global)]
POND_index_cspchef <- POND_index_cspchef[!is.na(indice_global)]

AGECHEF <- merge[, .(varc = unique(varc), AGECHEF = unique(AGECHEF), total_age = sum(wk * reponse)), by=.(varc, AGECHEF)][, c(1, 2, 5), with=FALSE]

index_age <- tri_a_plat[AGECHEF, on="varc", nomatch=0][, index := total_age / total]

POND_index_age <- index_age[index_univers, on=c(AGECHEF="LIBELLE"), nomatch=0][, indice_global := round((index / index_var) * 100, digits=0)]
POND_index_age <- POND_index_age[order(varc, -indice_global)]
POND_index_age <- POND_index_age[, .(varc, AGECHEF, indice_global)]
POND_index_age <- POND_index_age[!is.na(indice_global)]

REGION <- merge[, .(varc = unique(varc), REGION = unique(REGION), total_region = sum(wk * reponse)), by=.(varc, REGION)][, c(1, 2, 5), with=FALSE]

index_region <- tri_a_plat[REGION, on="varc", nomatch=0][, index := total_region / total]

POND_index_region <- index_region[index_univers, on=c(REGION="LIBELLE"), nomatch=0][, indice_global := round((index / index_var) * 100, digits=0)]
POND_index_region <- POND_index_region[order(varc, -indice_global)]
POND_index_region <- POND_index_region[, .(varc, REGION, indice_global)]
POND_index_region <- POND_index_region[!is.na(indice_global)]

NBPERS <- merge[, .(varc = unique(varc), NBPERS = unique(NBPERS), total_nbpers = sum(wk * reponse)), by=.(varc, NBPERS)][, c(1, 2, 5), with=FALSE]

index_nbpers <- tri_a_plat[NBPERS, on="varc", nomatch=0][, index := total_nbpers / total]

POND_index_nbpers <- index_nbpers[index_univers, on=c(NBPERS="LIBELLE"), nomatch=0][, indice_global := round((index / index_var) * 100, digits=0)]
POND_index_nbpers <- POND_index_nbpers[order(varc, -indice_global)]
POND_index_nbpers <- POND_index_nbpers[, .(varc, NBPERS, indice_global)]
POND_index_nbpers <- POND_index_nbpers[!is.na(indice_global)]

TU <- merge[, .(varc = unique(varc), TU = unique(TU), total_tu = sum(wk * reponse)), by=.(varc, TU)][, c(1, 2, 5), with=FALSE]

index_tu <- tri_a_plat[TU, on="varc", nomatch=0][, index := total_tu / total]

POND_index_tu <- index_tu[index_univers, on=c(TU="LIBELLE"), nomatch=0][, indice_global := round((index / index_var) * 100, digits=0)]
POND_index_tu <- POND_index_tu[order(varc, -indice_global)]
POND_index_tu <- POND_index_tu[, .(varc, TU, indice_global)]
POND_index_tu <- POND_index_tu[!is.na(indice_global)]

tab_test_CSPCHEF <- CrossTable(merge$varc, merge$CSPCHEF, chisq=TRUE)
tab_test_AGECHEF <- CrossTable(merge$varc, merge$AGECHEF, chisq=TRUE)
tab_test_REGION <- CrossTable(merge$varc, merge$REGION, chisq=TRUE)
tab_test_NBPERS <- CrossTable(merge$varc, merge$NBPERS, chisq=TRUE)
tab_test_TU <- CrossTable(merge$varc, merge$TU, chisq=TRUE)

DM_C <- merge[, 2:7, with=FALSE]

hp <- rpart.control(minsplit=120, minbucket=50, maxdepth=20, cp=0)
arbre_rpart <- rpart(Reponse ~ ., DM_C, control=hp)

rpart.plot(arbre_rpart)

hat_y <- predict(arbre_rpart, DM_C, type="class")

DM_L <- merge[, c(2:6, 8), with=FALSE]
DM_L$Reponse_B <- as.factor(DM_L$Reponse_B)
model <- glm(Reponse_B ~ ., data=DM_L, family=binomial(link="logit"))

ypredit <- model$fitted.values
o <- order(ypredit)

hommes <- c(50, 70, 110, 60)
femmes <- c(80, 75, 100, 30)
tableau <- matrix(c(hommes, femmes), 2, 4, byrow=TRUE)
khi_test <- chisq.test(tableau)

DM_P <- DM[, c(2:6, 8)]
hp_ctree <- ctree_control(minsplit=50, minbucket=30, maxdepth=20, mincriterion=0.3)
arbre_ctree <- ctree(Reponse_B ~ ., DM_P, control=hp_ctree)

plot(arbre_ctree)

hat_y_ctree <- predict(arbre_ctree, don)