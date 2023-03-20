# ========================= 1st DATASET =========================
library(foreach)
library(meta)

dat <- read.csv('C:/Users/PapaBaGAYE/Desktop/KaggleV2-May-2016.csv',
                header =TRUE,sep=",")
dim(dat)
head(dat)
summary(dat)

table(dat$Alcoholism)
table(dat$Mo.show)

dat$No.show_bis <- ifelse(dat$No.show=="Yes",1, 0)
unique(dat$Neighbourhood)
mod <-glm( No.show_bis ~ Age + Gender + Scholarship+Hipertension+Diabetes+Alcoholism+Handcap+SMS_received,family=binomial(),data=dat)
summary(mod)
#stp <- stepAIC(mod)
#summary(stp)
mod_e <- glm(No.show_bis~SMS_received,family = binomial,data=dat)
summary(mod_e)
Neighbouhood <- c()
k <-1
id <- unique(dat$Neighbourhood)
res <- foreach(i=1:length(id), .combine='rbind')%do%{
  sub<-subset(dat,Neighbourhood==as.character(id[i]))
  if (dim(sub)[1]<10){
    NULL
  }
  else
  {
    mod_e <- glm(No.show_bis~SMS_received,family = binomial(),data=sub)
    Neighbouhood[k]<-as.character(id[i])
    k <- k+1
    summary(mod_e)[[12]][2,1:2]
                 
  }
}
dat_res<-as.data.frame(res)
dat_res <-cbind(dat_res,Neighbouhood)
colnames(dat_res) <- c("Estimate","StandardError","Neighbouhood")

m = metagen(Estimate, StandardError, studlab = Neighbouhood, data=dat_res)
forest(m)

# ========================= 2nd DATASET =========================
library(foreach)
library(meta)

dat <- read.csv('C:/Users/PapaBaGAYE/Desktop/clinvar_conflicting.csv', header = TRUE, sep = ",")
head(dat)
dim(dat)
summary(dat)
dat$CHROM <- as.numeric(dat$CHROM)
table(dat$CHROM)

for (i in c(3, 4, 8, 16, 20, 21, 28, 31, 32, 35, 39, 40, 41)){
  dat[,i] = as.factor(dat[,i])
}
summary(dat)

dat = subset(dat, select = c(CHROM, POS, REF, ALT, AF_ESP, AF_EXAC, AF_TGP, CLNSIGINCL, CLNVC,
                             MC, CLASS, Consequence, IMPACT, SYMBOL, Feature_type, Feature, SIFT,
                             PolyPhen, MOTIF_SCORE_CHANGE, LoFtool, CADD_PHRED, CADD_RAW, BLOSUM62)) 

summary(dat)

dat$Impact_bis <- as.factor(ifelse(dat$IMPACT %in% c("HIGH", "LOW"), 1, 0))

summary(dat$Impact_bis)

# Meta analyse par chromosome : étudier l'effet de CADD_PHRED sur Impact_bis

mod_e <- glm(Impact_bis ~ CADD_PHRED, family = binomial, data=dat)
summary(mod_e)

dat <- subset(dat, !is.na(CHROM))
table(dat$CHROM)
CHROM <- c()
k <-1
id <- unique(dat$CHROM)
res <- foreach(i=1:length(id), .combine='rbind')%do%{
  sub <- subset(dat, CHROM == as.character(id[i]))
  if (dim(sub)[1]<10){
    NULL
  }
  else
  {
    mod_e <- glm(Impact_bis ~ CADD_PHRED, family = binomial(), data = sub)
    CHROM[k] <- as.character(id[i])
    k <- k+1
    summary(mod_e)[[13]][2,1:2]
  }
}

dat_res <- as.data.frame(res)
dat_res <- cbind(dat_res, CHROM)
colnames(dat_res) <- c("Estimate", "StandardError", "CADD_PHRED")
dat_res

m = metagen(Estimate, StandardError, studlab = CHROM, data = dat_res)
forest(m)

# Meta analyse par chromosome : étudier l'effet de AF_EXAC sur Impact_bis

mod_e <- glm(Impact_bis ~ AF_EXAC, family = binomial, data=dat)
summary(mod_e)

dat <- subset(dat, !is.na(CHROM))
table(dat$CHROM)
CHROM <- c()
k <-1
id <- unique(dat$CHROM)
res <- foreach(i=1:length(id), .combine='rbind')%do%{
  sub <- subset(dat, CHROM == as.character(id[i]))
  if (dim(sub)[1]<10){
    NULL
  }
  else
  {
    mod_e <- glm(Impact_bis ~ AF_EXAC, family = binomial(), data = sub)
    CHROM[k] <- as.character(id[i])
    k <- k+1
    summary(mod_e)[[12]][2, 1:2]
  }
}

dat_res <- as.data.frame(res)
dat_res <- cbind(dat_res, CHROM)
colnames(dat_res) <- c("Estimate", "StandardError", "AF_EXAC")
dat_res

m = metagen(Estimate, StandardError, studlab = CHROM, data = dat_res)
forest(m)

# Meta analyse par chromosome : étudier l'effet de AF_EXAC sur Impact_bis

# Meta analyse par chromosome : étudier l'effet de AF_EXAC sur Impact_bis

#dat$PolyPhen_bis = as.factor(ifelse(dat$PolyPhen %in% c(' ', 'unknown'), 0, ifelse(dat$PolyPhen == 1,
#                                                                                  ifelse(dat$PolyPhen == "probably_damaging",
#                                                                                         2,
#                                                                                         3))))

#summary(dat$PolyPhen_bis)
#tabe(dat$PolyPhen, dat$PolyPhen_bis)

#model_vglm = vglm(PolyPhen_bis ~ CADD_PHRED, family = multinomial(), data = dat)
#summary(model_vglm)

