library(readr)
library(tidyverse)
library(matrixStats)
####### Udregn vinklen til hver tid ######

benpersondata <- function(person, individ, benet){
  vb <- read.table(person, header = TRUE)
  ###Udtag de punkter der skal bruges for at regne vinklen i knaeet og tiden.
  vb <- vb[,c(1,5,7,8,10,14,16)]

  ###omdoeb vektorerne for at goere det nemmere
  colnames(vb) <- c("tid", "ankelv", "ankelAP", "knev", "knevAP", "hoftev", "hofteAP")
  
  ## AP er er x koordinatet, og v er y koordinatet
  
    
  # Indsaet hvilket ben det er:
  vb <- transform(vb, ben = benet)
  
  # Omrok?r s?jler, s? det bliver nemmere at arbejde med
  vb <- vb[,c(1,8,2,3,4,5,6,7)]
  
  ### Udregn laengden af underbenet
  vb <- transform(vb, ankelknaeL = sqrt((`ankelAP` - `knevAP`)^2 + (`ankelv` - `knev`)^2))
  ### Udregn laengden af laaet
  vb <- transform(vb, knaehofteL = sqrt((knevAP - hofteAP)^2 + (knev - hoftev)^2))
  ### Udregn laengden fra anklen til hoften.
  vb <- transform(vb, ankelhofteL = sqrt((ankelAP - hofteAP)^2 + (ankelv - hoftev)^2))
  
  
  ### Udregn haeldningen paa den linje der foelger laaret og den linje der foelger underbenet
  ## Bruger formlen (y1-y2) / (x1-x2) for at finde h?ldningen p? den rette linje
  
  vbO <- (vb$hoftev - vb$knev) / (vb$hofteAP - vb$knevAP) ## Overbenet
  vbU <- (vb$knev - vb$ankelv) /(vb$knevAP - vb$ankelAP) ## Underbenet

  
  ## S?rger her for at vi ikke dividere med 0
  
  for (i in 1:7000){
    if(vbO[i] == Inf){
      vbO[i] <- 0
    }
  }
  
  for (i in 1:7000){
    if(vbU[i] == Inf){
      vbU[i] <- 0
    }
  }
  
  ###udregn sk?ring med y-aksen linjen der f?lger overbenet

  ## F?lgende formel bruges: b = y - a * x
  
  b1 <- (vb$hoftev - vbO * vb$hofteAP)
  
  ## Udregn sk?ring med y-aksen for linjen der f?lger underbenet
  
  b2 <- (vb$ankelv - vbU * vb$ankelAP)

  
  ### Lav en tom matrix til at inds?tte de x-v?rdier hvor linjerne sk?r -4000 (-4000 er tilf?ldigt) p? y-aksen.
  XO <- NULL
  XO <- (-(b1-(-4000)))/vbO
  XU <- NULL
  XU <- (-(b2-(-4000)))/vbU

  
  ### Tjek for haeldning der er 0, hvis den er 0, s? ins?ttes x-koordinatet for hoften i stedet:
  for (i in 1:7000){
    if(XU[i] == -Inf){
      XU[i] <- vb$hofteAP[i]
    }
  }
  
  for (i in 1:7000){
    if(XO[i] == -Inf){
      XO[i] <- vb$hofteAP[i]
    }
  }
  ###lav en tom vektor i datasaettet, hvor vinklerne fra knaeet skal fyldes ind.
  vb <- transform(vb, knevinkel = NA)
  
  ### Lav et for loop, hvor vinklen for knaeet bliver udregnet, b?de hvor vinklen er under 180 grader og hvor vinklen for knaeet er over 180 grader
  ## f = (Pi-arccos((a^2-b^2-c^2)/(2*b*c)))*(180/Pi), Hvis der overeksponeres (XO[i]<XU[i], saa tages 360 - f)
  for (i in 1:length(XO)){
    if(XO[i] > XU[i]){
      vb$knevinkel[i] <- ((pi - acos((vb$ankelhofteL[i]^2 - vb$knaehofteL[i]^2 - vb$ankelknaeL[i]^2)/(2 *vb$knaehofteL[i] * vb$ankelknaeL[i])))*180/pi)}
    else {
      vb$knevinkel[i] <- 360 - ((pi - acos((vb$ankelhofteL[i]^2 - vb$knaehofteL[i]^2 - vb$ankelknaeL[i]^2)/(2 *vb$knaehofteL[i] * vb$ankelknaeL[i])))*180/pi)}}
  
  
  #Find alle lokale minima, som der findes hvor vinklen p? knaeet er mindre end 140
  findpeaks <- function(vec,bw = 1, x.coo = c(1 : length(vec)))
  {
    pos.x.min <- NULL
    pos.y.min <- NULL 	
    
    for(i in 1:(length(vec)-1)) 	{ 		
      
      if((i+1+bw)>length(vec))
      {sup.stop <- length(vec)}
      
      else
      {sup.stop <- i+1+bw}
      
      if((i-bw)<1)
      {inf.stop <- 1}
      
      else{inf.stop <- i-bw}
      
      subset.sup <- vec[(i+1):sup.stop]
      subset.inf <- vec[inf.stop:(i-1)]
      
      is.max   <- sum(subset.inf > vec[i]) == 0
      is.nomin <- sum(subset.sup > vec[i]) == 0
      no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
      no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
      
      if(no.max & no.nomin){
        pos.x.min <- c(pos.x.min,x.coo[i])
        pos.y.min <- c(pos.y.min,vec[i])
      }
    }
    U <- NULL
    for(i in 1:length(pos.x.min)){
      if(pos.y.min[i] < 140){
        U[i] <- pos.x.min[i]
      }
      U <- U[!is.na(U)]
    }
    return(U)
  }
  ### Lav en vector og inds?t de tider, hvor der er lokalt minima:
  P <- NULL
  P <- findpeaks((vb$knevinkel))
  
  ### Reduc?r data til 10 skridt. De 10 skridt der benyttes er de af 10 skridt der kommer f?r midten af vektoren, s? fors?gspersonen har v?ret i gang i noget tid
  vbny <- subset(vb, tid >= P[round(length(P)/2, digits = 0)-10] & tid <= P[round(length(P)/2, digits = 0)])
  
  ### Udpluk de 10 tider (punkter), hvor det lokale minima sker
  Y <- NULL
  Y <- P[(round(length(P)/2, digits = 0) -10) : (round(length(P)/2, digits = 0))]
  
  ###Skub minima tidsaksen, saa den starter i 1
  Y <- Y-(min(vbny$tid)-1)
  
  ###skub tidsaksen i datasaettet, saa den starter i 1.
  vbny$tid <- vbny$tid - (min(vbny$tid)-1)
  
  ###Hiv i datasaettet saa hvert skridt er lige langt i forhold til tiden
  vbny$tid[Y[1] :(Y[2]-1)] <- seq(1:(Y[2] - Y[1]))/(Y[2] - Y[1])
  vbny$tid[Y[2] :(Y[3]-1)] <- 1 + seq(1:(Y[3] - Y[2]))/(Y[3] - Y[2])
  vbny$tid[Y[3] :(Y[4]-1)] <- 2 + seq(1:(Y[4] - Y[3]))/(Y[4] - Y[3])
  vbny$tid[Y[4] :(Y[5]-1)] <- 3 + seq(1:(Y[5] - Y[4]))/(Y[5] - Y[4])
  vbny$tid[Y[5] :(Y[6]-1)] <- 4 + seq(1:(Y[6] - Y[5]))/(Y[6] - Y[5])
  vbny$tid[Y[6] :(Y[7]-1)] <- 5 + seq(1:(Y[7] - Y[6]))/(Y[7] - Y[6])
  vbny$tid[Y[7] :(Y[8]-1)] <- 6 + seq(1:(Y[8] - Y[7]))/(Y[8] - Y[7])
  vbny$tid[Y[8] :(Y[9]-1)] <- 7 + seq(1:(Y[9] - Y[8]))/(Y[9] - Y[8])
  vbny$tid[Y[9] :(Y[10]-1)] <- 8 + seq(1:(Y[10] - Y[9]))/(Y[10] - Y[9])
  vbny$tid[Y[10] :Y[11]] <- 9 + seq(1:((Y[11] - Y[10])+1))/((Y[11] - Y[10])+1)
  
  ###Indsaet nummeret p? personen
  vbny <- transform(vbny, id = individ)
  ###ryk rundt p? faktorerne som man gerne vil have det st?r.
  id_vinkeldata <- vbny[,c(13,2,1,12)]
  rownames(id_vinkeldata) <- 1:nrow(id_vinkeldata)
  
  gennemsnit <- function(x, id, b){
    B <- sapply(x$tid, function(t){c(sin(2 * pi * seq(1 : 10) * t), cos(2 * pi * seq(1 : 10) * t))})
    B <- t(B)
    model2alt <- lm(knevinkel ~ B, data = x)
    B_new <- cbind(1, t(sapply(seq(0, 1, length = 101), function(t){c(sin(2*pi*seq(1:10)*t), cos(2*pi*seq(1:10)*t))})))
    y_new <- B_new %*% matrix(ncol = 1, data  = coef(model2alt))
    y_new <- t(y_new)
    colnames(y_new) <- 1:ncol(y_new)
    y_new <- cbind.data.frame(id,b, y_new)
    return(y_new)
  }
  Y <- NULL
  Y <- gennemsnit(id_vinkeldata, individ, benet)
  return(Y)
}


## Indl?s al data, og udregn vinklerne ##

data1h <-benpersondata("id1h.tsv", "1", "dominerende")
data1v <-benpersondata("id1v.tsv", "1", "ikke dominerende")

data2h <-benpersondata("id2h.tsv", "2", "ikke dominerende")
data2v <-benpersondata("id2v.tsv", "2", "dominerende")

data3h <-benpersondata("id3h.tsv", "3", "dominerende")
data3v <-benpersondata("id3v.tsv", "3", "ikke dominerende")

data4h <-benpersondata("id4h.tsv", "4", "dominerende")
data4v <-benpersondata("id4v.tsv", "4", "ikke dominerende")

data5h <-benpersondata("id5h.tsv", "5", "ikke dominerende")
data5v <-benpersondata("id5v.tsv", "5", "dominerende")

data6h <-benpersondata("id6h.tsv", "6", "dominerende")
data6v <-benpersondata("id6v.tsv", "6", "ikke dominerende")

data7h <-benpersondata("id7h.tsv", "7", "dominerende")
data7v <-benpersondata("id7v.tsv", "7", "ikke dominerende")

data8h <-benpersondata("id8h.tsv", "8", "dominerende")
data8v <-benpersondata("id8v.tsv", "8", "ikke dominerende")

data9h <-benpersondata("id9h.tsv", "9", "dominerende")
data9v <-benpersondata("id9v.tsv", "9", "ikke dominerende")

data10h <-benpersondata("id10h.tsv", "10", "dominerende")
data10v <-benpersondata("id10v.tsv", "10", "ikke dominerende")

data11h <-benpersondata("id11h.tsv", "11", "dominerende")
data11v <-benpersondata("id11v.tsv", "11", "ikke dominerende")

data12h <-benpersondata("id12h.tsv", "12", "dominerende")
data12v <-benpersondata("id12v.tsv", "12", "ikke dominerende")

data13h <-benpersondata("id13h.tsv", "13", "dominerende")
data13v <-benpersondata("id13v.tsv", "13", "ikke dominerende")

data14h <-benpersondata("id14h.tsv", "14", "dominerende")
data14v <-benpersondata("id14v.tsv", "14", "ikke dominerende")

data15h <-benpersondata("id15h.tsv", "15", "dominerende")
data15v <-benpersondata("id15v.tsv", "15", "ikke dominerende")

data16h <-benpersondata("id16h.tsv", "16", "dominerende")
data16v <-benpersondata("id16v.tsv", "16", "ikke dominerende")


#bind alle personer sammen
A <- rbind(data1h, data1v, data2h, data2v, data3h, data3v, data4h, data4v, data5h, data5v,data6h, data6v, data7h, data7v, data8h, data8v, data9h, data9v, data10h, data10v, data11h, data11v, data12h, data12v, data13h, data13v,data14h, data14v, data15h, data15v, data16h, data16v)
gskridt <- as.data.frame(A)
#roter dataen, saa det kan plottes i ggplot
gskridt1 <- gather(gskridt, key = "time", value = "vinkel", 3:103)
gskridt2 <- mutate(gskridt1, time = parse_double(time))
gskridt3 <- arrange(gskridt2, id, b, time)

#### Plot alle individers vinkel under et enkelt skridt
ggplot(filter(gskridt3,id %in% c("1", "2", "8", "16") ) , aes(time, vinkel, col = b)) + geom_line(aes(group = id)) + facet_grid(b~id)
View(gskridt3)
gskridt3



#### PERMUTATIONSTEST FOR OM DER ER FORSKEL P? DOMINERENDE OG IKKE DOMINERENDE BEN ####


### Udtr?k v?rdier for det dominante og det ikke dominante ben
domi <- rep(NA, 101)
i_domi <- rep(NA, 101)

for (i in 1:101){
domi[i] <- tapply(t(gskridt[i+2]),gskridt$b, mean)[1]
i_domi[i] <- tapply(t(gskridt[i+2]),gskridt$b, mean)[2]
}

##### Laver et plot der viser gennemsnittet af alle personer p? deres henholdsvis dominerende og ikke dominerende ben

domi_snit <- cbind.data.frame(1:length(domi), domi, "Dominerende")
colnames(domi_snit) <- c("tid", "gennemsnits vinkel af alle personer", "ben")
i_domi_snit <- cbind.data.frame(1:length(i_domi), i_domi, "Ikke dominerende")
colnames(i_domi_snit) <- c("tid", "gennemsnits vinkel af alle personer", "ben")
begge_snit <- rbind(domi_snit,i_domi_snit)
ggplot(begge_snit, aes(tid, `gennemsnits vinkel af alle personer`, col = ben)) + geom_line()



### Den observerede forskel p? kurverne mellem dominant og ikke dominant ben
f <- function(x,y) {(x-y)^2}
diff_obs <- sqrt(sum(f(domi, i_domi))/101)
diff_obs


### F?rst udtr?kkes hvert individ:

id1 <- subset(gskridt, gskridt$id == 1)
id2 <- subset(gskridt, gskridt$id == 2)
id3 <- subset(gskridt, gskridt$id == 3)
id4 <- subset(gskridt, gskridt$id == 4)
id5 <- subset(gskridt, gskridt$id == 5)
id6 <- subset(gskridt, gskridt$id == 6)
id7 <- subset(gskridt, gskridt$id == 7)
id8 <- subset(gskridt, gskridt$id == 8)
id9 <- subset(gskridt, gskridt$id == 9)
id10 <- subset(gskridt, gskridt$id == 10)
id11 <- subset(gskridt, gskridt$id == 11)
id12 <- subset(gskridt, gskridt$id == 12)
id13 <- subset(gskridt, gskridt$id == 13)
id14 <- subset(gskridt, gskridt$id == 14)
id15 <- subset(gskridt, gskridt$id == 15)
id16 <- subset(gskridt, gskridt$id == 16)




### Lav en vektor hvor de to ben ender (Et af hvert ben for hvert individ i hver kop)

d0 <- NULL
d1 <- NULL
d2 <- NULL
ML <- NULL
ct <- function(g){
  flips <- sample(c(0,1), 
                  size = 1, 
                  replace = TRUE, 
                  prob = c(0.5,0.5))
  
  {if (flips == 0)
    {d0 <- sample_n(g, 2)}
    else{d0 <- g} 
    return(d0)}}


for (i in 1:5000){
#Lav et coinflip, for at bestemme hvilken vektor den skal ind i
flips <- sample(c(0,1), 
                size = 1, 
                replace = TRUE, 
                prob = c(0.5,0.5))

{if (flips == 0){
  d0 <- sample_n(id1, 2)
}
else{
  d0 <- id1
}}

d1 <- d0[1,]
d2 <- d0[2,]

k1 <- ct(id2)
d1 <- rbind(d1, k1[1,])
d2 <- rbind(d2, k1[2,])

k2 <- ct(id3)
d1 <- rbind(d1, k2[1,])
d2 <- rbind(d2, k2[2,])

k3 <- ct(id4)
d1 <- rbind(d1, k3[1,])
d2 <- rbind(d2, k3[2,])

k4 <- ct(id4)
d1 <- rbind(d1, k4[1,])
d2 <- rbind(d2, k4[2,])

k5 <- ct(id5)
d1 <- rbind(d1, k5[1,])
d2 <- rbind(d2, k5[2,])

k6 <- ct(id6)
d1 <- rbind(d1, k6[1,])
d2 <- rbind(d2, k6[2,])

k7 <- ct(id7)
d1 <- rbind(d1, k7[1,])
d2 <- rbind(d2, k7[2,])

k8 <- ct(id8)
d1 <- rbind(d1, k8[1,])
d2 <- rbind(d2, k8[2,])

k9 <- ct(id9)
d1 <- rbind(d1, k9[1,])
d2 <- rbind(d2, k9[2,])

k10 <- ct(id10)
d1 <- rbind(d1, k10[1,])
d2 <- rbind(d2, k10[2,])

k11 <- ct(id11)
d1 <- rbind(d1, k11[1,])
d2 <- rbind(d2, k11[2,])

k12 <- ct(id12)
d1 <- rbind(d1, k12[1,])
d2 <- rbind(d2, k12[2,])

k13 <- ct(id13)
d1 <- rbind(d1, k13[1,])
d2 <- rbind(d2, k13[2,])

k14 <- ct(id14)
d1 <- rbind(d1, k14[1,])
d2 <- rbind(d2, k14[2,])

k15 <- ct(id15)
d1 <- rbind(d1, k15[1,])
d2 <- rbind(d2, k15[2,])

k16 <- ct(id16)
d1 <- rbind(d1, k16[1,])
d2 <- rbind(d2, k16[2,])

mean_d1 <- colMeans(d1[,3:103])
mean_d2 <- colMeans(d2[,3:103])

ML[i] <- sqrt(sum(f(mean_d1, mean_d2))/101)
print(i)
}



### Tegner et histogram over resultaterne af simulationen, og indsaetter den observerede difference ###


hist(ML, main = "Værdierne udregnet ved  permutationstest", col = "blue", border = TRUE, xlab = "Værdier af teststørrelsen", ylab = "frekvens")
abline(v=diff_obs, col = "red")

### Udregner p-vaerdien
K <- (ML > diff_obs)
A_KV <- sum(K, na.rm=TRUE)
K
A_KV
pvalue <- 1-A_KV/length(K)


pvalue
### Da p-vaerdien er over 0.05, accepterer vi 0 hypotesen, altsaa er der ingen forskel paa vinklen paa dominant og ikke dominant ben

#### Find forskel paa kurverne og proev at lav et 95 % praediktionsinterval paa forskellen (Ogsaa kaldet limit of agreement) ####


f1 <- data1h[3:103] - data1v[3:103]
f2 <- data2h[3:103] - data2v[3:103]
f3 <- data3h[3:103] - data3v[3:103]
f4 <- data4h[3:103] - data4v[3:103]
f5 <- data5h[3:103] - data5v[3:103]
f6 <- data6h[3:103] - data6v[3:103]
f7 <- data7h[3:103] - data7v[3:103]
f8 <- data8h[3:103] - data8v[3:103]
f9 <- data9h[3:103] - data9v[3:103]
f10 <- data10h[3:103] - data10v[3:103]
f11 <- data11h[3:103] - data11v[3:103]
f12 <- data12h[3:103] - data12v[3:103]
f13 <- data13h[3:103] - data13v[3:103]
f14 <- data14h[3:103] - data14v[3:103]
f15 <- data15h[3:103] - data15v[3:103]
f16 <- data16h[3:103] - data16v[3:103]

af <- rbind(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16)


##### Limits of Agreement ######
LoAm <- NULL
LoAp <- NULL
for (i in 1:101){
  LoAp[i] <- gaf[i] + 1.96 * sdaf[i]
  LoAm[i] <- gaf[i] - 1.96 * sdaf[i]
}


LoA <- rbind (LoAp, LoAm)
matplot(t(af), type = "l", col = "black", xlim = c(0,101), ylim = c(-25,12), xlab = "Tid", ylab = "Vinkelforskel");par(new=TRUE);matplot(LoAm, type = "l", col = "red", xlim = c(0,101), ylim = c(-25, 12), ylab = " ");par(new=TRUE);matplot(LoAp, type = "l", col = "red", xlim = c(0,101), ylim = c(-25, 12), ylab = " ")



##### Konfidensinterval af gennemsnittet af forskellen #####

gaf <- colMeans(af)
gaf1 <- transform(gaf, tid = c(1:101))
colnames(gaf1) <- c("Gennemsnitsforskel", "tid")

## Konfidensinterval for observeret data ##

sdaf <- colSds(as.matrix(af))
sdaf1 <- transform(sdaf, tid = c(1:101))
sdaf1 <- sdaf1[,c(2,1)]
colnames(sdaf1) <- c("tid", "sd")

KIp <- NULL
KIm <- NULL
for (i in 1:101){
  KIp[i] <- gaf[i] + 1.96 * sdaf[i] * 1/sqrt(16)
  KIm[i] <- gaf[i] - 1.96 * sdaf[i] * 1/sqrt(16)
}


KI <- NULL
KI <- c(KIp, KIm, gaf)
KI <- transform(KI, tid = 1:101)

Kurve <- c(rep("Øvre", each = 101), rep("Nedre", each = 101), rep("Gennemsnitlige forskel", each = 101))
KI <- transform(KI, Kurve = Kurve)
colnames(KI) <- c("Gennemsnit af forskel", "tid", "Kurve")


##### Plotter gennemsnittet af forskellen af de 16 personer til tiden t




## Punktvis baand ##

PB <- NULL
PB1 <- NULL
for (i in 1:5000){
  sampl1 <- sample(16, replace = T)
  g1 <- af[sampl1, ]
  gag <- colMeans(g1)
  for (j in 1:101){
    PB[j] <- KIm[j] < gag[j] && gag[j] < KIp[j]
  }
  PB1[i] <- (sum(PB, na.rm=TRUE) / 101)
}
mean(PB1)
########

ggplot(KI, aes(tid, `Gennemsnit af forskel`, col = Kurve)) + geom_line() 








##Joint baand##
### Traekker 16 tilfaeldige datasaet ud, med tilbagelÃ¦gning:

S <- NULL
TF <- NULL
for (i in 1:5000){
sampl1 <- sample(16, replace = T)
g1 <- af[sampl1, ]
gag <- colMeans(g1)


  for (j in 1:101){
  TF[j] <- KIm[j] < gag[j] && gag[j] < KIp[j]
  TF1 <- sum(TF, na.rm=TRUE)
  }
  if (TF1 == 101) {
  S[i] <- 1
  } else {
  S[i] <- 0
  }
}
sum(S)/5000


## Det ses at det kun er omkring 80% der ligger inden for intervallet. Vinder ud af hvad konstanten skal vÃ¦re for at det bliver et 95% konfidensinterval
## JOINT B?ND
KIp1 <- NULL
KIm1 <- NULL
for (i in 1:101){
  KIp1[i] <- gaf[i] + 2.8 * sdaf[i] * 1/sqrt(16)
  KIm1[i] <- gaf[i] - 2.8 * sdaf[i] * 1/sqrt(16)
}

P <- NULL
X <- NULL
for (i in 1:5000){
  a <- sample(16, replace = T)
  b <- af[a,]
  c <- colMeans(b)
  for (j in 1:101){
    X[j] <- KIm1[j] < c[j] && c[j] < KIp1[j]
    d <- sum(X, na.rm = TRUE)
  }
  if (d == 101){
    P[i] <- 1
  } else {
    P[i] <- 0
  }
}

sum(P)/5000
1.96 = 80
2.1 = 85
2.3 = 90
2.5 = 94
2.57 = 95

c <- c(1.96, 2.1, 2.3, 2.5, 2.57, 2.8)
p <- c(80, 85, 90, 94, 95, 98)
data <- data.frame(c, p)

ggplot(data, aes(c, p)) + geom_line() + geom_point() + labs(x = "Værdier for c", y = "Procent der liger helt indenfor intervallet")
plot(c, p, type = "l")
### Laver et plot med gennemsnittet af forskellen for vores personer med bÃ¥de det standard 95 % KI og LoA KI


KI1 <- c(KIp1, KIm1, gaf)

FinalKI <- cbind(KI, KI1)

ggplot(FinalKI, aes(tid,`Gennemsnit af forskel`, col = Kurve)) + geom_line() + geom_line(aes(tid, KI1))
View(FinalKI)

