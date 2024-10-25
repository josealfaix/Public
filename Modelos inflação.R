library(BETS)
library(xts)
library(tstools)
library(tidyverse)
library(forecast)
library(lmtest)

acum_p <- function(data,n){
  factor <- (1+(data/100))
  prod <- RcppRoll::roll_prodr(factor, n=n)
  final <- (prod-1)*100
  
  return(final)
}

inicio <- Sys.Date()-3660

ipcadata <- BETSget(433, from = paste(inicio), to = paste(Sys.Date()), frequency = 12, data.frame = TRUE)
ipcadata <- data.frame(ipcadata)
ipcaxts <- xts(ipcadata$value, order.by = ipcadata$date)
colnames(ipcaxts) <- "ipca"

ipcaxts$ipcatrim <- acum_p(ipcaxts$ipca,3)
ipcaxts$ipcaanual <- acum_p(ipcaxts$ipca,12)

ipcaxts$ipcaanual4 <- lag.xts(ipcaxts$ipcaanual, -4)

ipcaxts$y <- ipcaxts$ipcaanual4 - ipcaxts$ipcaanual
ipcaxts$x1 <-  ipcaxts$ipcaanual - lag.xts(ipcaxts$ipcaanual, 1)
ipcaxts$x2 <- lag.xts(ipcaxts$x1,1)
ipcaxts$x3 <- lag.xts(ipcaxts$x2,1)
ipcaxts$x4 <- lag.xts(ipcaxts$x3,1)
ipcaxts <- ipcaxts%>%na.omit()

plot(ipcaxts[,1:5], legend.loc = "topright", main = "Medidas de inflação \n e variável de interesse, 'y'.")


################################

# Modelo 1 -> Recursivo
#Prevê a inflação um ano à frente usando lags do próprio IPCA. (Stock and Watson, 1999)

# π4t+4−πt=α+β(L)(πt−πt−1)+ϵt

modelo1 = lm(y ~x1+x2+x3+x4, data = ipcaxts)
summary(modelo1)

newdata <- xts(ipcadata$value, order.by = ipcadata$date)
colnames(newdata) <- "ipca"

newdata$ipcaanual <- acum_p(newdata$ipca,12)
newdata$x1 <- newdata$ipcaanual - lag.xts(newdata$ipcaanual,1)

predata <- tail(newdata$x1,6)

pred.modelo1 <- predict.lm(lm(y ~x1, data = ipcaxts), newdata = predata, se.fit=TRUE)

pred.modelo1

accuracy(modelo1)

coeftest(modelo1)

#Previsão usando o pacote forecast

arima_ipca <- auto.arima(ipcaxts$y, xreg = ipcaxts[,6:9])
summary(arima_ipca)

checkresiduals(arima_ipca)


################################

# Modelo 2 --> Ingênuo
# A previsão da inflação um ano à frente é simplesmente a taxa de crescimento anualizada dos mesmos quatro trimestres no ano anterior

#π4t+4=π4t=(∏i=14(1+πt−i))1/4−1

modelo2 <- tibble(Data = time(tail(ipcaxts$ipcaanual))+120, 'Previsão Ingênua' = tail(ipcaxts$ipcaanual))

print(modelo2)

RMSE_ingenuo <- sum((ipcaxts$ipcaanual4 - ipcaxts$ipcaanual)^2)^(1/2)

print(RMSE_ingenuo)

###############################

#Modelo 3 --> Curva de Phillips 

#Consiste no Modelo 1 adicionando medidas de atividade econômica e lags de inflação escolhendo o " lag ótimo" por meio do Critério Bayesiano de Informação (BIC). a equação de regressão é:
  
  #π4t+4−πt=α+β(L)(πt−πt−1)+γ(L)xt+ϵt

#Em que , xt
#reprsentam as medidas de atividade econômica. Neste caso, usarei o hiato PIB como variável de atividade econômica. Como proxy para o PIB mensal usarei o indice IBC-Br do Banco Central. Para calcular o hiato usarei o resíduo de um modelo ARIMA da série.

ibc <- BETSget(24364, from = paste(inicio), to=paste(Sys.Date()), frequency = 12, data.frame =TRUE)

ibcxts <- xts(ibc[,2], order.by = ibc[,1], frequency = 12)
colnames(ibcxts) <- "ibc"

autoplot(ibcxts, color = 'red')

ibc_arima <- auto.arima(ibcxts, seasonal = TRUE)

summary(ibc_arima)

checkresiduals(ibc_arima)

#Hiato do IBC

ibc_hiato <- xts(ibc_arima[["residuals"]], order.by = ibc$date)
autoplot(ibc_hiato)

ibc_hiato <- cbind(ibc_hiato, lag.xts(ibc_hiato, c(1:4)))
colnames(ibc_hiato) <- c("ibc_hiato",paste("ibc_hiato lag", c(1:4)))

ipcaxts <- merge.xts(ipcaxts,ibc_hiato)%>%na.omit()

ipca_phillips <- auto.arima(ipcaxts$y, xreg = ipcaxts[,6:ncol(ipcaxts)])

summary(ipca_phillips)
coeftest(ipca_phillips)

ipca_phillips <- Arima(ipcaxts$y, order = c(2,0,0), xreg = ipcaxts[,6:ncol(ipcaxts)])
summary(ipca_phillips)

##############################

#Modelo 4 - Núcleos de inflação

#Este modelo usa as diferentes medidas de inflação núcleo como regressores, como indicado em Pasaogullari and Meyer (2010) e fica:
  
  #π4t+4=α+λ(L)π∗t+ϵt

#Em que π∗t
#é o tipo de medida de inflação núcleo. Foram usados no máximo 4 defasagens e a quantidade ótima de defasagens foi determinada pelo critério Bayesiano de informação.

## Pegar núcleos

codes = c(4466, 11426, 11427,16121, 16122, 27838, 27839)

nucleos = BETSget(codes, from='2006-07-01', data.frame=T)

data_nucleos = matrix(NA, nrow=nrow(nucleos[[1]]),
                      ncol=length(codes))
for(i in 1:length(codes)){
  data_nucleos[,i] = t(nucleos[[i]]$value)
}

colnames(data_nucleos) = c('ipca_ms', 'ipca_ma', 'ipca_ex0',
                           'ipca_ex1', 'ipca_dp', 'ipca_ex2',
                           'ipca_ex3')

nucleos_vm =
  data_nucleos %>%
  as_tibble() %>%
  mutate(date = nucleos[[1]]$date) %>%
  select(date, everything())

nucleos_12m = 
  data_nucleos %>%
  ts(start=c(2006,07), freq=12) %>%
  acum_p(12) %>%
  as_tibble() %>%
  mutate(date = nucleos[[1]]$date) %>%
  select(date, everything()) %>%
  drop_na()


nucleoxts <- xts(nucleos_12m[,-1], order.by = nucleos_12m$date)

nucleodata <- merge.xts(ipcaxts$ipcaanual4, nucleoxts, lag.xts(nucleoxts, c(1:4)))%>%na.omit()
modelo4 <- auto.arima(nucleodata[,1], xreg = nucleodata[,-1])

summary(modelo4)

plot(nucleoxts)

coeftest(modelo4)

rbind(c(accuracy(modelo1), "-"), accuracy(ipca_phillips), accuracy(modelo4))








