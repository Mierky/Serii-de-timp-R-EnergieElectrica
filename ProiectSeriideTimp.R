
# Vector cu toate pachetele (fara duplicate)
packages <- unique(c(
  "tsibble", "readxl", "fpp3", "fpp2", "vars", "tseries", "urca", "stats",
  "changepoint", "dplyr", "uroot", "TSA", "FinTS", "gt", "tidyverse",
  "fGarch", "forecast", "fDMA", "lmtest", "lubridate", "nortsTest", 
  "tsDyn", "dynlm", "aTSA"
))

# Instaleaza doar pachetele care nu sunt deja instalate
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
}

# Incarca toate pachetele
lapply(packages, library, character.only = TRUE)


# SERII UNIVARIATE ****************************************

# Incarcare set de date
energie <- read_excel("C:/Users/Mihai/Desktop/Scoala/ASE - INFO ECO/AN III/Sem2/Serii de timp/Proiect/EnergieEolianaRomania.xlsx")
energie_ts <- ts(energie$generation_twh, start = c(2015, 1), frequency = 12)

# Convertire in format DATE
energie <- energie %>%
  mutate(date = as.Date(date, format = "%d-%m-%y")) %>%
  as_tsibble(index = date)

energie <- energie %>%
  mutate(date = yearmonth(date)) %>%   # transforma data in yearmonth pentru date lunare
  as_tsibble(index = date)

# CRONOGRAMA - avem sezonalitate
autoplot(energie_ts) +
  labs(title = "Generare energie eoliana in Romania",
       y = "TWh", x = "An") +
  theme_minimal()

# GRAFIC SEZONALITATE
ggsubseriesplot(energie_ts) +
  ylab("TWh") +
  ggtitle("Seasonal subseries plot: Wind energy generation in Romania")

# GRAFIC AUTOCORELARE
ggAcf(energie_ts) +
  labs(title = "Autocorrelation: Wind energy generation in Romania",
       y = "ACF", x = "Lag")

# DESCOMPUNERE ADITIVA
energie %>%
  model(classical_decomposition(generation_twh, type = "additive")) %>%
  components() %>%
  autoplot() +
  labs(title = "Classical additive decomposition: Wind energy generation in Romania",
       x = "Date", y = "TWh")


# STATIONALITATE

# Aplica testul ADF
adf_result <- adf.test(energie_ts)
print(adf_result) # p-value < 0.05 → respingem H0 => seria este stationara in medie

#checkresiduals(serie)

ArchTest(energie_ts) # p-value < 0.05 => seria este nestationara in varianta

# Radacina unitara
rw_none <- ur.df(energie_ts, type = "none", selectlags = c("AIC"))
summary(rw_none) # => nu e stationara  (1.18 < 2.5/1.9/1.6)

rw_drift <- ur.df(energie_ts, type = "drift", selectlags = c("AIC"))
summary(rw_drift) # => seria e stationara in jurul unei medii

rw_trend <- ur.df(energie_ts, type = "trend", selectlags = c("AIC"))
summary(rw_trend) # => seria este stationara cu trendul

#KPSS
energie_ts %>% ur.kpss() %>% summary() #0.018 < valorile critice => seria e stationara

# Philips-Perron
PP.test(energie_ts) #p-value < 0.01 => serie stationara

# Verificam nr de diferentieri pentru a fi stationara
ndiffs(energie_ts)

# ARIMA - inceput -----------------------
ggtsdisplay(energie_ts)

energie_adj <- energie_ts %>%
  stl(s.window = "periodic") %>%
  seasadj()

autoplot(energie_adj) +
  labs(title = "Generare energie eoliana in Romania",
       y = "TWh", x = "Luna") +
  theme_minimal()

ggsubseriesplot(energie_adj) +
  ylab("TWh") +
  ggtitle("Seasonal subseries plot: Wind energy generation in Romania")

ggtsdisplay(energie_adj)

# SARIMA -------------------------------------------------------
ggtsdisplay(energie_ts)
# Testam daca avem nevoide de diferentiere sezoniera
# Hegy
hegy.test(energie_ts) # as putea aplica o diferentiere sezoniera

# Canova Hansen
ch.test(energie_ts) # avem pattern sezonier stabil

# Diferentiem sezonier
energie_dif <- energie_ts %>% diff(lag=12) 
ggtsdisplay(energie_dif)

# ADF
rw_none <- ur.df(energie_dif, type = "none", selectlags = c("AIC"))
summary(rw_none) # => seria e stationara  (6.69 > 2.5/1.9/1.6)

rw_drift <- ur.df(energie_dif, type = "drift", selectlags = c("AIC"))
summary(rw_drift) # => seria e stationara in jurul unei medii

rw_trend <- ur.df(energie_dif, type = "trend", selectlags = c("AIC"))
summary(rw_trend) # => seria este stationara cu trendul

#KPSS
energie_dif %>% ur.kpss() %>% summary() #0.069 < valorile critice => seria e stationara

# Philips-Perron
PP.test(energie_dif) #p-value < 0.01 => serie stationara

# => seria este stationara la 99%

# AVEM SAR2 , SMA1 , AR3, MA3

fit1 <- Arima(energie_ts, order=c(3,0,3), seasonal=c(1,1,1))
coeftest(fit1) #potential

fit2 <- Arima(energie_ts, order=c(2,0,2), seasonal=c(1,1,1))
coeftest(fit2) #potential

fit3 <- Arima(energie_ts, order=c(2,0,2), seasonal=c(2,1,1))
coeftest(fit3) #nesemnificativ

fit4 <- Arima(energie_ts, order=c(3,0,3), seasonal=c(2,1,1))
coeftest(fit4) #nesemnificativ

fit5 <- Arima(energie_ts, order=c(2,0,2), seasonal=c(0,1,1))
coeftest(fit5) #potential

fit6 <- Arima(energie_ts, order=c(1,0,1), seasonal=c(0,1,0))
coeftest(fit6) #potential

summary(fit1)
summary(fit2) # are AIC AICC BIC mai mici
summary(fit5)
summary(fit6)

checkresiduals(fit2)

#Autocorelare
Box.test(residuals(fit2),lag=1,type='Lj')
Box.test(residuals(fit2),lag=2,type='Lj')
Box.test(residuals(fit2),lag=3,type='Lj')
Box.test(residuals(fit2),lag=4,type='Lj')
Box.test(residuals(fit2),lag=12,type='Lj')
Box.test(residuals(fit2),lag=24,type='Lj') #=> nu prezinta autocorelare

#Normalitate
jarque.bera.test(residuals(fit2)) # reziduuri normal distribuite, p-value > 0.1

#Heteroschedasticitate
ArchTest(residuals(fit2),lags=1)
ArchTest(residuals(fit2),lags=2)
ArchTest(residuals(fit2),lags=3)
ArchTest(residuals(fit2),lags=4)
ArchTest(residuals(fit2),lags=12)
ArchTest(residuals(fit2),lags=24) #=> nu avem heteroschedasticitate

# Prognoza
fit2 %>% forecast::forecast(h=12) %>% autoplot() + ylab('Twh') +
  theme_bw() +theme(plot.title = element_text(hjust = 0.5))

# Estimare ETS
fit_ets <- ets(energie_ts)
summary(fit_ets) 

fit_ets %>% forecast::forecast(h=12) %>%
  autoplot() +
  ylab("TWH") + 
  theme_bw() +theme(plot.title = element_text(hjust = 0.5))


# SARIMA forecast
fc_sarima <- forecast::forecast(fit2, h = 12)

# ETS forecast
fc_ets <- forecast::forecast(fit_ets, h = 12)

# Grafic combinat
autoplot(energie_ts) +
  autolayer(fc_sarima, series = "SARIMA", PI = FALSE) +
  autolayer(fc_ets, series = "ETS", PI = FALSE) +
  xlab("Timp") +
  ylab("TWh") +
  ggtitle("Prognoza combinată: SARIMA și ETS") +
  guides(colour = guide_legend(title = "Model")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))



# SERII MULTIVARIATE ****************************************
set_multi <- read_excel("C:/Users/Mihai/Desktop/Scoala/ASE - INFO ECO/AN III/Sem2/Serii de timp/Proiect/Set-Multivariat-EnergieCarbune.xlsx")

# Scatterplot
ggplot(data = set_multi) + 
  geom_point(mapping = aes(x = generation_twh, y = Valoare)) +
  xlab('Twh') +
  ylab('IPI Extractie carbune') + 
  ggtitle('Norul de puncte dintre indicele productiei industriale pentru extractia carbunelui si productia electricitatii pe baza de carbune')+
  theme_bw()
# avem corelatie direct proportionala


# Declaram variabilele de tip ts
twh_carbune <- ts(set_multi$generation_twh, start = c(2020,1), frequency = 12)
ipi_carbune <-ts(set_multi$Valoare, start = c(2020,1), frequency = 12)

ggtsdisplay(twh_carbune)
ggtsdisplay(ipi_carbune)

# Test Engle-Granger
# H0: seriile nu sunt cointegrate
# H1: seriile sunt cointegrate

# Pasul 1 - Testam stationaritatea seriilor
summary(ur.df(twh_carbune, type = "trend", selectlags = "AIC")) # serie stationara la 90%
summary(ur.df(diff(twh_carbune), type = "trend", selectlags = "AIC")) #serie stationara la 99%
summary(ur.df(ipi_carbune, type = "trend", selectlags = "AIC")) # serie stationara la 95%
summary(ur.df(diff(ipi_carbune), type = "trend", selectlags = "AIC")) # serie stationara la 99%

# Pas 2: aplicam testul de cointegrare
coint.test(y = twh_carbune,X = ipi_carbune,d = 0) # seriile nu sunt cointegrate
coint.test(y = ipi_carbune,X = twh_carbune,d = 0) # seriile nu sunt cointegrate

dset <- cbind(twh_carbune,ipi_carbune)
lagselect <- VARselect(dset, lag.max = 12, type = 'const')
lagselect$selection

# => APLICAM VAR

ggtsdisplay(diff(twh_carbune))
ggtsdisplay(diff(ipi_carbune))

dset_diff <- cbind(diff(twh_carbune),diff(ipi_carbune))
lagselect <- VARselect(dset_diff, lag.max = 12, type = 'const')
lagselect$selection

# Implementare VAR
model1 <- VAR(dset, p=1, type='const', season=NULL, exog=NULL) #radacinile sunt mai mici decat 1
summary(model1) #p-value < 0.05 => modelul este valid din punct de vedere statistic

model2 <- VAR(dset_diff, p=2, type='const', season=NULL, exog=NULL) #radacinile sunt mai mici decat 1
summary(model2) #p-value < 0.05 => modelul este valid din punct de vedere statistic


# Autocorelarea
Serial1 <- serial.test(model2, lags.pt = 12, type = 'PT.asymptotic')
Serial1 # pvalue > 0.1 nu avem autocorelare in reziduuri

Serial1 <- serial.test(model1, lags.pt = 12, type = 'PT.asymptotic')
Serial1 # pvalue > 0.1 nu avem autocorelare in reziduuri

# Heteroscedasticitate
Arch1 <- vars::arch.test(model2,lags.multi = 12,multivariate.only = TRUE)
Arch1 # pvalue > 0.1 modelul nu prezinta heteroschedasticitate la 99%

Arch1 <- vars::arch.test(model1,lags.multi = 12,multivariate.only = TRUE)
Arch1 # pvalue > 0.1 modelul nu prezinta heteroschedasticitate la 99%

# Normalitatea reziduurilor
Norm1 <- normality.test(model2, multivariate.only = TRUE)
Norm1 # pvalue JB > 0.1 reziduurile sunt normal distribuite

Norm1 <- normality.test(model1, multivariate.only = TRUE)
Norm1 # pvalue JB < 0.1 reziduurile nu sunt normal distribuite

# Testarea pentru rupturi in serie
Stability1 <- stability(model2,type = 'OLS-CUSUM')
plot(Stability1) # model stabil deoarece seriile noastre nu depasesc intervalul rosu

Stability1 <- stability(model1,type = 'OLS-CUSUM')
plot(Stability1) # model stabil deoarece seriile noastre nu depasesc intervalul rosu

# Cauzalitate Granger
# H0: valorile cu lag ale lui X, nu explica variatia in Y 
# H1: valorile cu lag ale lui X, explica variatia in Y

GrangerTWH <- causality(model2, cause='diff.twh_carbune.')
GrangerTWH # p-value < 0.1 => energia generata prin carbuni prezinta cauzalitate cu ipi pentru extractia carbunilor, 
# nu avem cauzalitate Granger

GrangerIPI <- causality(model2, cause='diff.ipi_carbune.')
GrangerIPI  # p-value < 0.1 => ipi pt extractia carbunilor prezinta cauzalitate pentru energia generata prin carbuni, 
# avem cauzalitate Granger
# => schimbarile in productia industriala pentru extractia carbunelui determina schimbarile in energia produsa din carbune


GrangerTWH <- causality(model1, cause='twh_carbune')
GrangerTWH # avem cauzalitate granger

GrangerIPI <- causality(model1, cause='ipi_carbune')
GrangerIPI # avem cauzalitate granger

# Functia de raspuns la impuls (IRF)
TWHirf <- irf(model2, impulse = 'diff.ipi_carbune.', response = 'diff.twh_carbune.', 
              n.ahead = 12, boot = TRUE, ci=0.90)
plot(TWHirf, ylab = 'TWH', main = 'Raspunsul energiei generate de carbuni la socurile indicelui de productie industriale pentru extractai carbunelui')

IPIirf <- irf(model2, impulse = 'diff.twh_carbune.', response = 'diff.ipi_carbune.', 
              n.ahead = 12, boot = TRUE, ci=0.90)
plot(IPIirf, ylab = 'IPI', main = 'Raspunsul indicelui de productie industriale pentru extractai carbunelui la socurile energiei generate de carbuni')


#---
TWHirf <- irf(model1, impulse = 'ipi_carbune', response = 'twh_carbune', 
              n.ahead = 12, boot = TRUE, ci=0.90)
plot(TWHirf, ylab = 'TWH', main = 'Raspunsul energiei generate de carbuni la socurile indicelui de productie industriale pentru extractai carbunelui')

IPIirf <- irf(model1, impulse = 'twh_carbune', response = 'ipi_carbune', 
              n.ahead = 12, boot = TRUE, ci=0.90)
plot(IPIirf, ylab = 'IPI', main = 'Raspunsul indicelui de productie industriale pentru extractai carbunelui la socurile energiei generate de carbuni')


# Descompunerea variantei
FEVD <- fevd(model2,n.ahead=12)
plot(FEVD)
FEVD

FEVD <- fevd(model1,n.ahead=12)
plot(FEVD)
FEVD

# Prognoza VAR
forecast <- predict(model2, n.ahead = 6, ci = 0.90)
plot(forecast, name = 'diff.ipi_carbune.')
plot(forecast, name = 'diff.twh_carbune.')

fanchart(forecast, names='diff.ipi_carbune.')
fanchart(forecast, names='diff.twh_carbune.')


forecast <- predict(model1, n.ahead = 6, ci = 0.90)
plot(forecast, name = 'ipi_carbune')
plot(forecast, name = 'twh_carbune')

fanchart(forecast, names='ipi_carbune')
fanchart(forecast, names='twh_carbune')

#legea lui okun se valideaza => atunci cand una creste si cealalta creste