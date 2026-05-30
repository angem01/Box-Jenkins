#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                           UNIVERSIDAD NACIONAL DE COLOMBIA
#                   Facultad de Ciencias Económicas | 2026 - 01
#
#      Metodología Box-Jenkins para la identificación, estimación y pronóstico de
#                           series de tiempo univariadas
#                                  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
install.packages("pacman")
install.packages("mFilter")
install.packages("urca")
install.packages("FinTS")
install.packages("fable")
install.packages("fabletools")
install.packages("tsibble")

# Limpiamos el entorno 

rm(list = ls())
dev.off()

#_____________________________________________________________________________________#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 0. Instalación de Paquetes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(pacman)
 

pacman::p_load(
  
  forecast,   # Para hacer pronósticos con modelos arima
  lmtest,     # Significancia individual de los coeficientes ARIMA
  urca,       # Prueba de raíz unitaria
  tseries,    # Para estimar modelos de series de tiempo y hacer pruebas de supuestos
  stargazer,  # Para presentar resultados más estéticos
  psych,      # Para hacer estadísticas descriptiva
  seasonal,   # Para desestacionalizar series
  aTSA,       # Para hacer la prueba de efectos ARCH
  astsa,      # Para estimar, validar y hacer pronósticos para modelos ARIMA/SARIMA
  xts,        # Para utilizar objetos xts 
  tidyverse,  # Conjunto de paquetes (incluye dplyr y ggplot2)
  readxl,     # Para leer archivos excel 
  car,        # Para usar la función qqPlot
  mFilter,    # Para aplicar el Filtro Hodrick-Prescott
  quantmod,    
  sandwich,
  FinTS,
  # Paquetes del tidyverts
  
  fable,      # Forma moderna de hacer pronóstiocs en R (se recomienda su uso)  
  tsibble,    # Para poder emplear objetos de series de tiempo tsibble
  feasts      # Provee una colección de herramientas para el análisis de datos de series de tiempo 
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                         METODOLOGÍA BOX-JENKINS                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 1. Primer paso: Identificación ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~ CARGAR DATOS ~~~#

# Se cargan las series de tiempo
base_fred <- read_excel("Datos/COLPRMNVG01IXOBSAM_1_.xlsx")

glimpse(base_fred)
View(base_fred)

#~~~ TRANSFORMACIÓN A DATOS TS y XTS ~~~#

# Convertir la serie en un objeto tipo ts

proman = ts(base_fred$COLPRMNVG01IXOBSAM, start = 2011, frequency = 12) # 4 observaciones/año

t = as.vector(t(base_fred$COLPRMNVG01IXOBSAM))
ts = ts(t[121:216], start = c(2011), frequency = 12)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### 1.1. Análisis gráfico ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~ GRÁFICAS DE LAS SERIE ~~~#


x11() 

plot.ts(proman, xlab="Año",ylab="index",
        main="Producción: Manufactura: Bienes de inversión: Total para Colombia",
        sub = "2011-2018"
        ,lty=1, lwd=1, col="blue")



#~~~ APLICACIÓN DEL FILTRO HODRICK-PRESCOTT ~~~#



hpf = hpfilter(base_fred$COLPRMNVG01IXOBSAM, freq = 14400) 

#Componentes de la prueba
names(hpf)

#La parte ciclica de la serie de tiempo
hpf$cycle
ciclo = ts(hpf$cycle, start = 2011, frequency = 12)

#La parte tendencial de la serie de tiempo
hpf$trend
trend = ts(hpf$trend, start = 2011, frequency = 12)

#El valor de lambda correspondiente a la periodicidad de la ts
hpf$lambda

#La serie original
base_fred$COLPRMNVG01IXOBSAM


#Grafico con los componentes ciclico y tendencial de la serie 
x11()
#Viendo el componente tendencial Claramente la tendencia no es plana y no es independiente  
#del tiempo
plot.ts(trend)

plot.ts(ciclo)
plot(hpf)


#~~~ GRÁFICOS DE LAS FAC Y FACP ~~~#

lags=25

x11()
par(mfrow=c(1,2))
acf(proman, lag.max = lags, plot=T, lwd=2,xlab='',main='FAC') 
pacf(proman,lag.max=lags,plot=T,lwd=2,xlab='',main='FACP')
par(mfrow=c(1,1))

# Proceso no  estacionario

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### 1.2. Prueba Dickey Fuller aumentada (ADF) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Prueba con trend

adf.trend_proman = ur.df(proman, type="trend", lags = 5)
plot(adf.trend_proman)
summary(adf.trend_proman) #No es estacionaria en tendenicia

# Prueba con drift

adf.drift_proman= ur.df(proman, type="drift", lags = 5)
plot(adf.drift_proman)
summary(adf.drift_proman) #No es estacionaria sin tendencia y con deriva


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### 1.3. Transformación para volver estacionaria la serie #### 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

d.proman = diff(proman)
l.proman = log(proman)
dl.proman = diff(log(proman))*100

x11()
par(mfrow=c(2,2))
plot.ts(proman, xlab="",ylab="", main="En nivel",lty=1, lwd=2, col="blue")
plot.ts(d.proman, xlab="",ylab="", main="Diferenciada",lty=1, lwd=2, col="red")
plot.ts(l.proman, xlab="",ylab="", main="Logaritmo",lty=1, lwd=2, col="green")
plot.ts(dl.proman, xlab="",ylab="", main="Logaritmo en diferencia",lty=1, lwd=2, col="pink")

adf.trend_d.proman = ur.df(d.proman, type="trend", lags = 5)
plot(adf.trend_d.proman)
summary(adf.trend_d.proman)

#Es suficiente con la diferencia 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### 1.4. Identificación Modelo Arima ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


AR.m <- 6 #Supondremos que el rezago autorregresivo máximo es 6 (pmax)
MA.m <- 6 #Supondremos que el rezago de promedio móvil máximo es 6. (qmax)


#Esta función selecciona el modelo ARIMA con el menor criterio de información

# ¿Qué hace? realiza mediante permutación distintos modelos ARMA y sintetiza
# el calculo de los criterios de información en un Data Frame. (SOLO ML)


arma_seleccion_df = function(ts_object, AR.m, MA.m, d, bool_trend, metodo){
  
  index = 1
  df = data.frame(p = double(), d = double(), q = double(), AIC = double(), BIC = double())
  for (p in 0:AR.m) {
    for (q in 0:MA.m)  {
      fitp <- arima(ts_object, order = c(p, d, q), include.mean = bool_trend, 
                    method = metodo)
      df[index,] = c(p, d, q, AIC(fitp), BIC(fitp))
      index = index + 1
    }
  }  
  return(df)
}

#~~~ FUNCIÓN PARA SELECCIONAR ARIMA POR MENOR AIC ~~~#

arma_min_AIC = function(df){
  df2 = df %>% 
    filter(AIC == min(AIC))
  return(df2)
}

#~~~ FUNCIÓN PARA SELECCIONAR ARIMA POR MENOR BIC ~~~#


arma_min_BIC = function(df){
  df2 = df %>% 
    filter(BIC == min(BIC))
  return(df2)
}

# Para nuestro caso D = 0 ya que ya hemos diferenciado.

# Usaremos la función que hemos creado denominada arma_seleccion_df para escoger 
# el ARIMA a usar, con los máximos rezagos que hemos fijado (p = 6, q = 6).

mod_d1_proman = arma_seleccion_df(d.proman, AR.m, MA.m, d = 0, TRUE, "ML")

# Veamos los criterios.
View(mod_d1_proman)

# Selecciono el mejor modelo según menor valor de los criterios AIC y BIC.

min_aic_proman = arma_min_AIC(mod_d1_proman)
min_aic_proman # ARIMA (4,1,3)

min_bic_proman = arma_min_BIC(mod_d1_proman)
min_bic_proman # ARIMA (0,1,1)



# 2. Método automático (auto.arima) 

auto.arima(d.proman, method = "ML") #NO SE USA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 2. Segundo paso: Estimación ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Existen 3 métodos de estimación para la función arima:

## ML: Máxima verosimilitud (más preciso y la mejor opción para bases pequeñas)
## CSS: (más veloz generalmente, usado en bases de datos grandisimas).
## CSS-ML: Una combinación de ambas.

# Siendo CSS  (Coherent Source Separation) 

# Sin embargo ML puede no converger y CSS puede no hacer estimaciones lo 
# suficientemente precisas y arrojar error. 


# La estimación de un modelo arima se puede realizar 3 funciones distintas: 

## arima: Paquete stats, es la más usual emplear y predicción con forecast
## Arima: Paquete forecast que es básicamente un wrapper de la función arima.
## sarima: Paquete astsa

# Tanto la función Arima como la función sarima están construidas sobre la 
# función arima(), por lo que es posible modelar estacionalidad en cada alternativa

# Estimamos mediante Maxima verosimilitud no incurrir en imprecisiones


#Colocamos la serie;orden; ver si colocamos la diferenciacion en el vector de orden (el de la mitad)
# ML- Maxima verosimilitud

arima_0.1.1_dproman = arima(d.proman, order = c(0,0,1), 
                             include.mean = F, method = "ML")


arima_0.1.1_proman = arima(proman, order = c(0,1,1), 
                            include.mean = T, method = "ML")
# Sintetizamos los resultados
summary(arima_0.1.1_proman) # modelamiento ARIMA(0,0,1)



arima_4.1.3_dproman = arima(d.proman, order = c(4,0,3), 
                            include.mean = F, method = "ML")


# Sintetizamos los resultados
summary(arima_4.1.3_dproman)

# Stargazer
#stargazer(arima_0.0.1_constperm, arima_1.0.0_constperm,
#          column.labels=c("ARIMA(0,0,1)", "ARIMA(1,0,0)"),
#          keep.stat=c("n","rsq"), 
#          type = "text", style = "aer") # Recordemos que podemos obtener salida
                                        # LaTeX.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 3. Tercer paso: Validación de supuestos ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Es importante verificar los supuestos de nuestro modelo ARIMA. se debe ver que
# los residuales estimados se comporten como un ruido blanco. Es decir, que la 
# media de los residuales sea cero, la varianza constante y la covarianza sea cero.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### 3.1. No autocorrelación de los errores ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Se dice que la cantidad "ideal" de lags es un cuarto de la muestra
lags.test = length(d.proman)/4;lags.test

#--> ARIMA(1,1,2)

# Argumento gráfico

X11()

res_arima_0.1.1_dproman = residuals(arima_0.1.1_dproman)
par(mfrow=c(1,2))



acf(res_arima_0.1.1_dproman,lag.max=24,plot=T,lwd=1,xlab='',
    main='ACF residuales (0,1,1)') 
pacf(res_arima_0.1.1_dproman,lag.max=24,plot=T,lwd=1,xlab='',
     main='ACF al cuadrado residuales (0,1,1)')
par(mfrow=c(1,1))

# Hay rezagos significativos, mal indicio para el supuesto.

# Pruebas formales:

#~~ BOX-PIERCE TEST ~~# 

#Ho = No autocorrelación 
#Ha = Hay autocorrelación


Box.test(res_arima_0.1.1_dproman, lag=lags.test, type = c("Box-Pierce")) # Rechazamos H0
Box.test(res_arima_0.1.1_dproman, lag=20, type='Box-Pierce') 


#~~ LJUNG-BOX ~~#

Box.test(res_arima_0.1.1_dproman, lag=lags.test, type = c("Ljung-Box")) 
Box.test(res_arima_0.1.1_dproman, lag=20, type='Ljung-Box') 

# ARIMA(0,1,1) No hay autocorrelacion en 20 rezagos


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### 3.2. Homocedasticidad de los residuales ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# La prueba ARCH nos dice si los residuales son homocedasticos.


#Ho = Homocedasticidad 
#Ha = heterocedasticidad


# Hay dos formas de hacer la prueba: Un test Pormenteau y un Test tipo 
# multiplicadores de Lagrange.

#--> ARIMA(1,1,2)

arch_dproman_arima_0.1.1 = arch.test(arima_0.1.1_dproman, output=TRUE)



#Sin embargo, si se desea obtener un único p-value para un número de lags
#en especifico se puede utilizar: 

#Hallamos los residuos 
residuos <- residuals(arima_0.1.1_dproman)

#Realizamos la prueba

ArchTest(residuos, lags = 20)

#Como podemos observar, rechazamos la hipótesis nula a favor de heterocedasticidad

# Vamos a graficar la ACF y PACF de residuales al cuadrado del modelo ARIMA(2,0,4)

# Argumento gráfico.

X11()

par(mfrow=c(1,2))
acf(res_arima_0.1.1_dproman^2,lag.max=20,plot=T,lwd=2,xlab='',main='ACF residuales al cuadrado') 
pacf(res_arima_0.1.1_dproman^2,lag.max=20,plot=T,lwd=2,xlab='',main='PACF residuales al cuadrado')
par(mfrow=c(1,1))

# Con esto concluimos que no existe heterocedasticidad, se cumple el supuesto.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### 3.3. Normalidad en los residuales ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Veremos si los residuales se comportan de manera normal (distribución normal)
# para ello usaremos un argumento gráfico como el QQ-plot y la prueba de normalidad
# Jarque-Bera.


#--> ARIMA(1,0,0)
x11()
qqPlot(res_arima_0.1.1_dproman, ylab = "ARIMA(0,1,1)")

# Colas pesadas, no se cumple el supuesto.

# Prueba formal: Jarque-Bera Test


#Ho = Normalidad
#Ha = No hay normalidad


jarque.bera.test(res_arima_0.1.1_dproman) # Se rechaza H0, no hay normalidad. 
# No se cumple el supuesto de normalidad.
#Conclusiones
#ARIMA(1,1,2) Cumple todos los supuestos a excepción del de normalidad
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### 4.1. Pronósticos futuros ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Pronóstico 10 pasos adelantes, valores numéricos - modelo(0,1,1)
summary(arima_0.1.1_dproman)

forescast_arima_0.1.1_proman = arima_0.1.1_proman_fable %>% 
  forecast(h = 10, bootstrap = TRUE, times = 10000); forescast_arima_0.1.1_proman

pronósticos_0.1.1_proman = forescast_arima_0.1.1_proman %>% 
  hilo(level = c(80, 90, 95)); pronósticos_0.1.1_proman

View(pronósticos_0.1.1_proman)

x11()
forescast_arima_0.1.1_proman %>% 
  autoplot(proman_tsibble) + ggtitle("Pronósticos producción manufacturera de inversión - Colombia") + ylab("Indice (Base 2015 = 100") + xlab("mes") + theme_light()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                   FIN DEL CÓDIGO                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#