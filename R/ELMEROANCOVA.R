#'Analisis de covarianza.
#'
#'Realiza analisis de covarianza con una sola variable explicativa.
#'
#'@param var_resp (vector) datos respuesta.
#'@param var_explic (vector) datos explicativos.
#'@param bloque1:ny (vector) datos repuesta por bloque.
#'@param bloque1:nx (vector) datos explicativos por bloque.
#'@param totvar_explic1:nx suma de los datos explcativos por bloque.
#'@param totvar_resp1:nx suma de los datos respuesta por bloque.
#'@param trat1:ny (vector) datos repuesta por tratamiento.
#'@param trat1:nx (vector) datos explicativos por tratamiento.
#'@param totvar_explictrat1:nx suma de los datos explcativos por tratamiento.
#'@param totvar_resptrat1:ny suma de los datos respuesta por tratamiento.
#'@param gltotal_varexplic suma de los datos de la variable explicativa.
#'@param gtotal_varesp suma de los datos de la variable respuesta.
#'@param bloq numero de bloques.
#'@param trat numero de tratamientos.
#'@param fc_ cuadrado de las variable.
#'@param sc_bloq suma de cuadrados de los bloques
#'@param coeficientes coeficiente del modelo ajustado.
#'@param sc_trat suma de cuadrados de los tratamientos.
#'@param sctot_ suma de cuadrados totales.
#'@param sc_error suma de cuadrados de los errores totales.
#'@param scerror_ajust suma de cuadrados de los errores ajustada.
#'@param cf_cov obtencion del coeficiente de covarianza de la suma de cuadrados.
#'@param cuadmed_ cuadrados medios de las variables.
#'@param fcal_ obtencion de f calculada para tratamiento y covarianza.
#'@return el resultado del analisis de covarianza y el modelo ajustado de regresion lineal simple.
#'@export
ana_covarianza <- function(var_resp, var_explic) {

  var_resp <- Datos$tonha
  var_explic <- Datos$numplantas

  # Vectores de la variable respuesta por bloque
  bloque1y <- Datos$tonha[1:3]
  bloque2y <- Datos$tonha[4:6]
  bloque3y <- Datos$tonha[7:9]
  bloque4y <- Datos$tonha[10:12]

  # Vectores de la variable explicativa por bloque
  bloque1x <- Datos$numplantas[1:3]
  bloque2x <- Datos$numplantas[4:6]
  bloque3x <- Datos$numplantas[7:9]
  bloque4x <- Datos$numplantas[10:12]

  # Suma total de la variable explicativa por bloque
  totvar_explic1 <- sum(bloque1x)
  totvar_explic2 <- sum(bloque2x)
  totvar_explic3 <- sum(bloque3x)
  totvar_explic4 <- sum(bloque4x)

  # Suma total de la variable respuesta por bloque
  totvar_resp1 <- sum(bloque1y)
  totvar_resp2 <- sum(bloque2y)
  totvar_resp3 <- sum(bloque3y)
  totvar_resp4 <- sum(bloque4y)

  # Vectores de la variable respuesta por tratamiento
  trat1y <- Datos$tratIy[1:4]
  trat2y <- Datos$tratIIy[1:4]
  trat3y <- Datos$tratIIIy[1:4]

  # Vectores de la variable explicativa por bloque
  trat1x <- Datos$tratIx[1:4]
  trat2x <- Datos$tratIIx[1:4]
  trat3x <- Datos$tratIIIx[1:4]

  # Suma total de la variable explicativa por bloque
  totvar_explictrat1 <- sum(trat1x)
  totvar_explictrat2 <- sum(trat2x)
  totvar_explictrat3 <- sum(trat3x)

  # Suma total de la variable respuesta por bloque
  totvar_resptrat1 <- sum(trat1y)
  totvar_resptrat2 <- sum(trat2y)
  totvar_resptrat3 <- sum(trat3y)

  # Calcular el total de cada variable
  gtotal_varexplic <- sum(var_explic)
  gtotal_vareps <- sum(var_resp)

  # Calcular el gran total
  gtotal_varesp <- sum(bloque1y, bloque2y, bloque3y, bloque4y)
  gtotal_varexplic <- sum(bloque1x, bloque2x, bloque3x, bloque4x)

  bloq <- 4
  trat <- 3

  fc_xx <- (gtotal_varexplic)^2 / (bloq * trat)
  fc_yy <- (gtotal_varesp)^2 / (bloq * trat)
  fc_xy <- (gtotal_varexplic * gtotal_vareps) / (bloq * trat)

  scbloq_xx <- ((sum(totvar_explic1^2 + totvar_explic2^2 + totvar_explic3^2 + totvar_explic4^2)) / trat) - fc_xx
  scbloq_xy <- ((sum((totvar_explic1 * totvar_resp1) + (totvar_explic2 * totvar_resp2) + (totvar_explic3 * totvar_resp3) + totvar_explic4 * totvar_resp4)) / trat) - fc_xy
  scbloq_yy <- ((sum(totvar_resp1^2 + totvar_resp2^2 + totvar_resp3^2 + totvar_resp4^2)) / trat) - fc_yy

  sctrat_xx <- ((sum(totvar_explictrat1^2 + totvar_explictrat2^2 + totvar_explictrat3^2)) / bloq) - fc_xx
  sctrat_yy <- ((sum(totvar_resptrat1^2 + totvar_resptrat2^2 + totvar_resptrat3^2)) / bloq) - fc_yy
  sctrat_xy <- ((sum((totvar_explictrat1 * totvar_resptrat1) + (totvar_explictrat2 * totvar_resptrat2) + (totvar_explictrat3 * totvar_resptrat3))) / bloq) - fc_xy

  # Vectores que contienen el cuadrado de las variables
  cuad_xx <- var_explic^2
  cuad_yy <- var_resp^2
  cuad_xy <- var_explic * var_resp

  # Suma de cuadrados total
  sctot_xx <- sum(cuad_xx) - fc_xx
  sctot_yy <- sum(cuad_yy) - fc_yy
  sctot_xy <- sum(cuad_xy) - fc_xy

  # Suma de errores
  sc_errorxx <- (sctot_xx - scbloq_xx - sctrat_xx)
  sc_erroryy <- (sctot_yy - scbloq_yy - sctrat_yy)
  sc_errorxy <- (sctot_xy - scbloq_xy - sctrat_xy)

  sc_error <- sc_erroryy - (sc_errorxy^2 / sc_errorxx)

  # Suma de cuadrados ajustados
  scerror_ajustxx <- sctrat_xx + sc_errorxx
  scerror_ajustxy <- sctrat_xy + sc_errorxy
  scerror_ajustyy <- sctrat_yy + sc_erroryy

  # Obtencion del coeficiente de covariaza de la suma de cuadrados
  cf_cov <- sc_errorxy / sc_errorxx

  # Suma de cuadrados de error ajustado
  scerror_ajust <- scerror_ajustyy - (scerror_ajustxy^2 / scerror_ajustxx)

  # Suma de cuadrados de tratamientos
  sc_trat <- scerror_ajust-sc_error

  # Obtencion de cuadrados medios
  cuadmed_trat <- sc_trat/(trat-1)
  cuadmed_error <- sc_error/(((bloq-1)*(trat-1))-1)

  # Obtencion de f calculada de tratamiento y covarianza
  fcal_trat <- cuadmed_trat/cuadmed_error
  fcal_cov <- (sc_errorxy^2/sc_errorxx)/cuadmed_error
  return(list(
    fc_xx = fc_xx,
    fc_yy = fc_yy,
    fc_xy = fc_xy,
    scbloq_xx = scbloq_xx,
    scbloq_xy = scbloq_xy,
    scbloq_yy = scbloq_yy,
    sctrat_xx = sctrat_xx,
    sctrat_yy = sctrat_yy,
    sctrat_xy = sctrat_xy,
    sctot_xx = sctot_xx,
    sctot_yy = sctot_yy,
    sctot_xy = sctot_xy,
    sc_error = sc_error,
    scerror_ajustxx = scerror_ajustxx,
    scerror_ajustxy = scerror_ajustxy,
    scerror_ajustyy = scerror_ajustyy,
    y_ajust = y_ajust,
    scerror_ajust = scerror_ajust,
    sc_trat = sc_trat,
    cuadmed_trat = cuadmed_trat,
    cuadmed_error = cuadmed_error,
    fcal_trat = fcal_trat,
    fcal_cov = fcal_cov
  ))

}

# Directorio de trabajo
ruta <- "/Users/RODI ESAU/OneDrive/Documentos/EJEMP2-ANCOVA.CSV"
Datos <- read.csv(ruta)

library(ANCOVA)
# Llamada a la función ana_covarianza
resultados <- ana_covarianza(var_resp, var_explic)
