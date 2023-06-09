#' Cross-Over Desing
#'
#' Obtiene la Tabla de Analisis (ADEVA) para un Diseño Cross-over/Comuntativo o Cruzado.
#'
#' @param y (string) nombre de la variable respuesta.
#' @param trat (string) nombre de la variable que representan el tratamientos.
#' @param sujetos (string) nombre de la variable que representan los sujeto.
#' @param periodos (string) nombre de la variable que representan los periodos.
#' @param significancia (string) nombre de la variable que representa en grado de significancia.
#' @param data (\code{data.frame}) Tabla de datos en formato largo con los datos
#'      de los tratamientos, sujetos,periodos y de la variable respuesta.
#' @return Devuelve una tabla en formato \code{data.frame} con los calculos correspondientes
#' al analisis.
#' @export
#'
#' @examples
#' /dontrun{
#' #Directorio de trabajo
#' ruta <-"C:/Users/Usuario/Downloads/Prueba de mi paquete/DATOS 1.csv"
#'
#' #------------------------------------------------------------------
#' # Ejemplo 1:
#' df<-read.cvs(ruta)
#' # Cargo la libreria
#' library(DCruzado)
#' TablaAdeva(y ="respuesta",trat ="trat",sujetos = "sujetos",periodos = "periodos",significancia = 0.02,data = df)
#'
#' #----------------------------------------------------------------
#' #Ejemplo 2:
#' df<- data.frame(
#' vaca= c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4),
#' alimento= c("B","A","C","B","A","C","D","A","C","D","A","C","D","B","B","C"),
#' litros_leche= c(28,31,32.5,25,29,32,34,31,27,36,34,33,34,25,28,36),
#' periodos= c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
#'
#' TablaAdeva(y ="litros_leche",trat ="alimento",sujetos = "vaca",periodos = "periodos",significancia = 0.05,data = df)
#' }
TablaAdeva<-function(y, trat, sujetos, periodos,significancia =0.05, data){
  # Defino la variable respuesta y los tratamientos y sujetos y periodos  como factores

  y <- data[ , y]
  sujetos <- factor(data[, sujetos])
  trat <- factor(data[ , trat])
  periodos <- factor(data[,periodos])

  t <- nlevels(trat)
  a <- nlevels(sujetos)
  p <- nlevels(periodos)

  # factor de correccion
  suma_total <- sum(y)
  fc <-suma_total^2/(t*a)

  # suma de los cuadrados totales
  sc <- sum(y^2)-fc
  gl <- (t*a)-1


  # suma de los cuadrados por tratamiento
  sumaxtratamiento <- tapply(y, INDEX = trat, FUN=sum)
  n_sujetos <- tapply(y, INDEX = trat, FUN = length)
  sc_trat <-sum((sumaxtratamiento^2)/n_sujetos) - fc
  gl_trat <- t-1
  cm_trat  <- sc_trat/gl_trat

  #suma de los cuadrados por sujetos
  sumaxsujetos <-tapply (y, INDEX=sujetos, FUN=sum)
  n_trat <-tapply(y, INDEX=sujetos, FUN=length)
  sc_sujetos <-sum((sumaxsujetos^2)/n_trat) - fc
  gl_sujetos <-a-1
  cm_sujetos <-(sc_sujetos/(gl_sujetos))

  #suma de los cuadrados por periodos
  sumaxperiodos <-tapply (y, INDEX=periodos, FUN=sum)
  n_trat <-tapply(y, INDEX=periodos, FUN=length)
  sc_periodos <-sum((sumaxperiodos^2)/n_trat)-fc
  gl_periodos <-t-1
  cm_periodos <-sc_periodos/(gl_periodos)

  # calculo de SCe. exp
  SCe.exp <-sc-sc_trat-sc_sujetos-sc_periodos
  gl_SCe.exp <-gl-gl_trat-gl_sujetos-gl_periodos
  cm_SCe.exp <-SCe.exp/gl_SCe.exp

  #valores de F_trat
  f_trat <-cm_trat/cm_SCe.exp

  #valores de F_sujetos
  f_sujetos <-cm_sujetos/cm_SCe.exp

  #valores de F_periodos
  f_periodos <-cm_periodos/cm_SCe.exp


  #FTabulados
  ft_trat<- qf(1-significancia,gl_trat,gl_SCe.exp)
  ft_sujetos<- qf(1-significancia,gl_sujetos,gl_SCe.exp)
  ft_periodos<- qf(1-significancia,gl_periodos,gl_SCe.exp)

  #P_values
  p_valuetrat<-pf(f_trat,gl_trat,gl_SCe.exp,lower.tail = FALSE)
  p_valuesujetos<-pf(f_sujetos,gl_sujetos,gl_SCe.exp,lower.tail = FALSE)
  p_valueperiodos<-pf(f_periodos,gl_periodos,gl_SCe.exp,lower.tail = FALSE)

  #creamos un data.frame para ordenar los datos
  tabla <-data.frame(FV=c("Tratamientos","Sujetos","Periodos","Error Experimental","Total"),
                     GL=c(gl_trat,gl_sujetos,gl_periodos,gl_SCe.exp,gl),
                     SC=c(sc_trat,sc_sujetos,sc_periodos,SCe.exp,sc),
                     CM=c(cm_trat,cm_sujetos,cm_periodos,cm_SCe.exp,NA),
                     FC=c(f_trat,f_sujetos,f_periodos,NA,NA),
                     FT=c(ft_trat,ft_sujetos,ft_periodos,NA,NA),
                     P_Value=c(p_valuetrat,p_valuesujetos,p_valueperiodos,NA,NA))
  rownames(tabla) <-NULL
  ADEVA <-format(tabla)
  ADEVA[is.na(tabla)] <-""
  return(ADEVA)
}
