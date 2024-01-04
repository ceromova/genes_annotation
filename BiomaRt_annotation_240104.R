#ANOTACION CON BASE DE DATOS GOSLIM
#GOSlim hace referencia a las versiones reducidas de las ontologias de GO que
#contienen un subconjunto de terminos biologicos.

#Autora: Celeste Moya Valera
#Fecha de creacion: 18/05/2021

#Mensaje al usuario
print("Anotacion de listas de genes con la informacion de la base de datos seleccionada.")
print("Bienvenido al Script de Anotacion de GOSlim")

#Cargamos la libreria BiomaRt de Ensembl
library(biomaRt)

#Lectura de las lista de datos

# NOTA: recuerda cambiar el directorio
directorio = "C:/Users/celes/Desktop"
setwd(directorio)
lista = list.files(pattern = "\\.txt")

#Creamos un directorio en el que guardar los archivos
dir.create("./ficheros_anotados")

#Leemos el archivo y realizamos la anotacion
for (fichero in lista){
  #Guardamos el nombre del fichero
  temp_name=fichero
  
  #Comprobamos la codificacion de cada archivo
  cod = readr::guess_encoding(fichero)
  codificacion = toString(cod[1,1])
  con = file(fichero, encoding = codificacion)
  
  #Leemos el fichero de cada alumno
  raw_file = readLines(con, warn=FALSE)
  close(con)
  
  #Eliminamos las lineas en blanco
  empty_line=raw_file==""
  raw_file=raw_file[empty_line==FALSE]
  
  #Separamos la informacion del raw_file
  genes_list_raw=raw_file[0:length(raw_file)] #genes raw
  
  #Comprobamos si hay "///"
  if(TRUE %in% grepl("///", genes_list_raw)){
    f_genes_procesados=c()
    for (k in 1:length(genes_list_raw)){
      pattern=grepl("///", genes_list_raw[k])
      ##Si el nombre del gen contiene "///"
      if(pattern==TRUE){
        temp_split=strsplit(genes_list_raw[k], "///")
        f_genes_procesados=append(f_genes_procesados, temp_split[1])}
      else {
        ##Si el nombre del gen NO contiene "///" lo anyadimos directamente
        f_genes_procesados=append(f_genes_procesados,genes_list_raw[k])}}}
  else{
    ##Si NO hay ningun "///" no lo procesamos
    f_genes_procesados=na.omit(genes_list_raw)}
  
  #Pasamos a mayusculas los genes
  F_GENES=toupper(f_genes_procesados)

  ###Realizamos la anotacion con un DataSet de genes humanos (Homo sapiens)
  df_anotacion = data.frame()
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  for (gen in F_GENES){
      # Modificar los atributos para que salga el output requerido
      total_anotacion <- getBM(attributes=c('hgnc_symbol','goslim_goa_accession',
                                            'goslim_goa_description'), filters = 'hgnc_symbol', 
                               values = gen, mart = ensembl)
      df_anotacion = rbind(df_anotacion, total_anotacion)}
  
  #Guardamos el nombre del fichero del alumno anotado
  nombre_fichero = substr(fichero, 0, nchar(fichero)-4)
  nombre_salida = paste0("ficheros_anotados/",nombre_fichero,"_anotado.txt")
  fichero_salida=file(nombre_salida)
  fichero_salida = write.table(df_anotacion, file = nombre_salida, 
                               append=TRUE, col.names = TRUE, row.names = TRUE,
                               sep="\t",quote = FALSE)
}
print("Fin del Script.")

#*******************************************************************************
#* GENERACION DE GRAFICAS                                                     **
#*******************************************************************************

library(ggplot2)
library(dplyr)

vector = read.csv("C:/Users/celes/Desktop/gene_list_info/ficheros_anotados/A_B_anotado.txt",sep="\t",header=TRUE)
df = as.data.frame(vector[3])
contar = as.data.frame(table(df))
contar_filtro = subset(contar, contar$Freq > 0) #Añadir el corte de frecuencia de aparicion del termino

#For palette choices:
RColorBrewer::display.brewer.all()

ggplot(data = contar_filtro) + geom_col(mapping = aes(df, Freq, fill = df)) + scale_fill_brewer(palette = "Pastel1") + labs(title = "Número de genes por término de GOSlim", x = "Términos de GOSlim", y = "Número de genes", fill="Términos de GOSlim") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(angle = 70, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5, color = "black", size = 20, face = "bold"))