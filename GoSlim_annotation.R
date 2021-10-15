#ANOTACION CON BASE DE DATOS GOSLIM
#GOSlim hace referencia a las versiones reducidas de las ontologias de GO que
#contienen un subconjunto de terminos biologicos.

#Autora: Celeste Moya Valera
#Fecha de creacion: 18/05/2021
#Asignatura: Nociones basicas de bioinformatica y genomica

#Mensaje al usuario
print("Anotacion de listas de genes con la informacion de la base de datos seleccionada.")
print("Bienvenido al Script de Anotacion de GOSlim")

#Cargamos la libreria BiomaRt de Ensembl
library(biomaRt)

#Lectura de las lista de datos
setwd("/home/ceromova/Desktop")
directorio = getwd()
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
  header=raw_file[1:4] #cabecera
  genes_list_raw=raw_file[5:length(raw_file)] #genes raw
  
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
  
  #Posibles nombres de los organismos
  Organism_Names=c("#homo_sapiens","#homo-sapiens","#homosapiens","#Homo_sapiens")
  
  #Planteamos la anotacion de los tres organismos posibles
  ##Homo sapiens
  ###Realizamos la anotacion con un DataSet de genes humanos (Homo sapiens)
  df_anotacion = data.frame()
  if(raw_file[3] %in% Organism_Names){
    print(paste("El organismo del archivo",temp_name,"es Homo sapiens."))
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    for (gen in F_GENES){
      total_anotacion <- getBM(attributes=c('hgnc_symbol','goslim_goa_accession',
                               'goslim_goa_description'), filters = 'hgnc_symbol', 
                               values = gen, mart = ensembl)
      df_anotacion = rbind(df_anotacion, total_anotacion)}}

  #Guardamos el nombre del fichero del alumno anotado
  nombre_alumno = substr(fichero, 0, nchar(fichero)-4)
  nombre_salida = paste0("ficheros_anotados/",nombre_alumno,"_anotado.txt")
  fichero_salida=file(nombre_salida)
  #writeLines(header,fichero_salida) #Opcion para poner encabezado
  fichero_salida = write.table(df_anotacion, file = nombre_salida, 
                               append=TRUE, col.names = TRUE, row.names = TRUE,
                               sep="\t",quote = FALSE)
}
print("Fin del Script.")

#*******************************************************************************
#* GENERACION DE GRAFICAS                                                     **
#*******************************************************************************

#library(ggplot2)
#library(plyr)

#vector = read.csv("/home/ceromova/Escritorio/goslim.txt",sep="\t",header=TRUE)
#df = as.data.frame(vector[3])
#contar = count(df, vars = NULL, wt_var = NULL)
#contar_filtro = subset(contar, contar$freq > 19)

#ggplot(data = contar_filt) + geom_histogram(aes(x=contar_filtro$freq))

#ggplot(contar_filtro, aes(ID, Count, fill=ONTOLOGY), split='ONTOLOGY') + geom_col()  + 
#  theme(axis.text.x=element_text(angle=-40, hjust=0)) + facet_grid(.~ONTOLOGY, scale="free_x")

