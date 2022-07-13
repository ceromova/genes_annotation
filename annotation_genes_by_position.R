#Cargamos la libreria BiomaRt de Ensembl
library(biomaRt)

#Lectura de las lista de datos
setwd("C:/Users/celes/Desktop/")
directorio = getwd()
lista = list.files(pattern = "intersect")

#Creamos un directorio en el que guardar los archivos
dir.create("./ficheros_anotados")

for (fichero in lista){
  #Leemos el fichero de cada alumno
  raw_file = read.table(fichero,header=FALSE)
  colnames(raw_file) <- c("chr","start","end","log2FC","pvalue","exonid") #opcional
  
  ###Realizamos la anotacion con un DataSet de genes humanos (Homo sapiens)
  df_anotacion = data.frame()
  result = data.frame()
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  for (i in 1:nrow(raw_file)){
    total_anotacion <- getBM(attributes = c('hgnc_symbol','goslim_goa_accession',
                                            'goslim_goa_description'), 
                                 filters = c('chromosome_name','start','end'),
                                 values = list(raw_file[i,1],raw_file[i,2],raw_file[i,3]),
                                 mart = ensembl)
    #print(total_anotacion)
    
    if (dim(total_anotacion)[1] == 0){
      df_anotacion = rbind(df_anotacion, total_anotacion)
    }
    else{
    raw_info = data.frame(raw_file[i,1],raw_file[i,2],raw_file[i,3],raw_file[i,4],raw_file[i,5])
    colnames(raw_info) <- c("chr","start","end","log2FC","pvalue")
    
    total_anotacion = cbind(raw_info,total_anotacion)
    df_anotacion = rbind(df_anotacion, total_anotacion)
    }
    
  }
  
  #Guardamos el nombre del fichero del alumno anotado
  nombre_fichero = substr(fichero, 0, nchar(fichero)-4)
  nombre_salida = paste0("ficheros_anotados/",nombre_fichero,"_anotado.txt")
  fichero_salida=file(nombre_salida)
  #writeLines(header,fichero_salida) #Opcion para poner encabezado
  fichero_salida = write.table(df_anotacion, file = nombre_salida, 
                                 append=TRUE, col.names = TRUE, row.names = TRUE,
                                 sep="\t",quote = FALSE)
  
}
  
#*******************************************************************************
#* GENERACION DE GRAFICAS                                                     **
#*******************************************************************************

library(ggplot2)
library(dplyr)

#Directorio a los ficheros anotados
setwd(paste0(directorio,"/ficheros_anotados"))
vector = read.csv("C-D-intersect-sort_anotado.txt",sep="\t",header=TRUE)

vector = df_anotacion
df = as.data.frame(vector[8])
contar = as.data.frame(table(df))
contar_filtro = subset(contar, contar$Freq > 2) #Añadir el corte de frecuencia de aparicion del termino

#For palette choices:
#RColorBrewer::display.brewer.all()

#ggplot(data = contar_filtro) + geom_col(mapping = aes(df, Freq, fill = df)) + scale_fill_brewer(palette = "Pastel1") + labs(title = "Número de genes por término de GOSlim", x = "Términos de GOSlim", y = "Número de genes", fill="Términos de GOSlim") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(angle = 70, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5, color = "black", size = 20, face = "bold"))

ggplot(data = contar_filtro) + geom_col(mapping = aes(df, Freq, fill = df)) + labs(x = "GOSlim term", y = "Gene number", fill="GOSlim term") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(angle = 70, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5, color = "black", size = 20, face = "bold"))

