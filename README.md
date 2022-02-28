# genes_annotation
Gene annotation script in R using bioMart (ENSEMBL). The input of the program must be a gene list with no header in a TXT format. The output will depend on the attributes parameters used but the file format will be TXT. The final output will be keeped in a new folder called "ficheros_anotados". By default the output will show three different columns with the information of the header: hgnc_symbol	goslim_goa_accession	goslim_goa_description. Every row will have four columns, three according to the header and the first one which referes to the line number of the file.
