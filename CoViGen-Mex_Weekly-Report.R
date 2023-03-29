# UPDATE 2021-05-24: I added an optional xlsx output (uncomment if required)
# Started: 2021-05-20
# by Rodrigo García-López for Carlos Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 4.0.5 (2021-03-31) -- "Shake and Throw"
# The readxl library is required:
# install.packages("readxl")
# install.packages("writexl") # UPDATE 2021-05-24: Optional: required for xlsx output (uncomment if required)

# IMPORTANT NOTE: GISAID tables are expected to be present at sheet 1

# The script is used to merge the GISAID report table, the Patient metadata table and the Lineage data table into a single table for the InDRE and CoViGen-Mex reports.

# 3 inputs are mandatory
	# INPUT 1: PATIENT - .xlsx (expected at sheet1) table containing the patient metadata but no sample data. Most of the columns in the output are derived from this table. On 2021-05-20, it has the following columns:
		#Folio Interno	Edad (años)	Sexo	Estado	Municipio 	Fecha de toma	Tipo de muestra 	Tipo de paciente	Resultado	CT Gen 1	CT Gen 2	CT Gen 3	Embarazo	Semanas de gestación	Folio SINAVE	Folio SINOLAVE
	# INPUT 2: GISAID - .xlsx (expected at sheet2) table containing the following metadata (these should not vary) for each sample. On 2021-05-20, it has the following columns:
		# Submitter	FASTA filename	Virus name	Type	Passage details/history	Collection date	Location	Additional location information	Host	Additional host information	Sampling Strategy	Gender	Patient age	Patient status	Specimen source	Outbreak	Last vaccinated	Treatment	Sequencing technology	Assembly method	Coverage	Originating lab	Address	Sample ID given by originating laboratory	Submitting lab	Address	Sample ID given by the submitting laboratory	Authors	Comment	Comment Icon
	# INPUT 3: LINAJE - MexCov2 lineage .csv table containing mutations and predicted lineage for each sample. On 2021-05-20, it has the following columns:
		# Nombre de la secuencia","Longitud de la secuencia","Linaje Pangolin","Probabilidad","Clado Nextstrain","Mutaciones Nucleótidos Conocidas","Nuevas Mutaciones Nucleótidos","Mutaciones Aminoácido Conocidas","Nuevas Mutaciones Aminoácido","Deleciones Nucleótidos Conocidas","Nuevas Deleciones Nucleótidos","Deleciones Aminoácido Conocidas","Nuevas Deleciones Aminoácido","Inserciones Nucleótidos Conocidas","Nuevas Inserciones Nucleótidos

# MAIN OUTPUT (CoViGen-Mex report): A regular table depicting the following information per sample:
# Folio Interno	Edad (años)	Sexo	Estado	Municipio 	Fecha de toma	Tipo de muestra 	Tipo de paciente	Resultado	CT Gen 1	CT Gen 2	CT Gen 3	Embarazo	Semanas de gestación	Folio SINAVE	Folio SINOLAVE	Virus name	SampleIDProcess	Pangolin Clade	Clado Nextstrain	Mutaciones Aminoacidos Conocidas	Nuevas Mutaciones Aminoá¡cido	Deleciones conocidas	Deleciones nuevas

# SECONDARY OUTPUT (InDRE report): A regular table with the following items per sample:
# Folio Interno	Edad (años)	Sexo	Estado	Municipio 	Fecha de toma	Tipo de muestra 	Tipo de paciente	Resultado	CT Gen 1	CT Gen 3	Embarazo	Semanas de gestación	Folio SINAVE	Folio SINOLAVE	Nombre del genoma	Clado

# Mapa (en español)
# Para resumir, habría que generar un programa para unir tres tablas: GISAID+LINAJES+PACIENTES y sacar el reporte del InDRE y los datos para el reporte CoViGen-Mex
#
# Tenemos 3 tablas:
# 1.- PACIENTE - Tabla xlsx con los datos de paciente en español
# 	Xref:
# 		[Folio Interno] <-> [Sample ID given by originating laboratory] GISAID
# 	Export:
# 		Todos los campos
# 2.- GISAID - Tabla xlsx del reporte oficial de https://www.gisaid.org/ sobre los genomas
# 	Xref:
# 		[Virus name] <-> [Nombre de la secuencia] (LINAJES)
# 		[Sample ID given by originating laboratory] <-> [Folio Interno] (PACIENTES)
# 	Export:
# 		[Virus name] -> REP_InDRE
# 	Únicos útiles:
# 		[FASTA filename], [Originating lab], [Address] [Authors]
# 3.- LINAJES (MexCov2) - Tabla csv con predicción de linajes Pangolín, basados en las mutaciones encontradas
# 	Xref:
# 		[Nombre de la secuencia] <-> [Virus name] (GISAID)
# 	Export:
# 		[Linaje Pangolin] -> Rep_InDRE, Rep_CoViGen-Mex
# 		[Clado Nextrain] -> Rep_CoViGen-Mex
# 		[Mutaciones] -> Rep_CoViGen-Mex

# EXECUTION:
# Execute from bash as follows (please write full Rscript PATH if not an env variable):
# Rscript CoViGen-Mex_Weekly-Report.R <PATIENTE.xlsx> <GISAID.xlsx> <LINAJE.csv> <prefix_output> [optional:GISAID_sheet]

# Test in R:
# library("readxl")
# setwd("/home/rod/Documents/01_Projects/SARS/Vigilancia")
# in1 <- "InDRE/24Mar2021_EpiCoV_Lote5_INER_AnexoVariantes.xlsx"
# in2 <- "Metadata/GSAID/24Mar2021_EpiCoV_Lote5_INER_GISAID_metadata_97.xls"
# in3 <- "Metadata/MexCoV2_Resul/24Marzo2021_MexCov2.csv"
# prefix <- "Test/test"
# sheet <- 2

# Test in bash:
# Rscript CoViGen-Mex_Weekly-Report.R InDRE/24Mar2021_EpiCoV_Lote5_INER_AnexoVariantes.xlsx Metadata/GSAID/24Mar2021_EpiCoV_Lote5_INER_GISAID_metadata_97.xls Metadata/MexCoV2_Resul/24Marzo2021_MexCov2.csv Test2/output_test


 ### PRE-LOAD PARAMETERS AND DEFINE GLOBAL VARIABLES ###
library("readxl")
library("writexl") # UPDATE 2021-05-24: Optional: required for xlsx output (uncomment if required)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) { # at least, 4 arguments should be included: <prefix_output>  <set_name> <#_group_name_1> <#_group_name_n>
  stop("A minimum of 4 arguments are mandatory: Rscript
CoViGen-Mex_Weekly-Report.R <PATIENTE.xlsx> <GISAID.xlsx> <LINAJE.csv> <prefix_output> [optional:GISAID_sheet]", call.=FALSE)
}
in1 <- as.character(args[1]) # Get a string handle for the input xlsx Patient-metadata file
in2 <- as.character(args[2]) # Get a string handle for the input xlsx GSAID sample-metadata file
in3 <- as.character(args[3]) # Get a string handle for the input csv lineage per sample file
prefix <- as.character(args[4]) # Get a string handle for the output (it can contain a complete path but requires a seed name for the output file)
sheet <- as.numeric(args[5]) # Get a numeric handle for specifying if the GISAID data is in a different sheet
sheet <- ifelse(is.na(sheet),2,sheet) # set default GISAID value to sheet 2 unless the user specifies otherwise

 ### LOAD FUNCTIONS AND INPUTS ###
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	return(string)
}
get_column <- function(string, table){ # Gets a string and return a column number where header matches (not case-sensitive) or aborts if it not found
	query <- strip(string)
	col <- grep(query,strip(colnames(table)))
	if(length(col)!=1){stop(paste0("Abortando: No se encontro columna única en tabla ",substitute(table), " para query: ", substitute(query)), call.=FALSE)} # abort if none or > 1 matching item was found
	return(col)
}
get_uniq <- function(df1, df2){ # Get unique items in two vectors and print them and their table
	test <- setdiff(df1,df2);if(length(test)!=0){print(paste("ADVERTENCIA!!! Item",test,"de", paste(substitute(df1), collapse = '>'), ", no encontrado en", paste(substitute(df2), collapse = '>')))}
	test <- setdiff(df2,df1);if(length(test)!=0){print(paste("ADVERTENCIA!!! Item",test,"de", paste(substitute(df2), collapse = '>'), ", no encontrado en", paste(substitute(df1), collapse = '>')))}
}
check_dates <- function(inmat){ # Check for errors in the table
	dates <- inmat[,get_column("Fecha", inmat)] # get the date from patient data (as it may have errors)
	date_errors <- which(is.na(as.Date(dates)))
	if(length(date_errors)>0){print(paste("ADVERTENCIA!!! Error en Fecha de toma Tabla: PACIENTE| Línea:",date_errors))}
	return(date_errors)
}
# Load inputs according to sample format:
PACIENTE <- as.data.frame(read_excel(in1, sheet = 1))
check_dates(PACIENTE)
GISAID <- as.data.frame(read_excel(in2, sheet = sheet))
LINAJE <- read.csv(in3,header=T, skip=0, comment.char='',fill=F, check.names=FALSE)
subprefix <- basename(prefix)
# The original GISAID may have an additional header (two headers) in which case the table will result in a badly assignated types for columns. The following line fixes this in case it's present:
if(GISAID[1,1]=="Submitter"){colnames(GISAID) <- GISAID[1,];GISAID <- GISAID[-1,]; GISAID <- readr::type_convert(GISAID)} # If the second line is present, reassign header names, delete extra header and reassign column types
# I found some rare exception an extra column was added by mistake in the GISAID table. To remove these, an extra filter was added and a warning is printed to stdout:
 ### MAIN ###
# Define the x-ref columns:
xref_paciente_folio <- get_column("Folio Interno", PACIENTE)
xref_gisaid_folio <- get_column("Sample ID given by originating laboratory", GISAID)
xref_gisaid_sample <- get_column("Virus name", GISAID)
xref_linaje_sample <- get_column("Nombre de la secuencia", LINAJE)
linaje_clade <- get_column("Linaje Pangolin", LINAJE)

# Check if other required columns exist in table PACIENTE. Any that fail to be found will abort current execution, uncomment as required
# get_column("Edad", PACIENTE)
# get_column("Sexo", PACIENTE)
# get_column("Estado", PACIENTE)
# get_column("Municipio", PACIENTE)
# get_column("Fecha de toma", PACIENTE)
# get_column("Tipo de muestra", PACIENTE)
# get_column("Tipo de paciente", PACIENTE)
# get_column("Resultado", PACIENTE)
# get_column("CT Gen 1", PACIENTE)
# get_column("CT Gen 2", PACIENTE)
# get_column("CT Gen 3", PACIENTE)
# get_column("Embarazo", PACIENTE)
# get_column("Semanas de gestación", PACIENTE)
# get_column("Folio SINAVE", PACIENTE)
# get_column("Folio SINOLAVE", PACIENTE)

# Delete items with missing identifiers
PACIENTE <- PACIENTE[!is.na(PACIENTE[,xref_paciente_folio]),]
GISAID <- GISAID[!is.na(GISAID[,xref_gisaid_folio]),]
GISAID <- GISAID[!is.na(GISAID[,xref_gisaid_sample]),]
LINAJE <- LINAJE[!is.na(LINAJE[,xref_linaje_sample]),]
print("Items por tabla:")
print(paste(c("PACIENTE:","GISAID:","LINAJE:"),c(nrow(PACIENTE),nrow(GISAID),nrow(LINAJE))))
# PACIENTE[1,1]=1 #USED for testing, not required
# LINAJE[100,1]=100 #USED for testing, not required
# Test for missing items in cross-reference tables (unique items)
get_uniq(PACIENTE[,xref_paciente_folio], GISAID[,xref_gisaid_folio])
get_uniq(LINAJE[,xref_linaje_sample], GISAID[,xref_gisaid_sample])

# Subset the TABLES to keep only items present in both tables and place items in the same order.
# 1.- Start with the GISAID table (which has both cross-references)
GISAID <- GISAID[as.logical(GISAID[,xref_gisaid_folio]%in%PACIENTE[,xref_paciente_folio] * GISAID[,xref_gisaid_sample]%in%LINAJE[,xref_linaje_sample]),] # Subset items in all tables (order based on the first table)
print(paste("Items en las tres tablas:", nrow(GISAID)))
GISAID <- GISAID[order(GISAID[,xref_gisaid_folio]),] # now sort the table by the folio (patient number)
# 2.- Subset and sort table PACIENTE (depends on subsetting GISAID first)
PACIENTE <- PACIENTE[which(PACIENTE[,xref_paciente_folio]%in%GISAID[,xref_gisaid_folio]),] # Remove missing
PACIENTE <- PACIENTE[order(match(PACIENTE[,xref_paciente_folio],GISAID[,xref_gisaid_folio])),] # Sort
# 3.- Subset and sort table LINAJE (depends on subsetting GISAID first)
LINAJE <- LINAJE[which(LINAJE[,xref_linaje_sample]%in%GISAID[,xref_gisaid_sample]),] # Remove missing
LINAJE <- LINAJE[order(match(LINAJE[,xref_linaje_sample],GISAID[,xref_gisaid_sample])),] # Sort

# Assembly InDRE table
InDRE <- cbind(PACIENTE,"Nombre del genoma"=GISAID[,xref_gisaid_sample],"Clado"=LINAJE[,linaje_clade])
write.table(InDRE,paste(prefix,"AnexoVariantes.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) # output the InDRE table
write_xlsx(InDRE,paste(prefix,"AnexoVariantes.xlsx", sep="_")) # UPDATE 2021-05-24: Optional: required for xlsx output (uncomment if required)
# Assembly the CoViGen-Mex table
CoViGen <- cbind(PACIENTE,"Nombre del genoma"=GISAID[,xref_gisaid_sample],"Clado Pangolin"=LINAJE[,linaje_clade], "SampleIDProcess"=GISAID[,xref_gisaid_folio], LINAJE[,-c(1:4)])
write.table(CoViGen,paste(prefix,"Metadata_Pangolin.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write_xlsx(CoViGen,paste(prefix,"Metadata_Pangolin.xlsx", sep="_")) # UPDATE 2021-05-24: Optional: required for xlsx output (uncomment if required)
print("--- End of execution ---")
