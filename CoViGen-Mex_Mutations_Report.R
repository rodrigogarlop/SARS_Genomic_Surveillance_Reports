# Update 2021-05-27 3: The new tables had fewer genes. I added exceptions to avoid errors
# Update 2021-05-27 2: The new tables had empty lines. I added a workaround for this issue.
# Update 2021-05-27 1: Newer versions of the input report have at signs (@) instead of asterisks (*) to mark historical mutations, so I updated the script accordingly
# Started: 2021-05-27
# by Rodrigo García-López for Carlos Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
# It should run with R-base libraries

# The script was created as a complementary summary from the CoViGen-Mex mutation tables (total mutations in current set and in historical sets)

# Input tables must have these header names
# name	gene	position	count	percent	nt_position

# Test in R:
# setwd("/home/rod/Documents/01_Projects/SARS/Vigilancia")
# csv <- "Metadata/MexCoV2_Resul/21Abril2021_amChartsAminoMut.csv"
# prefix <- "Test2/output2"
# name_date <- "2021-05-19"
# subprefix <- basename(prefix)
# past_linages <- "Test/5_fechas_linaje.tsv"
# variants <- "Test/VOC_VOI.tsv"

# Test in bash:
# Rscript CoViGen-Mex_Mutations_Report.R Metadata/MexCoV2_Resul/24Marzo2021_MutationAmino.csv Test2/output3 2021-05-19

# Run as follows:
# Rscript CoViGen-Mex_Mutations_Report.R <input_xlsx> <prefix_output> <date>

 ### PRE-LOAD PARAMETERS AND DEFINE GLOBAL VARIABLES ###
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) { # at least, 3 arguments should be included: <input_csv> <prefix_output>
  stop("A minimum of 3 arguments are mandatory: Rscript CoViGen-Mex_Mutations_Report.R <input_csv> <prefix_output> <name_date>", call.=FALSE)
}
csv <- as.character(args[1]) # Get a string handle for the input csv file
prefix <- as.character(args[2]) # Get a string handle for the output (it can contain a complete path but requires a seed name for the output file)
name_date <- as.character(args[3]) # Get a string handle for the name of the set (for graphs)
test <- readLines(csv) # UPDATE 2021-05-27: Newer files had some issues when loading. It turned out to be a list of empty lines (having only ",,,,,")
df <- read.table(text=test[!grepl(",,,,,",test)], sep=",",header=T, skip = 0, comment.char='',quote="",fill=F, check.names=F, row.names=1)

 ### LOAD FUNCTIONS ###
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
only_first_upcase <- function(string) { # Make only the first letter an uppercase character
  string <- strip(string)
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  return(string)
}
old_and_new <- function(inmat,grp="All"){ # split the table into old and new mutations based on rownames with and w/o * IMPORTANT: This assumes that no other categories exist historical records have a *. Group variable is used to retrieve only one group. It defaults to "All", where all items are included.
	get_column("percent",inmat) # abort if no percentages are present
	vect <- inmat[,get_column("gene",inmat)] # extract the gene column if present (abort if not)
	subset <- which(strip(vect)==strip(grp)) # get rows to collect. Make it case-insensitive
	if(length(subset)>0){inmat <- inmat[subset,]}else{paste("ATENCION!!! No se encuentran coincidencias para el gen:",grp, "... ignorando");break} # subset if present, if not, abort the cycle (skip it)
	hist <- inmat[grep("\\@",rownames(inmat)),] # Get the cummulative muttations
	curr <- inmat[grep("\\@",rownames(inmat),invert=TRUE),] # as well as the current ones
	if(nrow(hist)!=nrow(curr)){print(paste("ATENCION!!! Numero discrepante de mutaciones actuales e históricas en grupo:", grp));print(paste(" Mutaciones en histórico:",nrow(hist))); print(paste(" Mutaciones actuales:",nrow(curr)))} # Print a warning if different current and historical mutation numbers differ (not an error)
	rownames(hist) <- gsub(" ", "", gsub("\\@","",rownames(hist))) # Remove * and spaces
	rownames(curr) <- gsub(" ", "", rownames(curr)) # Remove spaces
	outmat <- merge(curr,hist,by="row.names",suffixes=c("-curr","-hist"), all = TRUE)
	colnames(outmat)[1]="Mutación"
	plotdata <- outmat[,rev(grep("percent",strip(names(outmat))))] # Subset historic and current percentages
	rownames(plotdata) <- outmat[,1];plotdata[is.na(plotdata)]=0
	out <- list(outmat,plotdata)
	names(out) <- c("Merged", "Counts")
	return(out)
}
gene_barplot <- function(inmat,grp="ALL"){ # Create barplot from the percent two-column matrix produced by old_and_new() and the name that will be used
	color <- c("blue","steelblue2")
	pdf(paste(prefix,"Mutacion",grp,"barplot.pdf", sep= "-"), width = 18, height = 8) # make a landscape plot
		name_size <- ifelse(nrow(inmat)>17,0.8,1)
		par(las=1)
		nmax <- ceiling(max(inmat)*.1)*10 # ceiling function to the nearest multiple of 10 to plot complete axes
		barplot(as.matrix(t(inmat)),beside=TRUE,main=paste("Mutaciones en",grp),col=color,border=NA, cex.axis=1.5, cex.names=name_size, ylab="% Total", cex.lab=1.5, cex.main=1.5,ylim=c(0,nmax))
# 		axis(2,las=1,at=seq(0,100,10))
		abline(h=0)
		legend("topleft",pch=15, col=color, legend=c("% en BD",name_date),bty = "n",cex=1.5) # legend
	dev.off()
}

analyze_gene <- function(inmat,grp="ALL"){ # Create the corresponding plot and table for the target gene
	if(!grepl(grp,gene_list)){return(paste("El gen", grp, "no se encontró ... ignorando"))} # Update 2021-05-27: Some newer tables have only the S gene, I added this checkpoint to deal with this. It uses global variable gene_list
	gene <- old_and_new(inmat, grp)
	gene_barplot(gene$Counts, grp)
	write.table(gene$Merged,paste(prefix,"Mutacion",grp,"table.tsv", sep= "-"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) # output the date vs type of patient table
}

remove_double_at_signs <- function(inmat){ # This was created to deal with a special case where repeated items are found in historical data, one of which has double asterisks and is the one removed here
	outmat <- inmat
	error <- grep("\\@\\@", rownames(inmat)) # get a vector of rows with double at signs
	if(length(error)>0){
		print(paste("ATENCION!!! Se encontraron mutaciones con dobles asteriscos** ... ignorando:" ))
		print(paste(" Item", rownames(inmat)[error], "en línea: ", error))
		outmat <- inmat[-error,]
	}
	return(outmat)
}

 ### MAIN ###
gene_list <- unique(df[,get_column("gene",df)]) # Print the gene list
print("Lista de genes encontrados:")
print(gene_list)

# Special case: deal with ORF8
# Historical data has double Q27 mutation data, thus we need to sum them, discard them or skip this alltogether I first decided to remove those with double **
df <- remove_double_at_signs(df) # Comment if required

# Analyze each item separately and create outputs
analyze_gene(df, "S")
analyze_gene(df, "N")
analyze_gene(df, "M")
analyze_gene(df, "ORF1a")
analyze_gene(df, "ORF1b")
analyze_gene(df, "ORF8")
analyze_gene(df, "ORF3a")

print("--- End of execution ---")
