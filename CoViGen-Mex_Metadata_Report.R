# UPDATE 2021-05-27: Added an extra column with the total items per table for the LaTeX file.
# First complete version (v1.0): 2021-05-26
# Started: 2021-05-24
# by Rodrigo García-López for Carlos Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
# IMPORTANT: Libraries readxl, xtable are required

# The script is intended to create summary tables and graphs from the CoViGen-Mex lineage tables (one lineage per sample)
# Output summary:
# 	1.- Files _linaje_single.tsv - These contain the lineage analyses from the current table (see NOTE 1)
# 	2.- Files _linaje_multi.tsv - These contain analyses from larger tables with past lineages and per-column percentages (see NOTE 1).
# 3.- Files _edad_vs_sexo, _linaje_vs_Estados, _fecha_vs_hosp - These contain crossed metadata tables and corresponding barplots in landscape. For _linaje_vs_Estados, see NOTE 1.

# For items 1 and 2, LaTeX should be installed in system for tools::texi2pdf to work and create pdf files with the tables.

# NOTE 1: VOI and VOC are optionally identified in these tables if an identification table is provided (with header; must have lineage and type)
######### VOI and VOC variants must be update manually in a sepatate table
##########

# All input tables must have a header
# MAIN INPUT: A regular table depicting the following information for each sample:
# Folio Interno	Edad (años)	Sexo	Estado	Municipio 	Fecha de toma	Tipo de muestra 	Tipo de paciente	Resultado	CT Gen 1	CT Gen 2	CT Gen 3	Embarazo	Semanas de gestación	Folio SINAVE	Folio SINOLAVE	Virus name	SampleIDProcess	Pangolin Clade	Clado Nextstrain	Mutaciones Aminoacidos Conocidas	Nuevas Mutaciones AminoÃ¡cido	Deleciones conocidas	Deleciones nuevas

# Test in R:
# library("readxl")
# library("xtable")
# setwd("/home/rod/Documents/01_Projects/SARS/Vigilancia")
# xlsx <- "Test2/output_test_Metadata_Pangolin.xlsx"
# prefix <- "Test2/output2"
# name_date <- "2021-05-19"
# past_linages <- "Test/5_fechas_linaje.tsv"
# variants <- "Test/VOC_VOI.tsv"

# Test in bash:
# Rscript CoViGen-Mex_Metadata_Report.R Test2/output_test_Metadata_Pangolin.xlsx Test2/output2 2021-05-19 Test/5_fechas_linaje.tsv Test/VOC_VOI.tsv

# Run as follows:
# Rscript CoViGen-Mex_Metadata_Report.R <input_xlsx> <prefix_output> <date> [OPTIONAL: past-lineages_table] [OPTIONAL: VOC_VOI_table]

 ### PRE-LOAD PARAMETERS AND DEFINE GLOBAL VARIABLES ###
library("readxl")
library("xtable")
# library("knitr")
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) { # at least, 3 arguments should be included: <prefix_output>  <set_name> <#_group_name_1> <#_group_name_n>
  stop("A minimum of 3 arguments are mandatory: Rscript CoViGen-Mex_Metadata_Report.R <input_xlsx> <prefix_output> <date> [OPTIONAL: past-lineages_table] [OPTIONAL: VOC_VOI_table]", call.=FALSE)
}
xlsx <- as.character(args[1]) # Get a string handle for the input xlsx file
prefix <- as.character(args[2]) # Get a string handle for the output (it can contain a complete path but requires a seed name for the output file)
name_date <- as.character(args[3]) # Get a string handle for the output (it can contain a complete path but requires a seed name for the output file)
past_linages <- as.character(args[4]) # OPTIONAL: Get a string handle for the table with past lineages' frequencies
variants <- as.character(args[5]) # OPTIONAL: Get a string handle for the table with current VOI and VOC
df <- as.data.frame(read_excel(xlsx, sheet = 1)) # Load the input table

 ### LOAD FUNCTIONS ###
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	return(string)
}
only_first_upcase <- function(string) { # Make only the first letter an uppercase character
  string <- strip(string)
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  return(string)
}
all_first_upcase <- function(string) { # Apply uppercase transform to each word separated by spaces
	string <- gsub("^ *|(?<= ) | *$", "", string, perl = TRUE) # Merge multiple spaces
	split <- strsplit(string, " ") # split each string item, create a list
	out <- lapply(split,function(x){paste(sapply(x, only_first_upcase),collapse=" ")}) # go through the complete list to apply the only_first_upcase function (must be defined as well) and concatenating the output
	out <- unlist(out)
	return(out)
}
get_column <- function(string, table){ # Gets a string and return a column number where header matches (not case-sensitive) or aborts if it not found
	query <- strip(string)
	col <- grep(query,strip(colnames(table)))
	if(length(col)!=1){stop(paste0("Abortando: No se encontro columna única en tabla ",substitute(table), " para query: ", substitute(query)), call.=FALSE)} # abort if none or > 1 matching item was found
	return(col)
}
new_color_scheme <- function(number=433){ # randomize default color names, default = max items in grDevices::colors
	return(sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],number))
}
colors_by_order <- function(vect,n=347){ # Define the color scheme: used a fixed seed for reproducibility
	if(file.exists("color_set.rds")){
		print("Using predefined color scheme from color_set.rds") # print the seed
		color=readRDS("color_set.rds") # if preset exists, use it
	}else{
		set.seed(n) # used a fixed seed for reproducibility
		print(paste("Seed used for color:",n)) # print the seed
		color=new_color_scheme(length(vect)) # Create a new color scheme for the total items
# 		print(cbind("RAW"=sort(vect,decreasing=T),color)) # Print which colors would match which items
		color=color[order(order(vect, decreasing=TRUE))];
# 		barplot(as.matrix(rev(vect)),col=rev(color),border=NA) # Test barplot with the set
		}
	return(color)
}
tag_variants <- function(ref_var, inmat) { # Labels VOI/VOC variants (requires user-provided list as argument)
		var <- read.table(ref_var, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, check.names=F, row.names=1) # load the list of reference variants
		out <- sapply(match(inmat[,1],rownames(var)),function(x){ifelse(is.na(x),NA,var[x,])}) # add a label to all those that contain them
		return(out)
}
get_lineage_vector <- function(inmat) {
	col_pango <- grep(strip("angolin"),strip(names(inmat))) # Get the column where the pangolin lineages are stored
	current_lineage <- as.data.frame(table(inmat[,col_pango])) # Create a df with the current lineage
	current_lineage <- cbind(current_lineage, current_lineage[,2]*100/sum(current_lineage[,2])) # append %
	colnames(current_lineage) <- c("Variant",name_date,paste(name_date,"%")) # rename columns
	current_lineage <- cbind(current_lineage, "VOI/VOC"=rep(NA,nrow(current_lineage))) # Initialize as empty vector
	if(!is.na(variants)){current_lineage[length(current_lineage)]=tag_variants(variants, current_lineage)} # Label VOI/VOC variants if the file is provided (argument)
	return(current_lineage) # and export
}
compare_past_lineages <- function(inmat, fullmat, vect) {
	lineages <- read.table(fullmat, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, check.names=F, row.names=1)
	temp <- as.data.frame(vect);temp <- temp[,1:2] # cast vector into dataframe
	m <- merge(lineages, temp, by.x = "row.names", by.y = "Variant", all = TRUE);rownames(m) <- m[,1]  # merge the tables and replace rownames with lineages
	var=rep(NA,nrow(m)) # Initialize empty variant variable as a vector the size of nrow(m)
	if(!is.na(variants)){var <- tag_variants(variants, m)} # Label VOI/VOC variants if the file is provided (argument)
	m <- m[,-1] # remove the var column to treat the rest as a matrix
	m[is.na(m)] <- 0 # Replace NAs with 0s if present
	totals <- colSums(m)
	m_p <- t(apply(m,1,function(x){x*100/totals})) # get the percentages
	colnames(m_p) <- paste(colnames(m_p),"%")
	out <- list(m,m_p,var); names(out) <- c("counts","percentage","var") # output both tables as a list of dfs
	return(out)
}
rename_grps <- function(vect){
	maskcat <- c("Otras","VOC/VOI no abundantes","Linajes abundantes") # Pre-define the substitutes
	names(maskcat)=c("0","1","2") # Label the 3 expected groups (as they are when vector is passed)
	mask <- as.factor(vect)
	levels(mask) <- maskcat[levels(mask)]
	return(mask)
}
collate_others_multi <- function(lineage_list){ # Order and collate lineages based on the lastest date
	out <- cbind(lineage_list$counts,lineage_list$percentage,"VOI/VOC"=lineage_list$var) # Append the lineages' absolute counts, relative counts (per-column) and the variants (if present)
	mask <- (!is.na(out[,length(out)]))*1 # Initialize a mask for preserving VOI/VOC
	tss <- out[,grep("%", names(out))] # Subset TSS (Total sum scaling) relative abundances
	mask <- mask+(tss[,ncol(tss)]>=1)*2# The last column has the newest lineages, we'ĺl use this to collate items with <1% or not VOI/VOC. Those with 1 are VOC/VOI but are not >1%, those with 2 have >1% in the last column but are not VOI/VOC and those with 3 are both things simultaneously
	mask[mask>1]=2 # Make all >1% the same group
	out <- cbind(out,"Nota"=mask) # append the new category
	out <- out[order(out$Nota,decreasing=TRUE),] # Order by importance (first the abundant ones, then the rare VOC/VOI, then the rest
	out2 <- rbind(out[out$Nota>=1,],"Otras"=c(colSums(out[out$Nota==0,-c(ncol(out),ncol(out)-1)]),NA,0)) # Create a collated version with those in low abundance and not VOC_VOI being clustered together
# 	maskcat <- c("Otras","VOC/VOI no abundantes","Linajes abundantes");names(maskcat)=c("0","1","2")
# 	mask <- as.factor(out$Nota);levels(mask) <- maskcat[levels(mask)];out$Nota <- mask # Recycle the mask variable, now as a factor, rename them and replace the note
# 	mask <- as.factor(out2$Nota);levels(mask) <- maskcat[levels(mask)];out2$Nota <- mask
	out$Nota <- rename_grps(out$Nota); out2$Nota <- rename_grps(out2$Nota) # rename accordingly
	rownames(out2)[nrow(out2)]=paste0(rownames(out2)[nrow(out2)], " (",nrow(out)-nrow(out2),")")
	result <- list(out,out2)
	names(result) <- c("Raw","Collated")
	return(result)
}
single_stack_barplot <- function(inmat,n=347){ # Gets matrix, outputs figure+table Uses one global variable: prefix
	vect <- inmat[,3];names(vect) <- inmat[,1] #
	color <- colors_by_order(vect,n)
	pdf(paste(prefix,"linaje_single-barplot.pdf",sep="_"))
		par(mfrow=c(1,2),las=2)
		v <- cumsum(rev(vect))-rev(vect)/2 # get the midpoint vertical positions of each bar (cummulative + 1/2 the size of each bar)
		i <- which(rev(vect)>3) # keep only those larger than 5%
		bp <- barplot(as.matrix(rev(vect)),col=rev(color),border=NA,yaxt='n', ylab="% de linaje")
		axis(1, las=1, at=bp, labels=,name_date)
		axis(2, at=seq(0,100,10), labels=paste0(seq(0,100,10),"%"))
		text(bp, v[i], labels = rev(names(vect))[i],col="black")
		plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
		legend("top", legend = names(vect), pch=15, pt.cex=1, cex=0.8, bty='n', col = color, title="Linajes") # create a simple legend
	dev.off()
	outmat <- inmat[,c(1,3,2,4)];outmat[,2] <- paste(round(outmat[,2],1),"%"); colnames(outmat) <- c("Linaje","% Total","# Muestras","Tipo")
	outLatexTable(outmat,paste("Linajes",name_date), paste(prefix,"linaje_single-table", sep="_"),"rlrrl",c(-1,0,nrow(outmat)),1) # Output matrix to pdf using latex (requires latex installed in sys)
# 		latex <- print.xtable(xtable(inmat), print.results = FALSE)
# 		tools::texi2pdf(latex, clean = TRUE)
# 		print(xtable(inmat), include.rownames = FALSE)
}
multi_stack_barplot <- function(inmat,n=347){ # Gets matrix, outputs figure+table Uses one global variable: prefix
	sub <- inmat[,grep("%",colnames(inmat))]
	N <- ncol(inmat)
	names(sub) <- gsub(" %","",names(sub))
	if(ncol(sub)>5){sub <- sub[,(ncol(sub)-4):ncol(sub)]} # If there are more than 5 columns, subset to last 5
	rev <- sub[rev(rownames(sub)),] # Flip the matrix upside down
	vect <- sub[,ncol(sub)];names(vect) <- rownames(sub) #
	color <- colors_by_order(vect,n); color[length(color)]="gray" # use the abd aware color assignation and replace the "other"'s color with gray
	pdf(paste(prefix,"linaje_multi-barplot.pdf",sep="_"))
		par(las=2)
		bp <- barplot(as.matrix(rev), col=rev(color), border=NA, xaxt='n', yaxt='n', ylab="% de linaje")
		axis(2, at=seq(0,100,10), labels=paste0(seq(0,100,10),"%"))
		h <- rep(bp,colSums((sub>=3)*1)) # Get a coordinate for x
		v <- apply(sub, 2L, cumsum)- sub / 2;v <- v[sub>=3] # Get a coordinate for y
		l <- unlist(apply(sub>=3,2,function(x){rownames(sub)[x]})) # Get the actual labels to be printed
		text(h,100-v,labels=l)
		text(bp-0.3, par("usr")[1]-7, labels=colnames(sub), srt=45, pos=1, xpd=TRUE)
		plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
		legend("top", legend = rownames(sub), pch=15, pt.cex=2, cex=2, bty='n', col = color, title="Linajes") # create a simple legend
	dev.off()
	outmat <- cbind(inmat[N],"Linaje"=rownames(inmat),inmat[,(N-2):(N-1)])
	hline <- unlist(lapply(sapply(unique(outmat[,1]),function(x){grep(x,outmat[,1])}),`[[`, 1)) # Get the first item of each group based on the last column in the matrix
	outmat[,1][(1:nrow(outmat))[-hline]] <- NA # empty repeated categories
	names(outmat) <- c("Nota","Linaje",name_date,"Tipo") # Rename columns
	outmat[,3] <- paste0(round(outmat[,3],2),"%") # Round % to 2 decimals
	outLatexTable(outmat,paste("Linajes",name_date), paste(prefix,"linaje_multi-table", sep="_"),"lp{3cm}lrl",c(-1,hline-1,nrow(outmat)),1)
	return(outmat)
}
bold <- function(x) {
	paste('{\\textbf{',x,'}}', sep ='')
}
gray <- function(x) {
	paste('{\\textcolor{gray}{',x,'}}', sep ='')
}
outLatexTable <- function(inmat,title,output,align,lines,width=1){
	colnames(inmat) <- gsub("%","\\%",names(inmat),fixed=TRUE)
	colnames(inmat) <- gsub("#","\\#",names(inmat),fixed=TRUE)
	table <- xtable(inmat)
	align(table) <- align
	latex <- print.xtable(table, print.results = FALSE,include.rownames=FALSE, sanitize.colnames.function=bold, sanitize.rownames.function=gray ,booktabs=T,hline.after=lines)
	writeLines( # Prepare a basic LaTex template:
	c(
		"\\documentclass[12pt]{standalone}",
		"\\usepackage{booktabs,colortbl,xcolor}",
		"\\usepackage[T1]{fontenc}",
		"\\usepackage[english]{babel}",
		"\\usepackage[scale=.9]{tgheros} % or helvet",
		"\\renewcommand{\\familydefault}{\\sfdefault}",
		"\\begin{document}",
		paste("\\minipage{",width,"\\textwidth}"),
		"\\bigskip",
		"\\centering",
		title,
		latex,
		"\\bigskip",
		"\\endminipage",
		"\\end{document}"
	),
	paste0(output,".tex") # Save the tex file
	)
	tools::texi2pdf(paste0(output,".tex"), clean = TRUE) # This requires Latex to be installed in system
	file.rename(from=paste0(basename(output),".pdf"),to=paste0(output,".pdf")) # move pdf to the outdir
}
paired_metadata_table <- function(string1,string2,inmat){ # This gets two string and a dataframe and creates a crossed-metadata table (two-way)
	meta1 <- inmat[,get_column(string1,inmat)] # These two will halt execution if the column name is not found in the table. They work regardless of the actual position in the table and whether or not they use uppercase or rare characters
	meta2 <- inmat[,get_column(string2,inmat)]
	pairtab <- xtabs(rep(1,nrow(inmat))~meta1+meta2, data=inmat) # Create a crossed table of lineages per location (State)
	return(pairtab)
}
locationVSlineages <- function(inmat,width,n=347){ # Create a 2way table with two metadata columns and create plots/tables (variants is a global variable)
	loc_lin <- paired_metadata_table("Estado","Pangolin",inmat) # Create Estado vs Linaje table
	Estados <- rownames(loc_lin) # Create a Estadosorary vector for holding the names of States
	Estados <- all_first_upcase(Estados) # Change all items to lower case (only first letter of each word to uc)
	Estados <- sub("Baja California Norte", "Baja California",Estados) # Fix some known issues (mostly accents)
	Estados <- sub("Distrito Federal", "Ciudad de México",Estados)
	Estados <- sub("Mexico", "Ciudad de México",Estados)
	Estados <- sub("Michoacan", "Michoacán",Estados)
	Estados <- sub("Leon", "León",Estados)
	Estados <- sub("Queretaro", "Querétaro",Estados)
	Estados <- sub("Potosi", "Potosí",Estados)
	Estados <- sub("Yucatan", "Yucatán",Estados)
	rownames(loc_lin) <- Estados # Now, replace with the updated names
	tss <- t(apply(loc_lin,1, function(x) x/colSums(loc_lin))) # Create a total sum scaling table (relative abd)
	bar_max_size <- log10(colSums(loc_lin)+1) # Calculate the total size each lineage has in log10, a +1 offset is used (this will define the size of each bar)
	plotdata <- t(t(tss)*bar_max_size) # calculate how much each lineage contributes to the total height of each bar
	set.seed(n);color=new_color_scheme()
	y <- c(0,1,2,5,10,20,50,100,200)
	tss <- tss*100 # Convert to percetages
	colnames(tss) <- paste(colnames(tss),"%")
	out <- list(loc_lin,tss)
	names(out)=c("abs","rel")
	pdf(paste(prefix,"linaje_vs_Estados-barplot.pdf", sep= "_"), width = 18, height = 8)
		bp <- barplot(plotdata,col=color,border=NA,xaxt='n',yaxt='n',ylab="Frecuencia (log) - Altura total de barras",cex.names=1.2) # Plot stacked bars: height represents each lineage's total in log10, and colors represent the % that comprise that total (not in log)
# 		axis(1,las=2,at=bp,labels=colnames(loc_lin),cex.axis=1.1)
		axis(1,las=2,at=bp,labels=NA,cex.axis=1.1) # Create only label ticks
		text(bp, 0, labels=paste0("     ",colnames(plotdata)), srt=270, xpd=TRUE,adj=0) # This fix was used to flip vertical labels
		axis(2,las=1,at=log10(y+1),labels=y, cex.axis=1.2) # y axis is log, labels are linear with a +1 offset
		vloc <- voc_voi <- "" # initialize empty holders
		if(!is.na(variants)){ # identify present VOI/VOC variants if the file is provided (argument)
			voc_voi=tag_variants(variants, data.frame(colnames(loc_lin))) # Create a list of VOI/VOC
			vloc=which(!is.na(voc_voi)) # get columns matching these
			text(bp[vloc], bar_max_size[vloc], labels="-",srt=90,adj=0) # Add a small tick (dash printed vertically and centered)
			text(bp[vloc], bar_max_size[vloc], labels=voc_voi[vloc],adj=c(0.5,-1),col="red",cex=1.2) # (add a VOC/VOI label if it so)
			mini <- cbind(loc_lin[,vloc],"Otras"=rowSums(loc_lin[,-vloc])) # subset the matrix to keep a mini matrix of only VOC/VOI whereas the others are summed
			mini <- mini[rowSums(mini[,-ncol(mini)])>0,] # keep only locations that have at least one VOC/VOI
			mini_tss <- apply(mini,2, function(x) x*100/rowSums(mini))
			tot <- rowSums(mini) # UPDATE 2021-05-27: Added an extra column with the total items per table for the LaTeX file. Changes are applied downstream.
			tab <- data.frame(cbind(rownames(mini_tss), tot, apply(mini_tss,2,function(x){paste0(sprintf("%.2f", x),"%")})))
# 			tab <- data.frame(cbind(rownames(mini_tss),mini_tss)) # append names column (this will recast items to "character")
# 			tab <- readr::type_convert(tab) # Reset type of fields for table, to avoid having "character"
			colnames(mini_tss) <- colnames(mini) <- c(paste(colnames(mini)[1:4],voc_voi[vloc]),"Otras")
			names(tab) <- c("Entidad Federativa","# Muestras",colnames(mini)) # Change names to include VOC/VOI info
			aln <- paste(c("l","l","r",rep("r",ncol(mini))),collapse="") # Create the alignment string
# 			aln <- paste(c("l","l","r",rep("p{1cm}",ncol(mini))),collapse="") # Same but reduce columns (but cannot align right)
			outLatexTable(tab,"Proporciones de variantes VOI/VOC por Estado", paste(prefix,"linaje_vs_Estados-table", sep="_"),aln,c(-1,0,nrow(tab)),width)
			colnames(mini_tss) <- paste(colnames(mini_tss),"%") # Append % for a more readable output
			out[[3]]=mini;out[[4]]=mini_tss # Append to output if present
			names(out)[3]="var_abs";names(out)[4]="var_rel"
		}
		plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
		legend("top", legend = Estados, pch=15, pt.cex=1, cex=1, bty='n', col = color, title="Estado") # create a simple legend
	dev.off()
	return(out)
}
check_dates <- function(inmat){ # Check for errors in the table
	dates <- inmat[,get_column("Fecha", inmat)] # get the date from patient data (as it may have errors)
	date_errors <- which(is.na(as.Date(dates)))
	if(length(date_errors)>0){print(paste("ADVERTENCIA!!! Error en Fecha de toma Tabla: PACIENTE| Línea:",date_errors))}
	return(date_errors)
}
age_vs_sex <- function(inmat,bin=10){ # this need a matrix with both sex and age cols and a bin defined by the user
	age_sex <- paired_metadata_table("Edad","Sexo",inmat) # Create Estado vs Linaje table
	ages <- as.numeric(rownames(age_sex)) # Get a vector with the variable ages (non-repeated)
	grp <- seq(bin,150,bin) # define the groups
	F <- sapply(grp,function(x){sum(age_sex[ages<(x)&ages>=(x-bin),1])}) # sum per bins and sex
	M <- sapply(grp,function(x){sum(age_sex[ages<(x)&ages>=(x-bin),2])})
	tab <- cbind(F,M)# create the output matrix
	rownames(tab) <- paste0("<",seq(bin,150,bin)); tab <- tab[rowSums(tab)>0,]# rename categories and retain only informative ones
	color <- c("purple","cornflowerblue") # Define two colors
 	pdf(paste(prefix,"edad_vs_sexo-barplot.pdf", sep= "_"), width = 18, height = 8) # make a landscape plot
		par(las=1) # define orientation
		barplot(t(tab),col=color,border=NA, space=0.5, cex.axis=1.5, cex.names=1.5, main="Edades", xlab="Años", ylab="Número de muestras", cex.lab=1.5, cex.main=1.5) # plot
		abline(h=0)
		legend("topleft",title="Genero", pch=15, col=color, legend=c("Femenino","Masculino"),bty = "n",cex=1.5) # legend
	dev.off()
	return(tab)
}
date_vs_hosp <- function(inmat){ #
	errors <- check_dates(inmat) # Define errors
	if(length(errors)>0){print("ADVERTENCIA!!! Se ignorarán las líneas con fechas faltantes.");inmat <- inmat[-errors,]} # Remove
	date_hosp <- paired_metadata_table("Fecha de toma","Tipo de paciente",inmat) # Create Fecha vs Hosp table
	date_hosp <- date_hosp[order(format(as.Date(rownames(date_hosp),"%d/%m/%Y"),"%Y/%m/%d")),] # Reformat dates (before: 30/06/2021, after: 2021/06/30) and use them to sort the dates correctly
	dates <- format(as.Date(rownames(date_hosp),"%d/%m/%Y"),"%b %d") # Reformat dates (before: 30/06/2021, after: Jun 30)
	dates <- gsub("Jan","Ene",dates) # Fix some language issues
	dates <- gsub("Apr","Abr",dates)
	dates <- gsub("Aug","Ago",dates)
	dates <- gsub("Dec","Dic",dates)
	rownames(date_hosp) <- dates # Replace in table with the new format
	colnames(date_hosp) <- all_first_upcase(colnames(date_hosp)) # Fix if uppercase
	date_hosp <- date_hosp[,c(get_column("amb",date_hosp),get_column("hosp",date_hosp),get_column("def",date_hosp))] # This may be optional but it is useful to have hopitalized patients before defunct ones
	color <- c("blue","orange","red") # Define two colors
 	pdf(paste(prefix,"fecha_vs_hosp-barplot.pdf", sep= "_"), width = 18, height = 8) # make a landscape plot
		par(las=2) # define orientation
		barplot(t(date_hosp),col=color,border=NA, space=0.5, cex.axis=1.5, cex.names=1.5, main="Fecha del muestreo", ylab="Número de muestras", cex.lab=1.5, cex.main=1.5) # plot
		abline(h=0)
		legend("topleft",title="Tipo de paciente", pch=15, col=color, legend=c("Ambulatorio","Hospitalizado","Defunción"),bty = "n",cex=1.5) # legend
	dev.off()
	return(date_hosp)
}

 ### MAIN ###
# LINEAGE ANALYSES
lin <- get_lineage_vector(df) # Get a vector summarizing the total variants per lineage
if(exists("lin")){ # if the past line was executed and the vector was created
	single_stack_barplot(lin,347) # It gets the current lineage data, and the color scheme (347 is an example passed to set.seed to create a fixed randomized color set)
	write.table(lin,paste(prefix,"linaje_single.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) # output the lineages table
	if(!is.na(past_linages)){ # If the a table with past lineages is provided
		lin_hist <- compare_past_lineages(df, past_linages, lin) # if the past lineages were provided (table as argument), append the new lineage vector to the existing table
		lineage_multi <- collate_others_multi(lin_hist)
		write.table(lineage_multi$Raw, paste(prefix,"linaje_multi-full.tsv", sep= "_"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) # output the lineages table
		write.table(lineage_multi$Collated, paste(prefix,"linaje_multi-collated.tsv", sep= "_"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) # output the lineages table
		out <- multi_stack_barplot(lineage_multi$Collated,347)
		write.table(out, paste(prefix,"linaje_single-table-collated.tsv", sep= "_"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) # output the summarized single lineage table
	}
}
# LOCATION VS LIUNEAGES
loclin <- locationVSlineages(df,1.5,347) # analyze locations and lineages, outputs two tables in a list with abs and rel abundance. It gets the matrix, the width for latex output (defaults to 1, so 1.5 is 50% more space) and the color scheme (347 is an example passed to set.seed to create a fixed randomized color set)
write.table(cbind(loclin$abs,loclin$rel),paste(prefix,"linaje_vs_Estados-abs_rel.tsv", sep="_"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) # output the lineages table
if(length(loclin)==4){write.table(cbind(loclin$var_abs,loclin$var_rel),paste(prefix,"linaje_vs_Estados-abs_rel-vocvoi.tsv", sep="_"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)} # output the lineages table

# AGE
agesex <- age_vs_sex(df,10) # This imports the input matrix and a bin size (starts at 0)
write.table(agesex,paste(prefix,"edad_vs_sexo-table.tsv", sep="_"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) # output the age vs sex table

# DATES
datehosp <- date_vs_hosp(df) # This inputs the matrix
write.table(datehosp,paste(prefix,"fecha_vs_hosp-table.tsv", sep="_"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) # output the date vs type of patient table

print("--- End of execution ---")
