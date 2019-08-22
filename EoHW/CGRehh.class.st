Class {
	#name : #CGRehh,
	#superclass : #Object,
	#instVars : [
		'cwd',
		'popCount',
		'outputDirectory',
		'outputRehhScript',
		'outputShellScript',
		'chrRangeString',
		'popNames',
		'ies2rsbMethod',
		'minPercGenoHap',
		'minPercGenoSnp',
		'writeBifurcationDiagram',
		'writeHistogramFunction',
		'mrkFoc',
		'allFoc',
		'nmrkL',
		'nmrkR',
		'bifurcationDiagramTitle',
		'outputRehhManhattanPlotScript',
		'ancestralPopNr'
	],
	#category : #'EoHW-Core'
}

{ #category : #accessing }
CGRehh >> allFoc [
	^ allFoc
]

{ #category : #accessing }
CGRehh >> allFoc: anObject [
	allFoc := anObject
]

{ #category : #accessing }
CGRehh >> ancestralPopNr [
	^ ancestralPopNr
]

{ #category : #accessing }
CGRehh >> ancestralPopNr: anObject [
	ancestralPopNr := anObject
]

{ #category : #accessing }
CGRehh >> bifurcationDiagramTitle [
	^ bifurcationDiagramTitle
]

{ #category : #accessing }
CGRehh >> bifurcationDiagramTitle: anObject [
	bifurcationDiagramTitle := anObject
]

{ #category : #accessing }
CGRehh >> buildRehhFiles: arg1 [ 
	^ self subclassResponsibility
]

{ #category : #'code generation' }
CGRehh >> chrRangeR [

	^ self chrRangeString copyReplaceAll: ' ' with: ':'
]

{ #category : #accessing }
CGRehh >> chrRangeString [
	^ chrRangeString
		ifNil: [ chrRangeString := '1 29' ]
]

{ #category : #accessing }
CGRehh >> chrRangeString: anObject [
	chrRangeString := anObject
]

{ #category : #'code generation' }
CGRehh >> cwd [

	^ cwd
		ifNil: [ cwd := 'getwd()' ]
]

{ #category : #accessing }
CGRehh >> cwd: anObject [
	cwd := anObject
]

{ #category : #'rehh functions' }
CGRehh >> eHHScan: popNr [
	^ '

######################
# iHS {1}
######################

ihsFileName{1} <- paste0("chr-", args[1] , "-rehh_ihs_pop{1}.txt")  
ihs{1} <- ihh2ihs(ehhScan{1}, freqbin=0.5)
write.table(ihs{1}, file = ihsFileName{1})
write.csv(ihs{1}, file = "rehh_table_ihs_pop-{1}.csv")

ihsFileName{1}Png <- paste0("chr-", args[1] , "-rehh_ihs_pop-{1}.png") 
ihsTitleName{1} <- paste0("iHS Population {1} chr " , args[1])
png(filename = ihsFileName{1}Png)
ihsplot_mod(ihs{1}, plot.pval = TRUE, ylim.scan=2, main = ihsTitleName{1})
dev.off()

' format: { popNr }
]

{ #category : #'code generation' }
CGRehh >> generate [

	self 
		generateMainRscript;
		generateManhattanPlotRscript;
		generateRunAllShellScript.
	self inform: 'rehh script generated'.
]

{ #category : #'code generation' }
CGRehh >> generateMainRscript [
	" Private - Write files to disk "

	self outputRehhPatches asFileReference writeStreamDo: [ : stream |
		self writeRehhPatches: stream ].
	self outputRehhScript asFileReference writeStreamDo: [ : stream |
		self writeHeaderComment: stream.
		self writeLoadLibs: stream.
		self loadInputFiles: stream.
		self buildRehhFiles: stream.
		self writeRehhFunctions: stream.
		self writeHistogramFunctionTo: stream.
		self writeBifurcationDiagramFunction: stream.
		self writeRsbFunction: stream ].

]

{ #category : #'code generation' }
CGRehh >> generateManhattanPlotRscript [

	self outputRehhManhattanPlotScript asFileReference writeStreamDo: [ : stream |
		self writeManhattanPlotFunction: stream ].

]

{ #category : #'code generation' }
CGRehh >> generateRunAllShellScript [

	self outputShellScript asFileReference writeStreamDo: [ : stream |
		self writeShellScript: stream ].

]

{ #category : #accessing }
CGRehh >> ies2rsbMethod [
	^ ies2rsbMethod
]

{ #category : #accessing }
CGRehh >> ies2rsbMethod: anObject [
	ies2rsbMethod := anObject
]

{ #category : #'code generation' }
CGRehh >> installLibsString [

	^ '
list.of.packages <- c("dplyr", "gplots", "ggplot2", "rehh")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.us.r-project.org")'

]

{ #category : #'code generation' }
CGRehh >> loadInputFiles: stream [

	1 to: self popCount do: [ : popNr |  
		stream
			nextPutAll: ('haps{1} <- read.table(args[{2}], stringsAsFactors = FALSE)
sample{1} <- read.table(args[{3}], header = TRUE, stringsAsFactors = FALSE)
' format: { popNr . (2 * popNr) . (2 * popNr + 1)  }) ].

]

{ #category : #'code generation' }
CGRehh >> loadLibsString [

	^ 'library(dplyr)
library(rehh)
library(ggplot2)
source("rehh_patches.R")
'
]

{ #category : #accessing }
CGRehh >> minPercGenoHap [
	^ minPercGenoHap
]

{ #category : #accessing }
CGRehh >> minPercGenoHap: anObject [
	minPercGenoHap := anObject
]

{ #category : #accessing }
CGRehh >> minPercGenoSnp [
	^ minPercGenoSnp
]

{ #category : #accessing }
CGRehh >> minPercGenoSnp: anObject [
	minPercGenoSnp := anObject
]

{ #category : #accessing }
CGRehh >> mrkFoc [
	^ mrkFoc
]

{ #category : #accessing }
CGRehh >> mrkFoc: anObject [
	mrkFoc := anObject
]

{ #category : #'code generation' }
CGRehh >> multiChrMainShellScript [

	^ '#!/bin/sh

# set -x
#############################
# Run all chromosomes
#############################

mac2unix {3}

pop1="{1}"
pop2="{2}"
rscript="{3}"

chrrange=$(seq {4})

for i in $chrrange; do
	Rscript --verbose "$rscript" $i \\
		"$pop1$i".PHASED.haps \\
		"$pop1$i".PHASED.sample \\
		"$pop2$i".PHASED.haps \\
		"$pop2$i".PHASED.sample
done ' format: { 
			self popNames first . 
			self popNames second . 
			self outputRehhScript . 
			self chrRangeString }
]

{ #category : #'code generation' }
CGRehh >> multiChrManhattanPlotScript [

	^ String streamContents: [ : stream | 
			stream 
				cr;
				nextPutAll: 'mac2unix ';
				nextPutAll: self outputRehhManhattanPlotScript;
				cr;
				nextPutAll: 'Rscript --verbose ';
				nextPutAll: self outputRehhManhattanPlotScript ]

	

]

{ #category : #accessing }
CGRehh >> nmrkL [
	^ nmrkL
]

{ #category : #accessing }
CGRehh >> nmrkL: anObject [
	nmrkL := anObject
]

{ #category : #accessing }
CGRehh >> nmrkR [
	^ nmrkR
]

{ #category : #accessing }
CGRehh >> nmrkR: anObject [
	nmrkR := anObject
]

{ #category : #accessing }
CGRehh >> outputDirectory [
	^ outputDirectory
]

{ #category : #accessing }
CGRehh >> outputDirectory: anObject [
	outputDirectory := anObject
]

{ #category : #accessing }
CGRehh >> outputRehhManhattanPlotScript [
	^ outputRehhManhattanPlotScript
]

{ #category : #accessing }
CGRehh >> outputRehhManhattanPlotScript: anObject [
	outputRehhManhattanPlotScript := anObject
]

{ #category : #'code generation' }
CGRehh >> outputRehhPatches [
	" Answer a <String> referencing the file name for Rehh patches "

	^ 'rehh_patches.R'
]

{ #category : #accessing }
CGRehh >> outputRehhScript [
	^ outputRehhScript
]

{ #category : #accessing }
CGRehh >> outputRehhScript: anObject [
	outputRehhScript := anObject
]

{ #category : #accessing }
CGRehh >> outputShellScript [
	^ outputShellScript
]

{ #category : #accessing }
CGRehh >> outputShellScript: anObject [
	outputShellScript := anObject
]

{ #category : #'code generation' }
CGRehh >> popCount [

	^ popCount
		ifNil: [ popCount := 2 ]
]

{ #category : #accessing }
CGRehh >> popCount: anObject [
	popCount := anObject
]

{ #category : #accessing }
CGRehh >> popNames [
	^ popNames
]

{ #category : #accessing }
CGRehh >> popNames: anObject [
	popNames := anObject
]

{ #category : #'rehh functions' }
CGRehh >> scanHHPop: popNr [

	^ '
#######################
# Scan_hh() funtion {1}
#######################

# Haplotype are recoded (if recode.allele option is activated) according to the ancestral and derived allele 
# definition available in the map file (fourth and fifth columns) as :
# 0=missing data, 
# 1=ancestral allele, 
# 2=derived allele. 

hap{1} <- data2haplohh(
        hap_file = tableHap{1}, 
        map_file = mapFileName{4}, 
    		# haplotype.in.columns = TRUE , 
  			# recode.allele = TRUE, 
			# min_maf = 0.1,
      		# min_perc_geno.hap = {2},
        	# min_perc_geno.snp = {3},
        chr.name = args[1])
ehhScan{1} <- scan_hh(hap{1})
scanFileName{1} <- paste0("chr-", args[1] , "-rehh_scanhh_pop-{1}.txt") 
write.table(
    ehhScan{1},
    file = scanFileName{1})

' format: { popNr . self minPercGenoHap . self minPercGenoSnp . self ancestralPopNr }.
]

{ #category : #'code generation' }
CGRehh >> versionString [

	^ '1.0'
]

{ #category : #accessing }
CGRehh >> writeBifurcationDiagram [
	^ writeBifurcationDiagram
]

{ #category : #accessing }
CGRehh >> writeBifurcationDiagram: anObject [
	writeBifurcationDiagram := anObject
]

{ #category : #'code generation' }
CGRehh >> writeBifurcationDiagramFunction: stream [

	self writeBifurcationDiagram
		ifTrue: [ 
			1 to: self popCount do: [ : popNr |
				stream
					nextPutAll: ('########################
# Bifurcation.diagram {1}
########################

bifFileName{1}Png <- paste0("chr-", args[1] , "-rehh_bifurcation_pop-{1}.png") 
png(filename = bifFileName{1}Png)
bifurcation.diagram(test.brang,mrk_foc={2},all_foc={3},nmrk_l={4},nmrk_r={5},main="{6}")

dev.off
' format: { popNr . self mrkFoc . self allFoc . self nmrkL . self nmrkR . self bifurcationDiagramTitle }) ] ].
]

{ #category : #'code generation' }
CGRehh >> writeHeaderComment: stream [

	^ stream nextPutAll: ('# rehh Workflow Script version {1}

' format: { self versionString })
]

{ #category : #accessing }
CGRehh >> writeHistogramFunction [
	^ writeHistogramFunction
]

{ #category : #accessing }
CGRehh >> writeHistogramFunction: anObject [
	writeHistogramFunction := anObject
]

{ #category : #'code generation' }
CGRehh >> writeHistogramFunctionTo: stream [

	self writeHistogramFunction
		ifTrue: [ 
			stream nextPutAll: '
plotRehhHistogram <- function(binWidthParam) \{ 
	histogramFileName{1}Png <- paste0("chr-", args[1] , "-bin-", binWidth , "-rehh_histogram_pop{1}.png")  
	png(filename = histogramFileName{1}Png)
	ehhScan{1} %>% ggplot() +
		geom_histogram(aes(x=freq_A), binwidth = binWidthParam)
	dev.off() \}

sapply(seq(from = 10, to = 900, by = 10), plotRehhHistogram)

' ]
]

{ #category : #'code generation' }
CGRehh >> writeLoadLibs: stream [

	stream 
		nextPutAll: ('setwd("{1}"' format: { self cwd });
		nextPutAll: self installLibsString;
		nextPutAll: self loadLibsString;
		crlf;
		nextPutAll: 'args<-commandArgs(TRUE)'

]

{ #category : #'code generation' }
CGRehh >> writeManhattanPlotFunction: stream [

	stream 
		nextPutAll: self loadLibsString.
	1 to: self popCount do: [ : popNr |
		stream
			nextPutAll: ('
for(i in {2}) \{
		hapFile <- paste0("chr-", i , "-rehh_hap_pop-{1}.txt")
		mapFile <- paste0("chr-", i ,"-rehh_map_pop-{1}.txt")
		data <- data2haplohh(
				hap_file = hapFile, 
				mapFile, 
				chr.name=i
				# haplotype.in.columns = TRUE, 
				# recode.allele = TRUE, 
				# min_maf = 0.1,
				# min_perc_geno.hap = {3},
				# min_perc_geno.snp = {4}
				)
		res <- scan_hh(data)
		if (i == 1) 
			wg.res <- res
		else 
			wg.res <- rbind(wg.res, res) 
		\}
	wg.ihs <- ihh2ihs(wg.res, freqbin=0.5)
	
ihsFileNamePng <- paste0("wgs-rehh_ihs_pop-{1}.png") 
ihsTitleName <- paste0("WGS iHS Population {1}")
png(filename = ihsFileNamePng)
ihsplot_mod(wg.ihs, plot.pval = TRUE, ylim.scan=2, main = ihsTitleName)
dev.off()	
	
' format: { popNr . self chrRangeR . self minPercGenoHap . self minPercGenoSnp }) ].

]

{ #category : #'code generation' }
CGRehh >> writeRehhFunctions: stream [

	1 to: self popCount do: [ : popNr |
		stream
			nextPutAll: (self scanHHPop: popNr);
			nextPutAll: (self eHHScan: popNr) ].
]

{ #category : #'code generation' }
CGRehh >> writeRehhPatches: stream [

	stream 
		nextPutAll: '
ihsplot_mod<-function (ihsdata, plot.pval = TRUE, ylim.scan = 2, pch = 16,  cex = 0.5, cex.lab = 1.25, main = NA, cex.main = 1.5, cex.axis = 1)
{
  if (is.null(ihsdata$iHS) | is.null(ihsdata$frequency.class)) {
    stop("The ihsdata list is invalid and was not generated by the ihh2ihs() function...")
  }
  data = ihsdata$iHS
  lst_chrm <- unique(data$CHR)
  nbr_chrm <- length(lst_chrm)
  cum <- rep(0, (nbr_chrm + 1))
  if (nbr_chrm > 1) {
    for (i in 1:nbr_chrm) {
      cum[i + 1] <- tail(data$POSITION[data$CHR == lst_chrm[i]],
                         n = 1)
    }
  }
  cum <- cumsum(cum)
  pos <- rep(0, length(data$CHR))
  pos_labels <- rep(0, nbr_chrm)
  for (i in 1:nbr_chrm) {
    pos[data$CHR == lst_chrm[i]] <- data$POSITION[data$CHR ==
                                                    lst_chrm[i]] + cum[which(lst_chrm == lst_chrm[i])]
    pos_labels[i] <- (cum[i] + cum[i + 1])/2
  }
  col_chr <- data$CHR
  if (nbr_chrm > 1) {
    #     dev.new()
    #     plot.new()
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(pos, data$iHS, pch = pch, cex = cex, las = 1, col = col_chr,
         xaxt = "n", xlab = "Chromosome", ylab = expression(italic(iHS)),
         cex.lab = cex.lab, main = main, cex.main = cex.main,
         cex.axis = cex.axis)
    axis(1, at = pos_labels, labels = lst_chrm, las = 1,
         cex.lab = cex.lab, cex.axis = cex.axis)
    abline(h = c(-ylim.scan, ylim.scan), lty = 2)
    if (plot.pval) {
      #         dev.new()
      #         plot.new()
      par(mar = c(5, 5, 4, 2) + 0.1)
      plot(pos, data$"-log10(p-value)", pch = pch, cex = cex,
           las = 1, col = col_chr, xaxt = "n", xlab = "Chromosome",
           ylab = expression("-" * log[10] * "[1" ~ "-" ~
                               "2" ~ "|" ~ Phi[scriptstyle(italic(iHS))] ~
                               "-" ~ 0.5 * "|]"), cex.lab = cex.lab, main = main,
           cex.main = cex.main, cex.axis = cex.axis)
      axis(1, at = pos_labels, labels = lst_chrm, las = 1,
           cex.lab = cex.lab, cex.axis = cex.axis)
      abline(h = c(-ylim.scan, ylim.scan), lty = 2)
    }
  }
  else {
    if (max(pos) < 1000) {
      scale <- 1
      unit = "(bp)"
    }
    else if (max(pos) < 1e+06) {
      scale <- 1000
      unit = "(kb)"
    }
    else if (max(pos) < 1e+09) {
      scale <- 1e+06
      unit = "(Mb)"
    }
    else {
      scale <- 1e+09
      unit = "(Gb)"
    }
    #      dev.new()
    #       plot.new()
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(pos/scale, data$iHS, pch = pch, col = col_chr, xlab = paste("Position",
                                                                     unit), ylab = expression(italic(iHS)), cex.lab = cex.lab,
         main = main, cex.main = cex.main, cex.axis = cex.axis)
    abline(h = c(-ylim.scan, ylim.scan), lty = 2)
    if (plot.pval) {
      #         dev.new()
      #        plot.new()
      par(mar = c(5, 5, 4, 2) + 0.1)
      plot(pos/scale, data$"-log10(p-value)", pch = pch,
           col = col_chr, xlab = paste("Position", unit),
           ylab = expression("-" * log[10] * "[1" ~ "-" ~
                               "2" ~ "|" ~ Phi[scriptstyle(italic(iHS))] ~
                               "-" ~ 0.5 * "|]"), cex.lab = cex.lab, main = main,
           cex.main = cex.main, cex.axis = cex.axis)
      abline(h = c(-ylim.scan, ylim.scan), lty = 2)
    }
  }
  par(mar = c(5, 4, 4, 2) + 0.1)
}

distribplot_mod<-function (data, lty = 1, lwd = 1.5, col = c("blue", "red"), main = "Genome-wide distribution", 
          xlab = "", cex.main = 1.5, cex.lab = 1.25, qqplot = TRUE) 
{
#  dev.new()
#  plot.new()
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(density(data, na.rm = TRUE), main = main, xlab = xlab, 
       col = col[1], lty = lty, lwd = lwd, cex.main = cex.main, 
       cex.lab = cex.lab)
  curve(dnorm, col = col[2], add = TRUE)
  legend("topright", c("Observed", "Gaussian"), bty = "n", 
         col = col, lty = lty, lwd = lwd)
  if (qqplot) {
    #dev.new()
    #plot.new()
    par(mar = c(5, 5, 4, 2) + 0.1)
    qqnorm(data[!is.na(data)], cex.main = cex.main, cex.lab = cex.lab, 
           pch = 16, cex = 0.75)
    abline(a = 0, b = 1, lty = 2)
  }
  par(mar = c(5, 4, 4, 2) + 0.1)
}

rsbplot_mod<-function (data, plot.pval = TRUE, ylim.scan = 2, pch = 16, cex = 0.5, 
          cex.lab = 1.25, main = NA, cex.main = 1.5, cex.axis = 1) 
{
  if (!(grepl("Rsb", colnames(data)[3]))) {
    stop("Unrecognized column name for Rsb: the data is not a matrix generated by the ies2rsb() function...")
  }
  if (colnames(data)[4] == "-log10(p-value) [bilateral]") {
    yleg = expression("-" * log[10] * "[1" ~ "-" ~ "2" ~ 
                        "|" ~ Phi[scriptstyle(italic(Rsb))] ~ "-" ~ 0.5 * 
                        "|]")
  }
  else if (colnames(data)[4] == "-log10(p-value) [unilateral]") {
    yleg = expression("-" * log[10] * "[1" ~ "-" ~ Phi[scriptstyle(italic(Rsb))] * 
                        "|]")
  }
  else {
    warning("Unrecognized column name for p-values: plot.pval has been turned off")
    plot.pval = FALSE
  }
  lst_chrm <- unique(data$CHR)
  nbr_chrm <- length(lst_chrm)
  cum <- rep(0, (nbr_chrm + 1))
  if (nbr_chrm > 1) {
    for (i in 1:nbr_chrm) {
      cum[i + 1] <- tail(data$POSITION[data$CHR == lst_chrm[i]], 
                         n = 1)
    }
  }
  cum <- cumsum(cum)
  pos <- rep(0, length(data$CHR))
  pos_labels <- rep(0, nbr_chrm)
  for (i in 1:nbr_chrm) {
    pos[data$CHR == lst_chrm[i]] <- data$POSITION[data$CHR == 
                                                    lst_chrm[i]] + cum[which(lst_chrm == lst_chrm[i])]
    pos_labels[i] <- (cum[i] + cum[i + 1])/2
  }
  col_chr <- data$CHR
  if (nbr_chrm > 1) {
    #dev.new()
    #plot.new()
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(pos, data[, 3], pch = pch, cex = cex, las = 1, col = col_chr, 
         xaxt = "n", xlab = "Chromosome", ylab = expression(italic(Rsb)), 
         cex.lab = cex.lab, main = main, cex.main = cex.main, 
         cex.axis = cex.axis)
    axis(1, at = pos_labels, labels = lst_chrm, las = 1, 
         cex.lab = cex.lab, cex.axis = cex.axis)
    abline(h = c(-ylim.scan, ylim.scan), lty = 2)
    if (plot.pval) {
      # dev.new()
      # plot.new()
      par(mar = c(5, 5, 4, 2) + 0.1)
      plot(pos, data[, 4], pch = pch, cex = cex, las = 1, 
           col = col_chr, xaxt = "n", xlab = "Chromosome", 
           ylab = expression(-log[10] ~ "(" * italic(p) * 
                               "-value)"), cex.lab = cex.lab, main = main, 
           cex.main = cex.main, cex.axis = cex.axis)
      axis(1, at = pos_labels, labels = lst_chrm, las = 1, 
           cex.lab = cex.lab, cex.axis = cex.axis)
      abline(h = c(-ylim.scan, ylim.scan), lty = 2)
    }
  }
  else {
    if (max(pos) < 1000) {
      scale <- 1
      unit = "(bp)"
    }
    else if (max(pos) < 1e+06) {
      scale <- 1000
      unit = "(kb)"
    }
    else if (max(pos) < 1e+09) {
      scale <- 1e+06
      unit = "(Mb)"
    }
    else {
      scale <- 1e+09
      unit = "(Gb)"
    }
    # dev.new()
    # plot.new()
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(pos/scale, data[, 3], pch = pch, col = col_chr, 
         xlab = paste("Position", unit), ylab = expression(italic(Rsb)), 
         cex.lab = cex.lab, main = main, cex.main = cex.main, 
         cex.axis = cex.axis)
    abline(h = c(-ylim.scan, ylim.scan), lty = 2)
    if (plot.pval) {
      # dev.new()
      # plot.new()
      par(mar = c(5, 5, 4, 2) + 0.1)
      plot(pos/scale, data[, 4], pch = pch, col = col_chr, 
           xlab = paste("Position", unit), ylab = expression(-log[10] ~ 
                                                               "(" * italic(p) * "-value)"), cex.lab = cex.lab, 
           main = main, cex.main = cex.main, cex.axis = cex.axis)
      abline(h = c(-ylim.scan, ylim.scan), lty = 2)
    }
  }
  par(mar = c(5, 4, 4, 2) + 0.1)
}
'
]

{ #category : #'code generation' }
CGRehh >> writeRsbFunction: stream [

	stream 
		nextPutAll: ('
#######################
# Rsb
#######################

#dim(ehhScan1)
#dim(ehhScan2)

rsb <- ies2rsb(ehhScan1, ehhScan2, method="{2}")
rsbTableName <- paste0("chr-", args[1] , "-rehh_rsb-{2}.csv")
write.table(rsb, file = rsbTableName)
#View(rsb)

rsbPngName <- paste0("chr-", args[1] , "-rehh_rsb-{2}.png")
png(filename = rsbPngName)
rsbTitle <- paste0("Rsb chr " , args[1])
rsbplot_mod(rsb, main = rsbTitle)
dev.off()
' format: { self ies2rsbMethod . self ies2rsbMethod })

]

{ #category : #'code generation' }
CGRehh >> writeShellScript: stream [

	stream 
		nextPutAll: self multiChrMainShellScript;
		nextPutAll: self multiChrManhattanPlotScript

]
