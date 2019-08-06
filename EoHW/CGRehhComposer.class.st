Class {
	#name : #CGRehhComposer,
	#superclass : #CGRehh,
	#category : #'EoHW-Core'
}

{ #category : #'instance creation' }
CGRehhComposer class >> exampleGenerate04 [

	self new
		outputRehhScript: 'rehh_workflow-all-chr.R';
		outputShellScript: 'run-rehh-all-chr';
		outputRehhManhattanPlotScript: 'rehh_manhattan-wgs.R';		
		cwd: 'D:\\Huellas_de_Seleccion\\CBBA+CAl\\Chromosomes';
		chrRangeString: '1 29';
		ancestralPopNr: 1;		
		
		popCount: 2;
		popNames: {  
			'Angus_geno-maf_exclude_recode' . 
			'Brangus_geno-maf_exclude_recode' };
		
		minPercGenoHap: 95;
		minPercGenoSnp: 50;
		ies2rsbMethod: 'bilateral';

		writeHistogramFunction: false;

		writeBifurcationDiagram: false;
		mrkFoc: 2429;
		allFoc: 2;
		nmrkL: 10;
		nmrkR: 10;
		bifurcationDiagramTitle: 'Bifurcation diagram (XX SNP on BTA1): Derived Allele';

		generate.
		
	self inform: 'Script generated'.
]

{ #category : #'code generation' }
CGRehhComposer >> buildCompositionFormula: stream [

	stream
		nextPutAll:('############################
# Apply composition formula
############################

compPop1 <- 5/8
compPop2 <- 3/8

# Iterate simultaneously both Angus and Brahman frequency tables
f0_snps <- apply(data.frame(al_a1, al_b1), 1, function(row) (compPop1 * row[1]) + (compPop2 * row[2]))
f1_snps <- apply(data.frame(al_a2, al_b2), 1, function(row) (compPop1 * row[1]) + (compPop2 * row[2]))

maxrow <- function(row) { 
  if (max(row[1],row[2]) == row[1])
    c(1,0)
  else if (max(row[2],row[1]) == row[2])
        c(0,1)
      else
        c(1,0)
}
# Iterate simultaneously both frequency tables and assign 1 to maximum of each row
map <- apply(
    data.frame(f0_snps, f1_snps), 
    1,
    maxrow)
map <- t(map)

# Iterate simultaneously both frequency tables and assign 1 to maximum of each row
# map <- apply(data.frame(f0_snps, f1_snps), 1, function(row) ifelse(max(row[1],row[2]),1,0))').
]

{ #category : #'code generation' }
CGRehhComposer >> buildRehhFiles: stream [

	1 to: self popCount do: [ : popNr |
		stream
			nextPutAll: (self freqCountRCode: popNr) ].
	self buildCompositionFormula: stream.
	1 to: self popCount do: [ : popNr |
		stream 
			nextPutAll: (self mapGenRCode: popNr) ].
]

{ #category : #'code generation' }
CGRehhComposer >> freqCountRCode: popNr [

	^ '############################################
# Build MAP
############################################

########################
# Count Frequencies 1
########################

# Count frequencies of haps1 (0,1) by row and convert to matrix
count_frq <- apply(haps1[,-(1:5)], 1, function(x) \{ table(factor(x, levels=0:1)) \})
# Column 1 = 0
# Column 2 = 1
tr_freq <- t(as.data.frame(count_frq))

# Sum frequencies for each allele by row (SNP)
al_a1 <- tr_freq[,c(T,F)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])
al_b1 <- tr_freq[,c(F,T)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])

########################
# Count Frequencies 2
########################

# Count frequencies of haps2 (0,1) by row and convert to matrix
count_frq <- apply(haps2[,-(1:5)], 1, function(x) \{ table(factor(x, levels=0:1)) \})
# Column 1 = 0
# Column 2 = 1
tr_freq <- t(as.data.frame(count_frq))

# Sum frequencies for each allele by row (SNP)
al_a2 <- tr_freq[,c(T,F)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])
al_b2 <- tr_freq[,c(F,T)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])'  format: { popNr }

]

{ #category : #'code generation' }
CGRehhComposer >> mapGenRCode: popNr [

	^ '######################
# Write MAP {1}
######################

#View(map)
#str(map)

maps{1} <- data.frame(haps{1}[,2], haps{1}[,1], haps{1}[,3], anc = map[,1] , der = map[,2])

# maps{1} <- haps{1}[c(2,1,3)
finalMap{1} <- maps{1}

mapFileName{1} <- paste0("chr-", args[1] ,"rehh_map_pop-{1}.txt")

# write to file
write.table(
    finalMap{1}, 
    file = mapFileName{1}, 
    quote = FALSE, 
    col.names = FALSE, 
    row.names = FALSE)

############################################
# HAP Formatting {1}
############################################

# Remove header row from samples file
sample{1} <- sample{1}[-1,]

# Convert haplotype data frame into matrix
haps{1} <- haps{1}[,-(1:5)] %>% as.matrix

# Assign column/row names as positions/individual IDs
row.names(haps{1}) <- haps{1}[,3]
colnames(haps{1}) <- rep(sample{1}$ID_1, each = 2)

# Add 1 to each cell
haps{1} <- haps{1} + 1

# Transpose
haps{1} <- t(haps{1})

# Annotate rows
haplo.names <- rep(sample{1}$ID_2, each = 2)
first.haplos <- seq(from = 1, to = length(haplo.names), by = 2)
second.haplos <- seq(from = 2, to = length(haplo.names), by = 2)

haplo.names[first.haplos] <- paste0(haplo.names[first.haplos],"_1")
haplo.names[second.haplos] <- paste0(haplo.names[second.haplos],"_2")

haps{1} <- data.frame(haps{1}) %>%
  mutate(ind = haplo.names) %>%
  select(ind, everything())
# write to file (necessary for rehh for some reason)

tableHap{1} <- paste0("chr-", args[1] , "-rehh_hap_pop-{1}.txt") 

write.table(haps{1}, file = tableHap{1}, quote = FALSE, col.names = FALSE, row.names = FALSE)

' format: { popNr }
]
