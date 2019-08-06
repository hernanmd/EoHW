Class {
	#name : #CGRehhFrqCounter,
	#superclass : #CGRehh,
	#category : #'EoHW-Core'
}

{ #category : #'code generation' }
CGRehhFrqCounter >> buildRehhFiles: stream [

	1 to: self popCount do: [ : popNr |
		stream
			nextPutAll: ('############################################
# Build MAP {1}
############################################

########################
# Count Frequencies {1}
########################

# Count frequencies of haps{1} (0,1) by row and convert to matrix
#View(haps{1}[,-(1:5)])
count_frq <- apply(haps{1}[,-(1:5)], 1, function(x) \{ table(factor(x, levels=0:1)) \})
# Column 1 = 0
# Column 2 = 1
tr_freq <- t(as.data.frame(count_frq))

# Sum frequencies for each allele by row (SNP)
al_a <- tr_freq[,c(T,F)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])
al_b <- tr_freq[,c(F,T)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])
#View(data.frame(al_a,al_b))
# Assign 1s to major frequency for each row
freqsMatrix <- t(apply(data.frame(al_a, al_b), 1, function(row) 1 * (row == max(row))))
# Check number of snps (rows)
#dim(freqsMatrix)

# all 1s become 2s and 0s become 1s
freqsMatrix[freqsMatrix == 1] = 1
freqsMatrix[freqsMatrix == 0] = 2
#dim(freqsMatrix)
#View(freqsMatrix)

######################
# Write MAP {1}
######################

# Major Allele in Control population will be "ancestral"
maps{1} <- haps{1}[c(2,1,3)];
finalMap{1} <- as.matrix(cbind(maps{1}, freqsMatrix))
# haps{1}[which(FR < 0.5),]
#dim(finalMap)

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

' format: { popNr }) ].
]
