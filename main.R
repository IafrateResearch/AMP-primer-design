
# get initial values
system("sed 's/: /=/g' config.txt | sed 's:{::' | sed 's:}::' | sed 's/, /; /g' > config.R")
source("config.R")

## create working dir for the panel
paneldir = paste(pjdir, panel, sep='/')
dir.create(paneldir, showWarnings = FALSE, recursive = TRUE)

system(paste("cp config.txt ", paneldir, sep=''))
system(paste("cp config.R ", paneldir, sep=''))

setwd(paneldir) 

## get gene list
if (assaytype == 'bed'){
	gene.list = read.table(genelist, sep='\t', header=F, fill=T)
		if (ncol(gene.list) <=3){
			gene.list$target = paste('target', 1:nrow(gene.list), sep='')
		}
	names(gene.list) = c('chrom', 'start', 'end', 'target')
	gene.list$target = gsub(' ', '.', gene.list$target)
} else {
	gene.list = read.table(genelist, sep='\t', header=T, fill=T)

	## remove blank space, heading 0s
	gene.list$NM = gsub(' ', '', gene.list$NM)
	gene.list$exon = tolower(gsub(' ', '', gene.list$exon))
	gene.list$exon = gsub('^0+', '', gene.list$exon)
}


if (assaytype == 'fusion'){
	gene.list$sense = tolower(gsub(' ', '', gene.list$sense))
} else {
	## for compatability with mutation assay, assign sense 'both'
	gene.list$sense = 'both'
}

if (assaytype %in% c('fusion', 'mutation')){
	panel.NM = sort(unique(gene.list$NM))
	writeLines(panel.NM, 'panel.NM')
	system(paste("join panel.NM ", depdir, "/hg19.RefGene.NM > target.refseq.0", sep=''))
	system(paste("sort -k1,1 -u target.refseq.0 > target.refseq", sep=''))
} else {
	write.table(gene.list, '_input.bed', sep='\t', quote=F, row.names=F)
}

## tail
if (previousTailFile != 'NO'){
	system(paste("cp ", previousTailFile, " ./previous.panel.t12", sep=''))
}


######################################################################
# step 1. retreive sequence from genome 
cat("step 1. retreive sequence from genome...\n")
source(paste(ampdir, '/step.1.retreive.sequence.R', sep='')) 
	

######################################################################
# step 2. call primer3 - several iteration rounds
cat("step 2. call primer3 - several iteration rounds...\n")
source(paste(ampdir, '/step.2.call.primer3.R', sep='')) 


######################################################################
# setp 3. BLAT candidate primers against genome
cat("setp 3. BLAT candidate primers against genome...\n")
	t = 11
	s = 4
	n.srv = ceiling(ncpu/4)
	ClnPsrv = 3
	n.fa=n.srv*ClnPsrv
	repMatch=1024
	start.port=9000
	faDir = 'split'
	depdir
	pjdir
	db = 'hg19.2bit'
	minIden = 93
	minScore = 12
	maxIntron = 700001
	NAtype = 'dna'
source(paste(ampdir, '/step.3.blat.candidate.primer.R', sep='')) 

# organize balt results
system(paste('bash ', ampdir, '/scripts/afterBlat.sh', sep=''))


######################################################################
# step 4. pair GSP1-GPS2
cat("step 4. pair GSP1-GPS2...\n")
source(paste(ampdir, '/step.4.pairing.GSPs.R', sep='')) 

	## top candidates
	system('head -n1 ranked.all.candidate.pairs > ranked.1st.pairs')
	system("grep ' 1$' ranked.all.candidate.pairs >> ranked.1st.pairs")


######################################################################
# step 5. check uniqueness of 12 bases at 3' of all primers. 
#      - If not unique, candidate pairs from 'ranked.all.candidate.pairs'
#	    will be retreived, ranked 2/3/.., until all tail 12 bases
#           are unique.
#
#      - If all unique, proceed.
#
#      - Save final tail 12 bases of all primers for future use.
#          e.g, to check uniqueness when adding new primers 
#	   to an existing panel.
#
cat("step 5. check uniqueness of 12 bases at 3' of all primers 
    and final output...\n")
pairs1 = read.table('ranked.1st.pairs', header=T, stringsAsFactors=F)
pairs1$target = sub(":.*", "", pairs1$r1.qName)
source(paste(ampdir, '/step.5.check.uniqueness.R', sep='')) 

writeLines(all.t12, paste(panel, '.tail.12bases.'
	  , format(Sys.time(), "%Y.%b.%d"), sep=''))
	   

######################################################################
# Lastly,
#  - add GSP tags
#  - output 
#  - save primer.bed
######################################################################
pairs.keep$GSP1 = paste(GSP1tag, pairs.keep$r1.seq,sep='')
pairs.keep$GSP2 = paste(GSP2tag, pairs.keep$r2.seq,sep='')
pairs.keep$gene = sapply(strsplit(pairs.keep$target, '_'), "[[", 1)
pairs.keep$exon = sapply(strsplit(pairs.keep$target, '_'), "[[", 2)
pairs.keep$exonSize = sapply(strsplit(pairs.keep$target, '_'), "[[", 3)
if (assaytype == 'fusion'){
	pairs.keep$sense = sapply(strsplit(pairs.keep$target, '_'), "[[", 4)
} else {
	pairs.keep$sense = pairs.keep$exonSize
}


if (assaytype %in% c('fusion', 'mutation')){
	gene.NM = ref[,c('name2','name')]
	names(gene.NM) = c('gene', 'NM')
	final = merge(pairs.keep, gene.NM, by='gene')
	final$gsp1.name = paste(final$gene, '_ex', final$exon, '_', final$sense
			  , '.1 (', final$NM, ')', sep='')
	final$gsp2.name = paste(final$gene, '_ex', final$exon, '_', final$sense
			  , '.2 (', final$NM, ')', sep='')
} else {
	final = pairs.keep
	final$gsp1.name = paste(final$gene, '_', final$sense
			  , '.1', sep='')
	final$gsp2.name = paste(final$gene, '_', final$sense
			  , '.2', sep='')
}



# output all exons
	gsp1 = final[, c('gsp1.name', 'GSP1')]
	gsp2 = final[, c('gsp2.name', 'GSP2')]

	# primer bed - only GSP2 is sequenced/relevant
	primer.bed = final[, c('r1.chr', 'r2.tStart', 'r2.tEnd', 'gsp2.name')]

	write.csv(gsp1, paste(panel, '_all.gsp1.csv', sep=''), quote=F, row.names=F)
	write.csv(gsp2, paste(panel, '_all.gsp2.csv', sep=''), quote=F, row.names=F)
	write.table(primer.bed, paste(panel, '_all.gsp2.primer.bed', sep=''), sep='\t'
		    , quote=F, row.names=F, col.names=F)

# selected exon/sense
if (assaytype != 'bed'){

	if (assaytype == 'fusion'){
		final$nm.sense = paste(final$NM, tolower(final$sense), sep='_')
		gene.list$nm.sense = paste(gene.list$NM, gene.list$sense, sep='_')
		nm.sense = unique(gene.list$nm.sense)
		final$sense.select = 0
		final$sense.select[final$nm.sense %in% nm.sense] = 1
	}

	final$nm.exon = paste(final$NM, as.numeric(final$exon), sep='_')
	gene.list$nm.exon = paste(gene.list$NM, gene.list$exon, sep='_')
	nm.exon = unique(gene.list$nm.exon)
	final$exon.select = 0
	final$exon.select[final$nm.exon %in% nm.exon] = 1
	nm.all.exons = gene.list$NM[gene.list$exon=='all']
	final$exon.select[final$NM %in% nm.all.exons] = 1

	final$select = 0
	if (assaytype == 'fusion'){
		final$select[final$exon.select ==1 & final$sense.select ==1] = 1
	} else {
		final$select[final$exon.select ==1] = 1
	}

	final.s = subset(final, select==1)
} else {
	final.s = final
}
	## backup final data
	 write.table(final, paste(panel, '_intermediate.data.txt', sep=''), sep='\t'
		     , quote=F, row.names=F)
	

	##############################
	## final output
	##############################
	GSP1 = final.s[, c('gsp1.name', 'GSP1')]
	GSP2 = final.s[, c('gsp2.name', 'GSP2')]

	# primer bed - only GSP2 is sequenced/relevant
	primer.bed.s = final.s[, c('r2.chr', 'r2.tStart', 'r2.tEnd', 'gsp2.name')]

	write.csv(GSP1, paste(panel, '_GSP1.csv', sep=''), quote=F, row.names=F)
	write.csv(GSP2, paste(panel, '_GSP2.csv', sep=''), quote=F, row.names=F)
	write.table(primer.bed.s, paste(panel, '_GSP2.primer.bed', sep=''), sep='\t'
		    , quote=F, row.names=F, col.names=F)


## cleanup
	system('mv unpaired.targets  targets.failed.pairing.GSPs.txt')
	if (file.exists('missed.seq_seq_primer3.exome.setting.3')){
		system('mv missed.seq_seq_primer3.exome.setting.3  targets.failed.primer3.design.txt')
	}

## comment out below to keep intermediate results for troubleshooting.
	system('rm -rf split seq seq.noMask out all.psl.matt.sorted primer.pos.tm.sorted primer.psl missed.seq_* splitFa.Rout panel.NM primer.candidates.0 primer3.exome.setting*')

cat("==========\nAMP-primer-design completed successfully.\n==========\n")
## END

