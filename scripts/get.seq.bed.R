
## This script is to be called by 'step.1.retreive.sequence.R' to:
##   - retreive sequences from genome (default hg19/GRCh37) based on user
##       requested chromosomal coordinates (bed file).


# === Requires:
#  input.bed - the bed file
#  target - 4th column of bed file
#  tempsize - size of template sequences to retreive
#  depdir - dependency data path
#  subExsize - for exons larger than which, design titling primers
#  leadsize - the distance between retreive sequence template 
#               (normally intronic) and target (normally exons)
#               

# === Should have existed:
#  seq/
#  seq.noMask/
#  target.bed

## usage e.g.:
# Rscript --vanilla get.seq.bed.R  target tempsize depdir subExonSize leadsize

library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP.20111119)
bed = read.table('_input.bed', header=T, stringsAsFactors=F)

target = commandArgs(TRUE)[1]
tempsize = as.numeric(commandArgs(TRUE)[2])
depdir = commandArgs(TRUE)[3]
subExonSize = as.numeric(commandArgs(TRUE)[4]) 
leadsize = as.numeric(commandArgs(TRUE)[5])

system(paste("touch _running_", target, sep=''))

   print(target)
   # make label strings 001, 002...
   iii = substr(1001:1999, 2, 5)

   targetTab = bed[bed$target==target,]
   chrom = targetTab$chrom
   start = targetTab$start -1
   end = targetTab$end +1
   size = end - start
   n.exon = 1

   ## read dbsnp
   chrn = toupper(sub('chr', '', chrom))
   snp.chr = read.table(paste(depdir, '/dbsnp/snp.', chrn, sep=''), sep='\t', header=F, stringsAsFactors=F)
   

   ##################
           ex.sub.n = ceiling(size / subExonSize)
           ex.sub.size = ifelse(ex.sub.n == 1, 0, floor(size/ex.sub.n) - 1)

       for (ex.sub in ex.sub.n:1){
	   ########################################
	   ### 3' of exon on plus (fwd) strand
	   ########################################
	   end.ex.sub = end - ex.sub.size * (ex.sub - 1) 
			+ leadsize # at least 'leadsize' bp in intron
	   end.ex.sub.tempsize = end.ex.sub + tempsize
	   end.ex.sub2 = end.ex.sub.tempsize
	   end.ex.sub.TplSize = end.ex.sub2 - end.ex.sub


    	if (end.ex.sub.TplSize > 30){
	   ex.3f = data.frame('getSeq'=getSeq(Hsapiens, chrom
			, end.ex.sub, end.ex.sub2, as.character=T
			, strand='+')
			, 'size' = size
			, 'sub' = letters[ex.sub]
			, 'strand'='+')

	   ## before mask dbsnp
	   ex.3f$seq = paste('SEQUENCE_ID=', target, '_000_'
				, ex.3f$size, ex.3f$sub
				, ex.3f$strand
			, '\nSEQUENCE_TEMPLATE=', ex.3f$getSeq
			, '\n=', sep='')
	   write.table(ex.3f$seq, file=paste('seq.noMask/', target
			, '_000_', ex.3f$size, ex.3f$sub
			, ex.3f$strand, sep='')
			, row.names=F, quote=F, col.names=F)

	   ## mask dbsnp
	   snp = subset(snp.chr, V2 >= end.ex.sub 
				& V2 <= end.ex.sub + tempsize)
	   snpn = nrow(snp)
	   seq = ex.3f$getSeq
	   if (snpn>0){
		   for (i in 1:snpn){
		   pos = snp$V2[i] - end.ex.sub +1
		   Ns = nchar(snp$V3[i])
		   l.seq = substr(seq, 1, pos -1)
		   r.seq = substr(seq, pos + Ns, end.ex.sub +
				 tempsize + 1)
		   seq = paste(l.seq, paste(rep('N', Ns), collapse='')
				, r.seq, sep='')

			}
	   ex.3f$getSeq = seq
	   }

	   ex.3f$seq = paste('SEQUENCE_ID=', target, '_', '000'
			, '_', ex.3f$size, ex.3f$sub, ex.3f$strand
			, '\nSEQUENCE_TEMPLATE=', ex.3f$getSeq
			, '\n=', sep='')
	   write.table(ex.3f$seq, file = paste('seq/', target, '_'
			, '000', '_', ex.3f$size, ex.3f$sub
			, ex.3f$strand, sep='')
			, row.names=F, quote=F, col.names=F)

    	} else {
		cat(c(target, '000', strand, ex, ex.sub, end.ex.sub.TplSize
			, '\n'), file = 'smallTplSize', append=T)
               }

		
	   ########################################
      	   # 5' of exon on the minus (revers) strand
	   ########################################
	   start.ex.sub = start + ex.sub.size * (ex.sub - 1) 
				- leadsize 
	   start.ex.sub.tempsize = start.ex.sub - tempsize
	   start.ex.sub2 = start.ex.sub.tempsize
	   start.ex.sub.TplSize = start.ex.sub - start.ex.sub2

	   if (start.ex.sub.TplSize > 30){

	   ex.5r = data.frame('getSeq' = getSeq(Hsapiens, chrom
			, start.ex.sub2, start.ex.sub, as.character=T
			, strand = '-')
			, 'size' = size
			, 'sub' = letters[ex.sub]
			, 'strand' = '-')

	   ## before mask dbsnp
	   ex.5r$seq = paste('SEQUENCE_ID=', target, '_', '000'
			, '_', ex.5r$size, ex.5r$sub, ex.5r$strand
			, '\nSEQUENCE_TEMPLATE=', ex.5r$getSeq
			, '\n=', sep='')
	   write.table(ex.5r$seq, file = paste('seq.noMask/', target
			, '_', '000', '_', ex.5r$size, ex.5r$sub
			, ex.5r$strand, sep='')
			, row.names = F, quote = F, col.names=F)

	   ## mask dbsnp
	   snp5 = subset(snp.chr, V2 >= (start.ex.sub2) 
				& V2 <= start.ex.sub)
	   snp5n = nrow(snp5)
	   seq5 = ex.5r$getSeq
	   if (snp5n>0){
		   for (i in 1:snp5n){
		   pos = start.ex.sub - snp5$V2[i] +1
		   Ns = nchar(snp5$V3[i])
		   l.seq = substr(seq5, 1, pos - Ns)
		   r.seq = substr(seq5, pos + 1, start.ex.sub.TplSize)
		   seq5 = paste(l.seq, paste(rep('N', Ns), collapse='')
				, r.seq, sep='')
		   }
	   ex.5r$getSeq = seq5
	   }
	   ex.5r$seq = paste('SEQUENCE_ID=', target, '_', '000', '_'
			, ex.5r$size, ex.5r$sub, ex.5r$strand
			, '\nSEQUENCE_TEMPLATE=', ex.5r$getSeq
			, '\n=',sep='')
	   write.table(ex.5r$seq, file = paste('seq/', target, '_', '000'
			, '_', ex.5r$size, ex.5r$sub, ex.5r$strand, sep = '')
			, row.names = F, quote = F, col.names = F)

	   }else{
	   	cat(c(target, '000', strand, ex, ex.sub, end.ex.sub.TplSize
			, '\n'), file = 'smallTplSize', append=T)
	   }
}

######################################
### write exon bed + (leadsize - 2)
######################################

        bed0 = paste(chrom, start - leadsize + 1, end + leadsize - 1
			, target, size, ex.sub.n, sep = '\t')
        write.table(bed0, 'target.bed0', append = T, quote = F, col.names = F
			, row.names = F)


system(paste("rm _running_", target, sep=''))
