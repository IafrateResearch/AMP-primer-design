
## To generate dbSNP database.

## Note: reference version, GRCh37 or GRCh38, should be consistent between 
## 	dependency data, input BED file (if used) and pipeline (get.seq.###.R).
## 	The default version of this pipeline is GRCh37 (hg19).

# run this under dependency-AMP folder (e.g. ~/dependency-AMP).
cd ~/dependency-AMP

# download
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/clinical_vcf_set/clinvar_20150305.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/common_all_20150217.vcf.gz

# unzip
gunzip clinvar_20150305.vcf.gz
gunzip common_all_20150217.vcf.gz

# get nessesary info
cut -f1,2,4,5 clinvar_20150305.vcf | grep -v '^#' > _TMP.f4
cut -f1,2,4,5 common_all_20150217.vcf | grep -v '^#' >> _TMP.f4

sort -k1,1n -k2,2n -k3,4 -u _TMP.f4 > clinvar_20150305_plus_common_20150317.txt

rm -rf dbsnp ## in case exists
mkdir dbsnp
awk '{print $0 >> "dbsnp/snp."$1}' clinvar_20150305_plus_common_20150317.txt

# cleanup
rm _TMP*
