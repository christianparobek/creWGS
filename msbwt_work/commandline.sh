python ~/lib/python/count_reference_kmers.py --regions pfmdr_region.bed --genome /proj/julianog/refs/Pf3D7_v9.3/PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta --kmer 30 --step 5 --msbwt Pf_Hi_Pool_msbwt/ --verbose --revcomp  > Pf_HI_pfmdr.out.gene

msbwt cffq CRE007 data/CRE007_TCCGCGA-TATAGCC_L001_R1_001.fastq.gz data/CRE007_TCCGCGA-TATAGCC_L001_R2_001.fastq.gz -p 2 -u
msbwt cffq CRE011 data/CRE011_TCCGCGA-AGGCGAA_L001_R1_001.fastq.gz data/CRE011_TCCGCGA-AGGCGAA_L001_R2_001.fastq.gz -p 2 -u

python count_reference_kmers.py --region region.bed --genome ref/CRE07_fixed.fasta --kmer 30 --step 5 --msbwt CRE007/ > CRE007.results
python count_reference_kmers.py --region region.bed --genome ref/CRE07_fixed.fasta --kmer 30 --step 5 --msbwt CRE011/ > CRE011.results
	## CRE007 and CRE011 are pretty much identical so it should be ok to use CRE07 fasta for both
## just need to analyze the right region, potentially normalize, etc!


