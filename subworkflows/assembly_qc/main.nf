
include { MINIMAP2_ALIGN_ONT } from '../../modules/minimap2/align_ont'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_LONG  } from '../../modules/samtools/stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_SHORT } from '../../modules/samtools/stats'
include { BWA_MEM            } from '../../modules/bwa/mem'
include { BWA_INDEX          } from '../../modules/bwa/index'


/* Assembly stats
Number of contigs
Total assembly length
N50
GC content
Longest contig size
Shortest contig size
Short Reads
  - % mapped
  - % properly paired
  - Avg coverage
  - Number chimeric pair
  - Number of identified mutation in the VCF
Long Reads
  - % mapped
  - Avg coverage
  - Number of identified mutation in the VCF
*/

/* Contigs stats
Contig length
GC content
Short Reads
  - % mapped
  - % properly paired
  - Avg coverage
Long Reads
  - % mapped
  - Avg coverage
*/

/* Report
Short Reads
   - inter contigs links

*/


/*
process IGV_SCRIPT {
	input:
	  val(meta)
  output:
    tuple val(meta), path("load_in_igv.sh")
	script:
	"""
	file("${moduleDir}/assets/load_in_igv.sh")
	"""
}
*/

workflow ASSEMBLY_QC {
	take:
		fa_ch
		fql_ch
		fqs_ch
	main:
			BWA_MEM(BWA_INDEX(fa_ch).join(fqs_ch))
			MINIMAP2_ALIGN_ONT(fa_ch.join(fql_ch))
			SAMTOOLS_STATS_LONG(MINIMAP2_ALIGN_ONT.out.cram)
			SAMTOOLS_STATS_SHORT(BWA_MEM.out.cram)
			
			//RUN VCF_LONG
			//RUN VCF_SHORT
			//HTML_AND_JSON_QC_REPORT()
			//CHARACTERIZE_UNMAPPED_READS
			//QC_AGGREGATOR
	emit:
		long_cram        = MINIMAP2_ALIGN_ONT.out.cram
		long_crai        = MINIMAP2_ALIGN_ONT.out.crai
		long_cram_stats  = SAMTOOLS_STATS_LONG.out
		long_vcf         = Channel.empty()
		
		short_cram       = BWA_MEM.out.cram
		short_crai       = BWA_MEM.out.crai
		short_cram_stats = SAMTOOLS_STATS_SHORT.out
		short_vcf        = Channel.empty()
}



