
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
Long Reads
  - % mapped
  - Avg coverage

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
			//HTML_AND_JSON_QC_REPORT()
			//CHARACTERIZE_UNMAPPED_READS
			//QC_AGGREGATOR
	emit:
		cram_long = MINIMAP2_ALIGN_ONT.out.cram
		crai_long = MINIMAP2_ALIGN_ONT.out.crai
		cram_stats_long = SAMTOOLS_STATS_LONG.out
		
		cram_short = BWA_MEM.out.cram
		crai_short = BWA_MEM.out.crai
		cram_stats_short = SAMTOOLS_STATS_SHORT.out
}



