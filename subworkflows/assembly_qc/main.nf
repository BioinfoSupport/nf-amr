
include { MINIMAP2_ALIGN_ONT } from '../../modules/minimap2/align_ont'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_LONG  } from './modules/samtools/stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_SHORT } from './modules/samtools/stats'
include { BWA_MEM            } from './modules/bwa/mem'
include { BWA_INDEX          } from './modules/bwa/index'
include { ORGANIZE_FILES     } from '../../modules/organize_files'
include { RMD_RENDER         } from '../../modules/rmd/render'

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
			// Short reads alignment and statistics
			BWA_MEM(BWA_INDEX(fa_ch).join(fqs_ch))
			SAMTOOLS_STATS_SHORT(BWA_MEM.out.bam)
			
			// Long reads alignment and statistics
			MINIMAP2_ALIGN_ONT(fa_ch.join(fql_ch))
			SAMTOOLS_STATS_LONG(MINIMAP2_ALIGN_ONT.out.bam)
			
			//RUN VCF_LONG
			//RUN VCF_SHORT
			//HTML_AND_JSON_QC_REPORT()
			
			fa_ch
				.join(SAMTOOLS_STATS_LONG.out,remainder:true)
				.join(SAMTOOLS_STATS_SHORT.out,remainder:true)
				.map({meta,x1,x2,x3 -> [meta,[[x1,"assembly.fasta"],[x2,"long_reads.bam.stats"],[x3,"short_reads.bam.stats"]]]})
				| ORGANIZE_FILES
			RMD_RENDER(
				ORGANIZE_FILES.out.map({m,x -> [m,x,"isolate_dir='${x}'"]}),
				file("${moduleDir}/assets/isolate_assembly_qc.Rmd"),
				file("${moduleDir}/assets/isolate_assembly_qc.Rmd")
			)

			//CHARACTERIZE_UNMAPPED_READS
			//QC_AGGREGATOR
	emit:
		long_bam        = MINIMAP2_ALIGN_ONT.out.bam
		long_bai        = MINIMAP2_ALIGN_ONT.out.bai
		long_bam_stats  = SAMTOOLS_STATS_LONG.out
		long_vcf        = Channel.empty()
		
		short_bam       = BWA_MEM.out.bam
		short_bai       = BWA_MEM.out.bai
		short_bam_stats = SAMTOOLS_STATS_SHORT.out
		short_vcf       = Channel.empty()
		
		html            = RMD_RENDER.out.html
}



