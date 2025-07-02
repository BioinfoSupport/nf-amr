process KRAKEN2_CLASSIFY {
    container 'docker.io/staphb/kraken2:2.1.5'
    memory '5 GB'
    cpus 4
    time '30 min'
    input:
    		each path("db")
        tuple val(meta), path(seq)
    output:
    		tuple val(meta), path('kraken2',type:'dir'), emit: report
    script:
		    """
		    mkdir -p kraken2
		    k2 classify --db db ${task.ext.args?:''} --report kraken2/report.txt --report-zero-counts --threads ${task.cpus} ${seq.extension=='gz'?'--gzip-compressed':''} ${seq} > kraken2/output.tsv
		    """
}

// The output can be processed with this script also:
//cat output.tsv | cut -f5 | tr " :" "\n\t" | awk '{n[$1]+=$2;N+=$2}END{for(v in n){print v,n[v],n[v]/N}}' | sort -k2,2n | tail -n30
