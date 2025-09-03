
process SPECIATOR {
    container 'speciator-v4.0.0.sif'
    memory '10 GB'
    cpus 1
    time '15 min'
    input:
        tuple val(meta), path('assembly.fna')
    output:
    		tuple val(meta), path("speciator",type: 'dir'), emit: orgfinder
    script:
		    """
		    mkdir ./speciator
		    python3 /speciator.py assembly.fna /libraries /bactinspector/data/taxon_info.pqt > speciator/results.json
		    """
}

