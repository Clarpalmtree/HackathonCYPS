nextflow.enable.dsl=2

// setting params
params.project = "SRA062359" // sra project number
params.resultdir = 'mapping' // results output directory

process mapping {
        publishDir params.resultdir, mode: 'copy'

        input:
        tuple val (id), path (r1), path (r2)
        file ref

        output:
        path '*.bam'

        script :
        """
        STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 --genomeDir ${ref} --readFilesIn ${r1} ${r2} --runThreadN ${task.cpus} --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --genomeLoad NoSharedMemory > ${id}.bam
        """
}

process mapping_bai {

	publishDir params.resultdir, mode: 'copy'
	
	input:
	path bam // bam files

	output:
	file '*.bai'

	script:
	"""
	samtools index ${bam} nthreads=${task.cpus}
	"""
}

process featureCounts {

    publishDir params.resultdir, mode: 'copy'

    input:
    path bam // bam files
    file gen // genome

    output:
    tuple file ('output.counts'), file ('output.counts.summary') //recupere la matrice de comptage et un résumé de l’attribution des reads

    script:
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a ${gen} -o output.counts ${bam}
    """
}

workflow MAPPING {
    take:
    fastq_files
    index
    genome
    
    main:
    map=mapping(fastq_files, index)
    mapping_bai(map)
    matrix=featureCounts(map,genome)
}
