nextflow.enable.dsl=2

// setting params
params.project = "SRA062359" // sra project number
params.resultdir = 'mapping_bam' // bam directory
params.resultdir2 = 'mapping_bai' // bam.bai directory
params.resultdir3 = 'count' // couting matrix directory

process mapping {
    // Mapping reads with STAR
    publishDir params.resultdir, mode: 'copy'

    input:
    tuple val (id), path (r1), path (r2)
    file index_files

    output:
    path '*.bam'

    script :
    """
    STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 --genomeDir ${index_files} --readFilesIn ${r1} ${r2} --runThreadN ${task.cpus} --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --genomeLoad NoSharedMemory --limitBAMsortRAM 50000000000 > ${id}.bam
    """
}

process mapping_bai {
    // Indexing bam fils with samtools
	publishDir params.resultdir2, mode: 'copy'
	
	input:
	path bam // bam files

	output:
	file '*.bai'

	script:
	"""
	samtools index ${bam}
	"""
}

process featureCounts {
    // Getting counting matrix with feartureCounts
    publishDir params.resultdir3, mode: 'copy'

    input:
    path bam // bam files
    file annot // annotation
    path bai // getting indexed bam files

    output:
    tuple file ('output.counts'), file ('output.counts.summary') // counting matrix and read attribution summary

    script:
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a ${annot} -o output.counts ${bam}
    """
}

workflow MAPPING {

    take:
    fastq_files
    index_files
    annot
    
    main:
    mapping(fastq_files, index_files)
    mapping_bai(mapping.out)
    matrix=featureCounts(mapping.out.collect(),annot,mapping_bai.out.collect())
    
    emit:
    matrix[0]
}
