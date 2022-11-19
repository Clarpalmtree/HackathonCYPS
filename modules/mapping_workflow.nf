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
        
        STAR --outSAMstrandField intronMotif \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --genomeDir ${ref} \
        --readFilesIn ${r1} ${r2}
        --runThreadN ${task.cpus} \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 50000000000 \
        > ${id}.bam
        """
}

workflow MAPPING {
    take:
    fastq_files
    index
    main:
    mapping(fastq_files, index)
}
