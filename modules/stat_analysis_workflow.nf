nextflow.enable.dsl=2

// setting params
params.resultdir = 'stat_analysis' // results output directory


process mapping {
        publishDir params.resultdir, mode: 'copy'

        input:
        tuple val (id), path (r1), path (r2)
        file index_files

        output:
        path '*.bam'

        script :
        """
        STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 --genomeDir ${index_files} --readFilesIn ${r1} ${r2} --runThreadN ${task.cpus} --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --genomeLoad NoSharedMemory > ${id}.bam
        """
}
