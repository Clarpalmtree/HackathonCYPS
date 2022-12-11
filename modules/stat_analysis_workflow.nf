nextflow.enable.dsl=2

// setting params
params.resultdir = 'stat_analysis' // analysis outputs diretory

process stat_analysis {
    // Getting analysis files with bin/stat_analysis.r script
    publishDir params.resultdir, mode: 'copy'

    input:
    file 'output.counts'

    output:
    file 'PCA_DE.png'
    file 'heatmap_de.png'
    file 'DESeq_results.csv'
    file 'significative_DEgenes.csv'
    file 'top_10_de_genes.csv'
    file 'summary.csv'

    script:
    """
    stat_analysis.r ${'output.counts'} $PWD
    """
}

workflow STAT {
    
    take:
    matrix
    
    main:
    stat_analysis(matrix)
}
