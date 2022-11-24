nextflow.enable.dsl=2

// setting params
params.resultdir = 'stat_analysis' // results output directory

process stat_analysis {
    publishDir params.resultdir, mode: 'copy'

    input:
    file 'output.counts'

    output:
    tuple file('PCA_GraphOfIndividuals.pdf'), file('DESeq_results.txt'), file('plot_counts.pdf'), file ('heatmap_MostVariableGenes.pdf'), file ('MostVariableGenes.txt'),
          file ('Significative_DEgenes.txt'), file('Significative_DEgenes_Summary.txt') into ana_stat //Récupère les résultats de l'analyse statistique

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
