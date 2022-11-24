nextflow.enable.dsl=2

include {COLLECT} from './modules/collect_workflow.nf'
include {MAPPING} from './modules/mapping_workflow.nf'
include {STAT} from './modules/stat_analysis_workflow.nf'

workflow{
    
    COLLECT()
    fastq_files=COLLECT.out[0]
    index=COLLECT.out[1]
    annot=COLLECT.out[2]
    MAPPING(fastq_files, index, annot)
    matrix=MAPPING.out[1]
    STAT(matrix)

}
