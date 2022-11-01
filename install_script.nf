nextflow.enable.dsl=2
// setting params
params.project = "SRA062359" // sra project number
params.resultdir = 'results' // results output directory

process getSRAIDs {

        publishDir params.resultdir, mode: 'copy' //Les résultats sont copié dans le dossier 'params.resultdir'

        input:
        val projectid

        output:
        file 'sra.txt'//Récupération des numéros SRR dans un fichier txt

        script:
        """
        esearch -db sra -query $projectid  | efetch --format runinfo | grep TRANSCRIPTOMIC | cut -d ',' -f 1 > sra.txt
        """
}

workflow {
    getSRAIDs(params.project)
}