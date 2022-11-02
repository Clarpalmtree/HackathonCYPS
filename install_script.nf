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

process fastqDump {
	
	publishDir params.resultdir, mode: 'copy'

	input:
	val id //pour chaque numero SRR

	output:
	tuple val (id),file('*1.fastq.gz'), emit: reads_1  // tuple (sraID, read1)
    tuple val (id),file('*2.fastq.gz'), emit: reads_2  // tuple (sraID, read2)

	script:
	"""
    fasterq-dump ${id} --threads 4 --split-files
	"""	
}

workflow {
    projectID=params.project
    getSRAIDs(projectID)
    sraID = getSRAIDs.out.splitText().map { it -> it.trim() }
    sraID.view()
    fastqDump(sraID)
    readsTable = fastqDump.out.reads_2.join(fastqDump.out.reads_1) //ajoute les SRR aux tuples des reads pour former des tuples : (file reads1,  file reads2, val SRR)
}