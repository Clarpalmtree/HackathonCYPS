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

process chromosome {

    publishDir params.resultdir, mode: 'copy'

    input:
    val chr

    output:
    path 'Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz', emit: chrfasta //place tous les chromosomes telecharges dans 1 channel

    script: 
    """
    wget -o ${chr}.fa.gz "ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz"
    """
}


process mergechr {
	
    publishDir params.resultdir, mode: 'copy'

    input:
    file allchr from chrfasta.collect() //attend que tous les chromosomes soient telecharges

    output:
    file 'ref.fa' into fasta //unique fichier contenant tous les chromosomes

    script:
    """
    gunzip -c ${allchr}> ref.fa
    """
}

process gtf {
    
    output:
    file 'annot.gtf' into human_genome

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
    gunzip -c Homo_sapiens.GRCh38.104.chr.gtf.gz > annot.gtf
    """
}

process index{
	publishDir params.resultdir, mode: 'copy'

	input:
	file c from fasta
	file annot from human_genome

	output:
	file 'ref/' into index //renvoie un unique repertoire contenant tous les fichiers de l'index de reference

	script: 
	"""
	mkdir ref
	chmod +x ref
	STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ${c} --sjdbGTFfile ${annot}
	"""
}



workflow {
    projectID=params.project
    list = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'MT')
    getSRAIDs(projectID)
    sraID = getSRAIDs.out.splitText().map { it -> it.trim() }
    sraID.view()
    fastqDump(sraID)
    readsTable = fastqDump.out.reads_2.join(fastqDump.out.reads_1) //ajoute les SRR aux tuples des reads pour former des tuples : (file reads1,  file reads2, val SRR)
    chromosome(list)

}