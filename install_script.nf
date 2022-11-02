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
    file 'Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz' //place tous les chromosomes telecharges dans 1 channel

    script: 
    """
    wget -o ${chr}.fa.gz "ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz"
    """
}


process mergechr {
	
    publishDir params.resultdir, mode: 'copy'

    input:
    file allchr

    output:
    file 'ref.fa'//unique fichier contenant tous les chromosomes

    script:
    """
    gunzip -c ${allchr} > ref.fa
    """
}

process getAnnot {
    // Getting annotation

    publishDir params.resultdir, mode: 'copy'

    output:
    file 'Homo_sapiens.GRCh38.101.chr.gtf'

    script:
    """
    #Getting genome annotations
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gzip -d Homo_sapiens.GRCh38.101.chr.gtf.gz
    """
}

process index{
	publishDir params.resultdir, mode: 'copy'

	input:
	file c 
	file annot

	output:
	path 'ref/' //renvoie un unique repertoire contenant tous les fichiers de l'index de reference

	script: 
	"""
	mkdir ref
	chmod +x ref
	STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ${c} --sjdbGTFfile ${annot}
	"""
}



workflow {
    // Channels
    projectID=params.project
    list = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'MT')
    // sraid
    getSRAIDs(projectID)
    sraID = getSRAIDs.out.splitText().map { it -> it.trim() }
    sraID.view()
    //chr
    chromosome(list)
    mergechr(chromosome.out)
    //annot
    getAnnot()
    // Indexation
    index(mergechr.out, getAnnot.out)

}