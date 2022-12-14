nextflow.enable.dsl=2

// setting params
params.project = "SRA062359" // sra project number
params.resultdir = 'sra_info' // results output directory
params.resultdir2 = 'fastq' // results output directory
params.resultdir3 = 'genome' // results output directory
params.resultdir4 = 'annotation' // results output directory
params.resultdir5 = 'index' // results output directory

process getSRAIDs {
    // Getting SRA IDs of data of interest in a txt file

    publishDir params.resultdir, mode: 'copy' 

    input:
    val projectid

    output:
    file 'sra.txt'

    script:
    """
    esearch -db sra -query $projectid  | efetch --format runinfo | grep TRANSCRIPTOMIC | cut -d ',' -f 1 > sra.txt
    """
}

process getSRA {
    // Downloadind .sra files
    publishDir params.resultdir, mode: "copy"

    input:
    val id

    output:
    file '*.sra'

    script:
    """
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/${id}/${id} -O ${id}.sra
    """
}

process fastqDump {
    // Downloading fastq files

    publishDir params.resultdir2, mode: 'copy'

    input:
    val id
    file sra_file

    output:
    tuple val(id), path("*_1.fastq"), path("*_2.fastq")

    script:
    """
    fastq-dump --gzip --split-files ${sra_file}
    gunzip *.fastq.gz
    """
}

process getGenome {
    // Downloading each chromosome genome file
    publishDir params.resultdir3, mode: 'copy'
    output:
    file 'ref.fa' // unique file with all chromosomes
    script: 
    """
    wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz
    gunzip -c *.fa.gz > ref.fa
    rm *.fa.gz
    """
}

process getAnnot {
    // Getting annotation file and unzipping it for the index process

    publishDir params.resultdir4, mode: 'copy'

    output:
    file 'Homo_sapiens.GRCh38.101.chr.gtf'

    script:
    """
    #Getting genome annotations
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip -c Homo_sapiens.GRCh38.101.chr.gtf.gz > Homo_sapiens.GRCh38.101.chr.gtf
    rm Homo_sapiens.GRCh38.101.chr.gtf.gz
    """
}

process index{
    // indexing with STAR

	publishDir params.resultdir5, mode: 'copy'

	input:
	file gen
	file annot

	output:
	path 'index_files/' // indexing files are stored in ref/ directory

	script: 
	"""
	mkdir index_files
	STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir index_files --genomeFastaFiles ${gen} --sjdbGTFfile ${annot}
	"""
}


workflow COLLECT {
    
    main:
    // Getting sra ids
    getSRAIDs(params.project)
    sraID = getSRAIDs.out.splitText().map { it -> it.trim() }
    sraID.view()
    // Testing with only 2 samples
    // sraID=Channel.from('SRR628582','SRR628583')

    // Getting sra files
    getSRA(sraID)

    // Getting fastq files
    fastq=fastqDump(sraID,getSRA.out)

    // Getting chromosome files and reference genome
    getGenome()

    // Annotation file
    annot=getAnnot()

    // Indexation
    ind=index(getGenome.out, getAnnot.out)

    emit:
    fastq
    ind
    annot

}
