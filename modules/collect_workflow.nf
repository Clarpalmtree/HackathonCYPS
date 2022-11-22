nextflow.enable.dsl=2

// setting params
params.project = "SRA062359" // sra project number
params.resultdir = 'data' // results output directory

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

        publishDir params.resultdir, mode: 'copy'

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
    publishDir params.resultdir, mode: 'copy'
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
    // indexing with STAR

	publishDir params.resultdir, mode: 'copy'

	input:
	file gen
	file annot

	output:
	path 'ref/' // indexing files are stored in ref/ directory

	script: 
	"""
	mkdir ref
	STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ${gen} --sjdbGTFfile ${annot}
	"""
}


workflow COLLECT {
    

    main:
    	chr_list = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'MT')
	
        // Getting sra ids
        getSRAIDs(params.project)
        sraID = getSRAIDs.out.splitText().map { it -> it.trim() }
        sraID.view()

        // Testing with only 2 samples
        // sraID=Channel.from('SRR628582','SRR628583')

        // get sra files
        getSRA(sraID)

        // get fastq files
        fastq=fastqDump(sraID,getSRA.out)

        // get chromosome files and reference genome
        genome=getGenome()

        //annotation file
        getAnnot()

        // Indexation
        ind=index(getGenome.out, getAnnot.out)

        emit:
        fastq
        ind
	genome

}
