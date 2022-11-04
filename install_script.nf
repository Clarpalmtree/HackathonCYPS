// COLLECT WORKFLOW

nextflow.enable.dsl=2

// setting params
params.project = "SRA062359" // sra project number
params.resultdir = 'results' // results output directory

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

        publishDir params.resultdir, mode: 'copy'

        input:
        val id
        file sra_file

        output:
    tuple val(id), path("*_1.fastq.gz"), path("*_2.fastq.gz")

        script:
        """
    fastq-dump --gzip --split-files ${sra_file}
        """
}

process chromosome {
    // Downloading each chromosome genome file
    publishDir params.resultdir, mode: 'copy'

    input:
    val chr

    output:
    file 'Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz'

    script: 
    """
    wget -o ${chr}.fa.gz "ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz"
    """
}


process mergechr {
    // Merging all chromosome files into a single 'ref.fa' file (=the reference genome)
    publishDir params.resultdir, mode: 'copy'

    input:
    file allchr

    output:
    file 'ref.fa' // unique file with all chromosomes

    script:
    """
    gunzip -c ${allchr} > ref.fa
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
	path 'ref/' //renvoie un unique repertoire contenant tous les fichiers de l'index de reference

	script: 
	"""
	mkdir ref
	chmod +x ref
	STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ${gen} --sjdbGTFfile ${annot}
	"""
}

process mapping {
	publishDir params.resultdir, mode: 'copy'

	input:
	tuple val (id), file (r2), file (r1)  
	path ref 

	output:
	file '*.bam', emit: bamind		//recupere les fichiers bam pour l'indexation samtools
	file '*.bam', emit: bamcount	//recupere les fichiers bam pour le comptage

	script :
	"""
	STAR --outSAMstrandField intronMotif \
	--outFilterMismatchNmax 4 \
	--outFilterMultimapNmax 10 \
	--genomeDir ${ref}\
	--readFilesIn <(gunzip -c ${r1}) <(gunzip -c ${r2}) \
	--runThreadN 6 \
	--outSAMunmapped None \
	--outSAMtype BAM SortedByCoordinate \
	--outStd BAM_SortedByCoordinate \
	--genomeLoad NoSharedMemory \
	--limitBAMsortRAM 50000000000 \
	> ${id}.bam
	"""
}

process mappingbai{

	publishDir params.resultdir, mode: 'copy'
	
	input:
	file bamind

	output:
	file '*.bai', emit: map

	script:
	"""
	samtools index ${bamind}
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
    // get sra files
    getSRA(sraID)
    // get fastq files
    fastqDump(sraID,getSRA.out)
    //chr
    chromosome(list)
    mergechr(chromosome.out)
    //annot
    getAnnot()
    // Indexation
    index(mergechr.out, getAnnot.out)
    // mapping(sraID, reads_1, reads_2)
    // mappingbai(bamind)

}
