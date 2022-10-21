#Work in tmp
cd /tmp/
#Downloading FASTQ files
SRAID=SRR....
wget -O ${SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${SRAID}/${SRAID}.1
fastq-dump --gzip --split-files ./${SRAID}.sra
#Downloading chromosome sequences
wget -o <chromosome>.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.!{chr}.fa.gz
gunzip -c *.fa.gz > ref.fa
#Getting genome annotations
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
