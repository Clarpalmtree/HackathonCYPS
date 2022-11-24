
<div align="center"><h1>Projet Hackathon</h1></div>
<br>
<div align="justify">
  <p>
    Dans le cadre de notre formation en M2 AMI2B à l'Université Paris-Saclay, nous avons été amenés à réaliser un workflow d'analyses RNA-Seq pour l'UE Hackathon     Reproductible. L'objectif de ce projet consiste à reproduire les résultats des analyses décrites dans ces deux articles :
    
   * [Recurrent mutations at codon 625 of the splicing factor SF3B1 in uveal melanoma](https://pubmed.ncbi.nlm.nih.gov/23313955), Harbour et al. (2013)
   * [SF3B1 mutations are associated with alternative splicing in uveal melanoma](https://pubmed.ncbi.nlm.nih.gov/23861464), Marais et al. (2013)
  </p>
</div>

<div align="left"><h2>Utilisation du workflow</h2></div>

<div align="justify">
  <p>
    Le workflow s'exécute dans Nextflow et fait appel à Docker pour les conteneurs. Nextflow est lancé depuis Conda (Bioconda). Pour exécuter le workflow, il faut     donc au préalable avoir installé Conda, Nextflow et Docker sur sa machine. La configuration du workflow proposée nécessite d'avoir au minimum 16 CPUs et 50 GB     de mémoire vive.  <br>
    Procédure à suivre pour lancer le workflow :  <br>
    Commencer par récupérer les différents scripts et accorder les permissions au répertoire bin :
    
    $ git clone git@github.com:github.com/Clarpalmtree/HackathonCYPS.git
    $ sudo chmod -R a+rwx bin
  
   Puis activer Conda et lancer le workflow : 
    
    $ conda activate nextflow
    $ nextflow run main.nf -resume
  </p>
</div>

<div align="left"><h2>Résultats</h2></div>
  <p>
  Les fichiers issus de notre analyses des données sont observables dans le dossier `/analysis` : 
  </p>
</div>

<div align="left"><h2>Matériels</h2></div>
<div align="justify">
  <p> 
  Comme expliqué précédemment, le workflow que nous proposons utilise [Nextflow](https://nextflow.io/) comme Workflow Management System. Nextflow peut être lancé   depuis [Conda](https://conda.io).
    
  Les [Dockers](https://www.docker.com/en) utilisés sont : 
     
   * [clarpalmtree/samtools](https://hub.docker.com/r/clarpalmtree/samtools) (version 1.9) : créer l'index de référence et réaliser le mapping des reads
   * [clarpalmtree/starbis](https://hub.docker.com/r/clarpalmtree/starbis) (version 2.7.10a) : pour indexer le mapping
   * [clarpalmtree/sra](https://hub.docker.com/r/clarpalmtree/sra) (version current) : pour récuperer les fichiers fastq, associer les reads 1 et 2
   * [yanismadi/fastqc](https://hub.docker.com/r/yanismadi/fastqc) (version 0.11.9) : pour faire un contrôle qualité
   * [yanismadi/subread](https://hub.docker.com/r/yanismadi/subread) (version 2.0.0) : pour faire la matrice de comptage
   * [siwarhm/dseq2](https://hub.docker.com/r/siwarhm/dseq2) (version current) : pour faire l'analyse statistique
  
  Les données biologiques utilisées dans notre étude sont celles des SRR contenues dans le fichier [SraAccList_SRA062359.txt](https://github.com//Clarpalmtree/HackathonCYPS/blob/main/SraAccList_SRA062359.txt).
  </p>
</div>


<div align="left"><h2>Auteurs</h2></div>

<div align="justify">
  <p>
Clara Toussaint, Pauline Lim, Siwar Hammami, Yanis Madi
  </p>
</div>

# Projet Hackaton 
# Année 2022-2023 - Clara Toussaint, Pauline Lim, Siwar Hammami, Yanis Madi

## Description

## Visuals

## Installation

