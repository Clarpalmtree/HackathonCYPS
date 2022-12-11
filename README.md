
<div align="center"><h1>Projet Hackathon</h1></div>
<br>
<div align="justify">
  <p>
    Dans le cadre de notre formation en M2 AMI2B à l'Université Paris-Saclay, nous avons été amenés à réaliser un workflow d'analyses RNA-Seq pour l'UE Hackathon Reproductible. L'objectif de ce projet consiste à reproduire les résultats des analyses décrites dans ces deux articles :
    
   * [Recurrent mutations at codon 625 of the splicing factor SF3B1 in uveal melanoma](https://pubmed.ncbi.nlm.nih.gov/23313955), Harbour et al. (2013)
   * [SF3B1 mutations are associated with alternative splicing in uveal melanoma](https://pubmed.ncbi.nlm.nih.gov/23861464), Marais et al. (2013)
  </p>
</div>

<div align="left"><h2>Utilisation du workflow</h2></div>
<div align="justify">
  <p>

  **Prérequis** 
  La configuration actuelle du workflow proposée et la grande dimension du jeu de données nécessitent un noeud de calcul ou bien un machine avec au minimum 16 CPUs et 50 GB de mémoire vive.  
  
  **Installation**

  Le workflow s'exécute dans Nextflow et fait appel à Docker pour les conteneurs. Nextflow est lancé depuis Conda (Bioconda). Pour exécuter le workflow, il faut donc au préalable avoir installé Conda, Nextflow et Docker sur sa machine. <br>
  <br>
  **Procédure à suivre pour lancer le workflow** <br>

  Commencez par récupérer les différents scripts et accorder les permissions au répertoire bin :
  
    $ git clone https://github.com/Clarpalmtree/HackathonCYPS/
    $ cd HackathonCYPS/
    $ sudo chmod -R a+rwx bin

  Puis activez Conda et lancez le workflow : 
  
    $ conda activate nextflow
    $ nextflow run main.nf

  Attention, afin que l'analyse statistique se déroule correctement, veuillez vous assurer que le fichier `metaData.txt` se retrouve bien dans votre répertoire de travail.
  
  </p>
</div>

<div align="left"><h2>Résultats</h2></div>
  <p>
  
  Les fichiers issus de notre analyse des données sont retrouvables dans le dossier `/stat_analysis` : 
  * `PCA_DE.png` : graphe PCA des 8 échantillons
  * `heatmap_de.png` : heatmap des 8 échantillons
  * `DESeq_results.csv` : fichier csv avec les ID de gène, le log Fold Change et la p-value correspondante.
  * `significative_DEgenes.csv` : fichier csv avec les gènes significativement différentiellement exprimés ordonnés par ordre de p-value ajustée.
  * `top_10_de_genes.csv` : fichier csv avec les 10 gènes les plus significativement différentiellement exprimés ordonnés par ordre de p-value ajustée.
  * `summary.csv` : fichier csv résumant les gènes différentiellement exprimés chez les échantillons Wild Type ou mutants.
  </p>
</div>


<div align="left"><h2>Matériels</h2></div>
<div align="justify">
  <p> 

  **Outils informatiques**

  Ce pipeline est essentiellement codé en Bash et en R. Les différentes applications utilisées sont : STAR, FastQC, Samtools, FastQ-Dump, featureCounts. Nextflow (https://nextflow.io/) est utilisé comme Workflow Management System. Nextflow peut être lancé depuis [Conda](https://conda.io). Afin d'assurer la reproductibilité des résultats, les images des différentes applications ont été générées via Docker.
    
  Les images [Dockers](https://www.docker.com/en) utilisées sont : 
     
   * [yanismadi/samtools](https://hub.docker.com/r/yanismadi/samtools) (version 1.9) 
   * [clarpalmtree/starbis](https://hub.docker.com/r/clarpalmtree/starbis) (version 2.7.10a)
   * [clarpalmtree/sra](https://hub.docker.com/r/clarpalmtree/sra) (version current) 
   * [yanismadi/subread](https://hub.docker.com/r/yanismadi/subread) (version 2.0.0)
   * [evolbioinfo/deseq2](https://hub.docker.com/r/evolbioinfo/deseq2) (version 1.28.1)
  
  
  **Données biologiques**

  Les données biologiques utilisées dans notre étude sont accessibles depuis le site du NCBI (https://www.ncbi.nlm.nih.gov/sra?term=SRA062359). Les échantillons sélectionnés sont les échantillons de type "Transcriptome sequencing" et leurs ID sont répertoriés dans le fichier [SraAccList_SRA062359.txt](https://github.com//Clarpalmtree/HackathonCYPS/blob/main/SraAccList_SRA062359.txt). Il s'agit de traiter des données de séquençage transcriptomique de tumeur, de mélanome uvéal plus précisément. Le séquençage a été réalisé avec un Illumina HiSeq 2000.
  </p>
</div>


<div align="left"><h2>Auteurs</h2></div>

<div align="justify">
  <p>
  Ce projet a été réalisé en collaboration : 
  
  * Siwar Hammami (https://github.com/siwarHm)
  * Pauline Lim (https://github.com/plim2021)
  * Yanis Madi (https://github.com/YanisMadi)
  * Clara Toussaint (https://github.com/Clarpalmtree)

  </p>
</div>
