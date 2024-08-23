# metaviroma
Identifica fagos y profagos de DNA a partir de archivos FASTQ de muestras metagenomicas secuenciadas por Shotgun, a través de QA de las secuencias, ensamble, asignación taxonómica y finalmente identificación de fagos y profagos.

#checkV identifica "virus totales" y profagos
#PhaMer identifica bacteriofagos (OUT.csv)
#PHASTER anota fagos y profagos

#Quality Control, trimm and remove raw reads
fastp --overrepresentation_analysis --correction --thread 8 --html file.html --json file.json -i file_R1_001.fastq.gz -I file_R2_001.fastq.gz -o file_1.fq.gz -O file_2.fq.gz

#Ensamble del metagenoma
megahit -t 8 --kmin-1pass -1 489_1.fastq.gz,490_1.fastq.gz,491_1.fastq.gz,493_1.fastq.gz,494_1.fastq.gz -2 489_2.fastq.gz,490_2.fastq.gz,491_2.fastq.gz,493_2.fastq.gz,494_2.fastq.gz -o out

#Quality Control, assembly
~/programas/quast-5.0.2/quast.py -t 8 final.contigs.fa


########################################################################
#Identificación de genomas virales TOTALES y clasifica Profagos.

#CheckV, identifica genomas virales cerrados, estima la integridad de los fragmentos genómicos y elimina las regiones huésped flanqueantes de los provirus integrados.

#Using a single command to run the full pipeline
checkv end_to_end out/megahit_results/final.contigs.fa output_checkv -d ~/db/checkv-db-v1.1/ -t 4

#There are two ways to run CheckV:
#Using individual commands for each step in the pipeline in the following order:

#checkv contamination input_file.fna output_directory -t 16
#checkv completeness input_file.fna output_directory -t 16
#checkv complete_genomes input_file.fna output_directory
#checkv quality_summary input_file.fna output_directory

#######################
#INTERPRETAR RESULTADOS
#quality_summary.tsv, pestaña 'checkv_quality' DETERMINA SI CADA CONTIG, ES O NO UN VIRUS, ADEMÁS LA PESTAÑA 'provirus' identifica si es o no un provirus
#'Not-determined' This contig also has no viral genes identified, so there's a chance it may not even be a virus.
#'Low-quality' since its completeness is <50%. This is estimate is based on the 'AAI' method. Note that only either high- or medium-confidence estimates are reported in the quality_summary.tsv file. You can see 'completeness.tsv' for more details. This contig had a DTR, but it was flagged for some reason (see complete_genomes.tsv for details)
#'Medium-quality' since its completeness is estimated to be 80%, which is based on the 'HMM' method. This means that it was too novel to estimate completeness based on AAI, but shared an HMM with CheckV reference genomes. Note that this value represents a lower bound (meaning the true completeness may be higher but not lower than this value). Note that this contig is also classified as a provirus.
#High-quality based on a completness of >90%. However, note that value of 'kmer_freq' is 1.7. This indicates that the viral genome is represented multiple times in the contig. These cases are quite rare, but something to watch out for.
#Complete based on the presence of a direct terminal repeat (DTR) and has 100% completeness based on the AAI method. This sequence can condifently treated as a complete genome.


########################################################################
#PhaMer. Identifica bacteriofagos de datos metagenòmicos

#Input ensamble de MEGAHIT (final.contigs.fa) y con los virus de checkV (viruses.fna).
python PhaMer_preprocessing.py --threads 8 --contigs final.contigs.fa
python PhaMer.py --out OUTPUT.csv --threads 8
#Resultado: En el archivo OUTPUT.csv identifica contigs que contienen fagos, columna Pred y columna Score > 0.5.

########################################################################
#PHASTER (PHAge Search Tool Enhanced Release) es una actualización significativa del popular servidor web PHAST para la rápida identificación y anotación de secuencias de profagos dentro de genomas y plásmidos bacterianos.

#PHASTER (phaster_scripts.py) es necesario estar conectado a internet, ya que sube las secuencias a la plataforma https://phaster.ca/. La secuencia de ADN debe ajustarse a la codificación IUPAC. Un archivo metagenómico debe contener contigs, cada uno con su propio encabezado. A partir de esta entrada, sólo se procesarán los contigs de longitud >=2000. Si la opción metagenómica no está marcada, las 10 primeras secuencias se procesarán individualmente.

#EN PHASTER se deben de correr los resultados de checkV (virus totales y/o Profagos o solo de los Profagos para poder anotarlos).

#Someter un solo genoma
python phaster.py --fasta genome.fasta

#Someter varios contigs
python phaster.py --contigs --fasta contigs.fasta
#Saber el status de la búsuqeda, en el archivo "phaster_jobs.tsv" se encuentra la accesión se copia y se pega al final de https://phaster.ca/submissions/'accesion' o en https://phaster.ca/batches/'clave' si ya terminó
python phaster.py --get-status


########################################################################
#KRAKEN2, Taxonomic Assignment, solo de virus (RECOMIENDO HACERLO CON LOS CONTIGS EXTRAIDOS DE LOS RESULTADOS DE checkV, virus putativos y/o profagos) y los resultados de PhaMer y PHASTER

#construir base de datos
#kraken2-build --standard --threads 6 --db $DBNAME
#Clasification
#kraken2 --db $DBNAME --threads 6 seqs.fa --output 
#Custom database
#kraken2-build --download-taxonomy --db $DBNAME
#ejemplo: kraken2-build --download-library bacteria --db $DBNAME

#Asignación Taxonomica de reads metagenomicos
#kraken2 --db ~/programas/kraken2-2.0.8-beta/k2_viral_20210517/ --threads 4 --paired --fastq-input 1F.fq.gz 1R.fq.gz --output xxx_reads.kraken --report xxx_reads.report

#Asignación Taxonomica de contigs de un MAG
kraken2 --db ~/programas/kraken2-2.0.8-beta/k2_viral_20210517/ --threads 8 --output virus_MAG.kraken --report virus_MAG.report out/megahit_results/final.contigs.fa
cut -f2,3 virus_MAG.kraken > virus.krona.input
ktImportTaxonomy virus.krona.input -o virus.krona.html

#Revisar resultados de checkv
kraken2 --db ~/programas/kraken2-2.0.8-beta/k2_viral_20210517/ --threads 8 --output virus_checkv.kraken --report virus_checkv.report output_checkv/viruses.fna
kraken-biom virus_checkv.report --fmt json -o virus_checkv.biom

#Phagos
kraken2 --db ~/programas/kraken2-2.0.8-beta/k2_viral_20210517/ --threads 8 --output virus_phagos_checkv.kraken --report virus_phagos_checkv.report PhaMer/fagos.fasta
kraken-biom virus_phagos_checkv.report --fmt json -o virus_phagos_checkv.biom

#Prophagos
kraken2 --db ~/programas/kraken2-2.0.8-beta/k2_viral_20210517/ --threads 8 --output virus_prophagos_checkv.kraken --report virus_prophagos_checkv.report output_checkv/proviruses.fna
kraken-biom virus_prophagos_checkv.report --fmt json -o virus_prophagos_checkv.biom


#El resultado se sube al programa MEGAN6 y se graficó el árbol de abundancia y el word cloud




########################################################################

#Asignación taxonómica del hospedero para Fagos y Profagos, base de datos genoes/HUMAN_MICROBIOM/Bacteria/all.fna.tar.gz del NCBI (asignacion_tax_bacterias)

makeblastdb -in all_final/all.fna -dbtype nucl -out ../all -parse_seqids

blastn -query PhaMer/fagos.fasta -db all -max_target_seqs 1 -outfmt "6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore" -html -evalue 1e-23 -num_threads 8 > fagos.csv

blastn -query output_checkv/proviruses.fna -db all -max_target_seqs 1 -outfmt "6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore" -html -evalue 1e-23 -num_threads 8 > profagos.csv

#Copiar columna B de las accesiones de los genomas y crear lista para extraer headers
blastdbcmd -db all -entry_batch list.txt -out lista_seq_fago.fasta 
blastdbcmd -db all -entry_batch list_pro.txt -out lista_seq_pro.fasta

#Extraer headers de secuencias, borrar > copiar y pegar en la columna final de fagos.csv y profagos.csv (pegar despues la asignación de fagos y profagos en la última columna)
grep '>' lista_seq_fago.fasta > taxa_fago.txt
grep '>' lista_seq_pro.fasta > taxa_pro.txt

