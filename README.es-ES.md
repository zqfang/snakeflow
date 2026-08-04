# snakeflow
Mis flujos de trabajo de NGS basados en snakemake

**los flujos de trabajo incluyen**:
- RNA-seq (Salmon, hisat2)
- ChIP-seq (MACS2)
- ATAC-seq (MACS2)
- CITE-seq (Antibody captured, 10X genomics)
- Llamada de SNV germinal (GATK, BCFtools)
- Llamada de Variantes Estructurales germinales
    - lecturas cortas: Speedseq + svtools
    - lecturas largas:
        - Pacbio: ngmlr/minimap2 + sniffiles


## Dependencia

* General
  - samtools, deeptools, bedtools
  - fastqc, rseqc, multiqc, fastp
  - graphviz

* python 3
  - numpy
  - pandas
  - snakemake
  - matplotlib
  - seaborn
  - gseapy
  - macs2
  - rseqc

* R
  - DESeq2
  - tximport
  - readr
  - pheatmap
  - ggplot2
  - ggrepel
  - clusterProfiler
  - ChIPSeeker
  - EnsDb.Hsapiens.v86


* Llamada de variantes
  - GATK (> 4.0)
  - BCFtools
  - Speedseq + svtools
  - minimap2
  - ngmlr
  - sniffiles

* RNA-seq
  - hisat2, salmon
  - rMATS-turbo, rmats2sashimiplot

* Genómica de célula única
 - cellranger

## Instalación

### Configurar el entorno de ejecución. Este archivo de configuración creará un entorno python 3.x.

    bash snakeflow-enviroment-setup.sh
    
## uso
    
    # Paso 1: activar snakemake
    source activate snakeflow

    # Paso 2: clonar este repositorio
    
    # Paso 3: copiar todos tus archivos fastq al directorio fastq
    find . -name "*fastq.gz" | while read id; do cp $id fastq/; done;
    
    # Paso 4: modificar config.yml con tus propios parámetros
    # Nota: coloca config.yml en el mismo directorio que tu snakefile.
    vim  config.yml

    # Paso 5: ejecutar snakemake con la opción -np. Esto prueba que tus ``comandos`` se ejecutan sin errores.
    snakemake -s salmon-tximport-deseq2-v0.2.snakefile -np

    # Paso 6: exportar los diagramas del flujo de trabajo
    snakemake -s salmon-tximport-deseq2-v0.2.snakefile --dag | dot -Tpdf > dag.pdf

    # Paso 7: o usando el entorno predeterminado de snakemake que creaste anteriormente.
    snakemake -s salmon-tximport-deseq2-v0.1.snakefile -p -j 8
