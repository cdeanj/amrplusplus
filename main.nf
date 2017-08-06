#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

if( !nextflow.version.matches('0.25+') ) {
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    exit 1
}
if( params.host_index ) { 
    host_index = Channel.fromPath(params.host_index).toSortedList() 
    if( !host_index.exists() ) exit 1, "Host index files could not be found: ${params.host_index}"    
}
if( params.amr_index ) { 
    amr_index = Channel.fromPath(params.amr_index).toSortedList() 
    if( !amr_index.exists() ) exit 1, "AMR index files could not be found: ${params.amr_index}"
}
if( params.host ) { 
    host = file(params.host) 
    if( !host.exists() ) exit 1, "Host genome file could not be found: ${params.host}"
}
if( params.amr ) { 
    amr = file(params.amr) 
    if( !amr.exists() ) exit 1, "AMR database file could not be found: ${params.amr}" 
}
if( params.adapters ) { 
    adapters = file(params.adapters) 
    if( !adapters.exists() ) exit 1, "Adapter file could not be found: ${params.adapters}" 
}
if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) exit 1, "Annotation file could not be found: ${params.annotation}"
}

threads = params.threads
threshold = params.threshold

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .into { reads }

process RunQC {
    tag { sample_id }

    publishDir "${params.output}/Trimmomatic", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf("P.fastq") > 0) "Paired/$filename"
            else if(filename.indexOf("U.fastq") > 0) "Unpaired/$filename"
            else if(filename.indexOf(".log") > 0) "Log/$filename"
            else {}
        }
	
    input:
        set sample_id, file(forward), file(reverse) from reads

    output:
        set sample_id, file("${sample_id}.1P.fastq"), file("${sample_id}.2P.fastq") into (paired_fastq)
        set sample_id, file("${sample_id}.1U.fastq"), file("${sample_id}.2U.fastq") into (unpaired_fastq)
        set sample_id, file("${sample_id}.trimmomatic.stats.log") into (trimmomatic_logs)

    """
    java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar \
      PE \
      -threads ${threads} \
      $forward $reverse -baseout ${sample_id} \
      ILLUMINACLIP:${adapters}:2:30:10:3:TRUE \
      LEADING:${leading} \
      TRAILING:${trailing} \
      SLIDINGWINDOW:${slidingwindow} \
      MINLEN:${minlen} \
      2> ${sample_id}.trimmomatic.stats.log

    mv ${sample_id}_1P ${sample_id}.1P.fastq
    mv ${sample_id}_2P ${sample_id}.2P.fastq
    mv ${sample_id}_1U ${sample_id}.1U.fastq
    mv ${sample_id}_2U ${sample_id}.2U.fastq
    """
}

if( !params.host_index ) {
    process BuildHostIndex {
        tag { host.baseName }

        input:
            file(host)

        output:
            file '*' into host_index

        """
        bwa index ${host}
        """
    }
}

process AlignReadsToHost {
    tag { sample_id }
        
    publishDir "${params.output}/Host", mode: "copy"
        
    input:
        set sample_id, file(forward), file(reverse) from paired_fastq
        file index from host_index.first()
        file host
            
    output:
        set sample_id, file("${sample_id}.host.sam") into host_sam
            
    """ 
    bwa mem ${host} ${forward} ${reverse} -t ${threads} > ${sample_id}.host.sam
    """ 
}

process RemoveHostDNA {
    tag { sample_id }

    publishDir "${params.output}/Host", mode: "copy"

    input:
        set sample_id, file(sam) from host_sam

    output:
        set sample_id, file("${sample_id}.host.sorted.removed.bam") into non_host_bam

    """
    samtools view -bS ${sam} | samtools sort -@ ${threads} -o ${sample_id}.host.sorted.bam
    samtools view -h -f 4 -b ${sample_id}.host.sorted.bam -o ${sample_id}.host.sorted.removed.bam
    """
}

process BAMToFASTQ {
    tag { sample_id }

    publishDir "${params.output}/Host", mode: "copy"

    input:
        set sample_id, file(bam) from non_host_bam

    output:
        set sample_id, file("${sample_id}.non.host.R1.fastq"), file("${sample_id}.non.host.R2.fastq") into non_host_fastq

    """
    bedtools  \
       bamtofastq \
      -i ${bam} \
      -fq ${sample_id}.non.host.R1.fastq \
      -fq2 ${sample_id}.non.host.R2.fastq
    """
}

if( !params.amr_index ) {
    process BuildAMRIndex {
        tag { amr.baseName }

        input:
            file(amr)

        output:
            file '*' into amr_index

        """
        bwa index ${amr}
        """
    }
}

process AlignToAMR {
     tag { sample_id }

     publishDir "${params.output}/AMR", mode: "copy"

     input:
         set sample_id, file(forward), file(reverse) from non_host_fastq
         file index from amr_index.first()
         file amr

     output:
         set sample_id, file("${sample_id}.amr.alignment.sam") into amr_sam

     """
     bwa mem ${amr} ${forward} ${reverse} -t ${threads} > ${sample_id}.amr.alignment.sam
     """
}

process AnalyzeResistome {
    tag { sample_id }

    publishDir "${params.output}/AnalyzeResistome", mode: "copy"

    input:
        set sample_id, file(sam) from amr_sam
        file annotation
        file amr

    output:
        set sample_id, file("*.tsv") into resistome
    
    """
    resistome \
      -ref_fp ${amr} \
      -annot_fp ${annotation} \
      -sam_fp ${sam} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -class_fp ${sample_id}.class.tsv \
      -mech_fp ${sample_id}.mechanism.tsv \
      -t ${threshold}
    """
}
