#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

if( !nextflow.version.matches('0.25+') ) {
    return nextflow_version_error()
}
if( params.host_index ) { 
    host_index = Channel.fromPath(params.host_index).toSortedList() 
    if( !host_index.exists() ) return index_error(host_index)
}
if( params.host ) { 
    host = file(params.host) 
    if( !host.exists() ) return host_error(host)
}
if( params.amr ) { 
    amr = file(params.amr) 
    if( !amr.exists() ) return amr_error(amr)
}
if( params.adapters ) { 
    adapters = file(params.adapters) 
    if( !adapters.exists() ) return adapter_error(adapters)
}
if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) return annotation_error(annotation)
}

threads = params.threads

threshold = params.threshold

min = params.min
max = params.max
skip = params.skip
samples = params.samples

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

    publishDir "${params.output}/RunQC", mode: 'copy',
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
        publishDir "${params.output}/BuildHostIndex", mode: "copy"

        tag { host.baseName }

        input:
            file(host)

        output:
            file '*' into (host_index)

        """
        bwa index ${host}
        """
    }
}

process AlignReadsToHost {
    tag { sample_id }
        
    publishDir "${params.output}/AlignReadsToHost", mode: "copy"
        
    input:
        set sample_id, file(forward), file(reverse) from paired_fastq
        file index from host_index.first()
        file host
            
    output:
        set sample_id, file("${sample_id}.host.sam") into (host_sam)
            
    """ 
    bwa mem ${host} ${forward} ${reverse} -t ${threads} > ${sample_id}.host.sam
    """ 
}

process RemoveHostDNA {
    tag { sample_id }

    publishDir "${params.output}/RemoveHostDNA", mode: "copy"

    input:
        set sample_id, file(sam) from host_sam

    output:
        set sample_id, file("${sample_id}.host.sorted.removed.bam") into (non_host_bam)

    """
    samtools view -bS ${sam} | samtools sort -@ ${threads} -o ${sample_id}.host.sorted.bam
    samtools view -h -f 4 -b ${sample_id}.host.sorted.bam -o ${sample_id}.host.sorted.removed.bam
    """
}

process BAMToFASTQ {
    tag { sample_id }

    publishDir "${params.output}/BAMToFASTQ", mode: "copy"

    input:
        set sample_id, file(bam) from non_host_bam

    output:
        set sample_id, file("${sample_id}.non.host.R1.fastq"), file("${sample_id}.non.host.R2.fastq") into (non_host_fastq)

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
            file '*' into (amr_index)

        """
        bwa index ${amr}
        """
    }
}

process AlignToAMR {
     tag { sample_id }

     publishDir "${params.output}/AlignToAMR", mode: "copy"

     input:
         set sample_id, file(forward), file(reverse) from non_host_fastq
         file index from amr_index.first()
         file amr

     output:
         set sample_id, file("${sample_id}.amr.alignment.sam") into (resistome_sam, rarefaction_sam, snp_sam)

     """
     bwa mem ${amr} ${forward} ${reverse} -t ${threads} > ${sample_id}.amr.alignment.sam
     """
}

process RunResistome {
    tag { sample_id }

    publishDir "${params.output}/RunResistome", mode: "copy"

    input:
        set sample_id, file(sam) from resistome_sam
        file annotation
        file amr

    output:
        set sample_id, file("*.tsv") into (resistome)
    
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

process RunRarefaction {
    tag { sample_id }
   
    publishDir "${params.output}/RunRarefaction", mode: "copy"

    input:
        set sample_id, file(sam) from rarefaction_sam
        file annotation
        file amr

    output:
        set sample_id, file("*.tsv") into (rarefaction)

    """
    rarefaction \
      -ref_fp ${amr} \
      -sam_fp ${sam} \
      -annot_fp ${annotation} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -class_fp ${sample_id}.class.tsv \
      -mech_fp ${sample_id}.mech.tsv \
      -min ${min} \
      -max ${max} \
      -skip ${skip} \
      -samples ${samples} \
      -t ${threshold}
    """
}

process RunSNPFinder {
    tag { sample_id }

    publishDir "${params.output}/RunSNPFinder", mode: "copy"

    input:
        set sample_id, file(sam) from snp_sam
        file amr

    output:
        set sample_id, file("*.tsv") into (snp)

    """
    snpfinder \
      -amr_fp ${amr} \
      -sampe ${sam} \
      -out_fp ${sample_id}.tsv
    """
}

def nextflow_version_error() {
    println ""
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    println ""
    return 1
}

def adapter_error(def input) {
    println ""
    println "[params.adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def amr_error(def input) {
    println ""
    println "[params.amr] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def annotation_error(def input) {
    println ""
    println "[params.annotation] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastq_error(def input) {
    println ""
    println "[params.reads] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def host_error(def input) {
    println ""
    println "[params.host] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def index_error(def input) {
    println ""
    println "[params.host_index] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}
