#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

if( params.index ) {
    index = Channel.fromPath(params.index).toSortedList()
}

host = file(params.host)

threads = params.threads

adapters = file(params.adapters)

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
