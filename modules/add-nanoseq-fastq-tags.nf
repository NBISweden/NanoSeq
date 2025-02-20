process ADD_NANOSEQ_FASTQ_TAGS {

	// Directives

	debug true
	tag "${meta.id}_${meta.type}"
	label 'process_single'
	container 'oras://community.wave.seqera.io/library/python:3.13.1--9856f872fdeac74e'

	// I/O & script

	input:
	tuple val(meta), path(reads)

	output:
	tuple val(meta), path("*.fastq.gz"), emit: ch_tagged_fastqs
	tuple val(task.process), val('python'), eval('python --version | sed "s/.* //"'), topic: versions
	tuple val(task.process), val('runNanoSeq.py'), eval('runNanoSeq.py -v'), topic: versions

	script:

		"""

		# Get read length

			length=`zcat ${reads[0]} | head -2 | tail -1 | awk '{ print length }'`

		# Extract tags

			extract_tags.py -a ${reads[0]} -b ${reads[1]} -c ${meta.id}_${meta.type}_R1.fastq.gz -d ${meta.id}_${meta.type}_R2.fastq.gz -m ${params.fastq_tags_m} -s ${params.fastq_tags_s} -l \$length

		"""

}
