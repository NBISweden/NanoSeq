process BWA_MEM2_MAP {

	// Directives

	debug true
	tag "${meta.id}"
	label 'process_medium'
	container 'oras://community.wave.seqera.io/library/bwa-mem2:2.2.1--e269358148d98816'

	// I/O & script

	input:
	tuple val(meta), path(input_files)
	path (reference_fasta)
	path (bwa_indexes)

	output:
	//FIXME: path(), emit: 

	tuple val(task.process), val('bwa-mem2'), eval('bwa-mem2 version 2>/dev/null'), topic: versions

	script:

	// Always map fastqs
	if (meta.format == 'fastq')
		"""
		
		echo ${input_files[0]}
		echo ${input_files[1]}
		echo ${input_files[2]}
		echo ${input_files[3]}
		echo ${reference_fasta}
		
		"""
	
	// Map bams if specified in metadata
	else if (meta.format == 'bam' && mapping == true)
		"""

		echo ${input_files[0]}
		echo ${input_files[1]}
		echo ${input_files[2]}
		echo ${input_files[3]}
		echo ${reference_fasta}

		"""

	// Map crams if specified in metadata
	else if (meta.format == 'cram' && mapping == true)
		"""

		echo "do nothing yet"

		"""

}
