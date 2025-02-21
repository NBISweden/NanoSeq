process MARK_DUPLICATES {

	// Directives

	debug true
	tag "${meta.id}_${meta.type}"
	label 'process_low'
	//container 'oras://community.wave.seqera.io/library/bwa-mem2_htslib_samtools:ada17f4d0d757cc3'
	// samtools, bamsormadup, fragmergepar, bammarkduplicatesopt

	// I/O & script

	input:
	tuple val(meta), path(cram)
	path reference_fasta

	output:
	tuple val(meta), path(TODO:), emit: ch_cram
	//TODO: tuple val(task.process), val('bwa-mem2'), eval('bwa-mem2 version 2>/dev/null'), topic: versions


	script:
	"""

		mkdir -p nsorted
		mkdir -p optdup
	
		ln -s ../${cram} ./nsorted/${meta.name}.cram

		samtools view -H ${cram} | grep SO:queryname > /dev/null || \\
			( rm ./nsorted/${meta.name}.cram; samtools sort -@ ${task.cpus} -O cram -m 2G -n -o ./nsorted/${meta.name}.cram ${cram} )

		samtools view --no-PG -@ ${task.cpus} -h ./nsorted/${meta.name}.cram | bamsormadup inputformat=sam level=0 blocksortverbose=0 rcsupport=1 threads=${task.cpus} fragmergepar=${task.cpus} optminpixeldif=10 | \\
			bammarkduplicatesopt verbose=0 level=0 index=0 optminpixeldif=2500 | samtools view --no-PG -@ ${task.cpus} -C -o ./optdup/${meta.name}.cram
		samtools index -@ ${task.cpus} ./optdup/${meta.name}.cram

	"""

}
