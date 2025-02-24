process ADD_READ_BUNDLES {

	// Directives

	debug true
	tag "${meta.id}_${meta.type}"
	label 'process_single'
	container 'docker://cormackinsella/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(cram)
	path reference_fasta

	output:
	//tuple val(meta), path("out/${meta.name}.cram"), path("out/${meta.name}.cram.crai"), emit: cram
	//tuple val(task.process), val('samtools'), eval('samtools version | head -n 1 | sed "s/samtools //"'), topic: versions

	script:
	"""

	mkdir -p out
	NLINES=`samtools view $cram | head -1 | grep rb: | grep rc: | grep mb: | grep mc: | wc -l` || true
	## *MODIFIED* (ao7): fixed numeric comparison operator
	## if [ \$NLINES != 0 ]; then
	if [ \$NLINES -ne 0 ]; then
	##
		bamaddreadbundles -I $cram -O ./out/${meta.name}.cram
		samtools index ./out/${meta.name}.cram
	else
	ln -s ../$cram ./out/${meta.name}.cram
	ln -s ../$index ./out/${meta.name}.cram.crai
	fi
	"""

}
