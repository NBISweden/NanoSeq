process COVERAGE {

	// Directives

	debug true
	tag "${meta.id}_${meta.type}"
	label 'process_medium'
	container 'docker://ghcr.io/nbisweden/nanoseq-src:latest'

	// I/O & script

	input:
	tuple val(meta), path(crams)
	path reference_fasta
	path indexes

	//publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern:"cov/*", overwrite: true

	output:
	// tuple val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal) , emit : done
	// path 'cov/1.done'
	// path 'cov/1.cov.bed.gz'
	// path 'cov/gIntervals.dat'
	// path 'cov/args.json'
	// path 'cov/nfiles'

	script:
	"""

		runNanoSeq.py --ref ${reference_fasta} --normal ${crams[2]} --duplex ${crams[0]} --threads ${task.cpus} cov -Q ${params.minimum_duplex_mapq} --exclude \"${params.coverage_exclude}\" --include \"${params.coverage_include}\" --larger ${params.coverage_larger_than}

				# To test: exclude/include - do we need these characters?
						# --exclude \"exclude\"  --include \"include\" ;

				# Test leaving out --out
						#--out $baseDir/work/temps/$meta.id/

				# Solve dependencies

				# outputs

				# Include some of their post-processing steps - it's a weird script

	"""

}
