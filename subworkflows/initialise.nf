/* 
----------------------------------------------------------------------------------------
Initialisation subworkflow
	- Verify base dependencies
	- Ensure required inputs/params are provided
	- Import samplesheet
	- Import reference genome
----------------------------------------------------------------------------------------
*/

workflow INITIALISE {

	main:

		// Warn about current issues with topic channels

			println ("WARN: package versions are not currently reported for indexing modules, due to a Nextflow bug: https://github.com/nextflow-io/nextflow/issues/5785")

		// Help message

			if (params.help) {
				println ("A very helpful message")
				error ("--help							Print this message.")
			}

		// Nextflow version

			if (!nextflow.version.matches('>=24.10.4')) {
				error ("ERROR: This workflow requires Nextflow version '>=24.10.4' or later. You are running '${nextflow.version}'. Update with 'nextflow self-update' or use the provided pixi environment.")
			}

		// Apptainer executable

			if (!"apptainer".execute().text.trim()) {
				error ("ERROR: This workflow requires the Apptainer executable available in PATH.")
			}

		// Parse samplesheet

			// Initialise empty set to store unique sample ids

				def uniqueIds = new HashSet<String>()

			// Load samplesheet

				Channel
					.fromPath("${params.samplesheet}", checkIfExists: true)
					.splitCsv(header: true)
					.map { row ->

						// Populate metadata fields from samplesheet
						def meta = [
							id: row.id,
							format: row.input_format.toLowerCase(),
							mapping:row.mapping_required.toLowerCase(),
							normalMethod: row.normal_method.toLowerCase()
						]

						// Check id names are unique
						if (!uniqueIds.add(meta.id)) {
							error ("ERROR: Duplicate sample ID found: ${meta.id}. All sample IDs must be unique.")
						}

						// Check that valid input formats were supplied
						if (!['fastq', 'bam', 'cram'].contains(meta.format)) {
							error ("ERROR: Sample '${meta.id}' has an invalid input format '${meta.format}'. Please only supply 'fastq', 'bam', or 'cram' (case insensitive).")
						}

						// Check that mapping_required is true or false, convert to boolean
						if (!['true', 'false'].contains(meta.mapping)) {
							error ("ERROR: Sample '${meta.id}' has an invalid mapping_required value '${meta.mapping}'. Please only supply 'true' or 'false' (case insensitive).")
						}
						meta.mapping = meta.mapping.toBoolean()

						// Check that matched normal sequencing method is either duplex or standard
						if (!['duplex', 'standard'].contains(meta.normalMethod)) {
							error ("ERROR: Sample '${meta.id}' has an invalid 'normal_method' value '${meta.normalMethod}'. Please only supply 'duplex' or 'standard', depending on how your matched normal samples were sequenced (case insensitive).")
						}

						// If initial checks pass, add file paths specific to input format & check their existence

						// FASTQ
						if (meta.format == 'fastq') {
							def files = [
								file(row.duplex_1, checkIfExists: true),
								file(row.duplex_2, checkIfExists: true),
								file(row.normal_1, checkIfExists: true),
								file(row.normal_2, checkIfExists: true)
							]
							return [meta, files]

						// BAM (if remapping is required, don't add indexes to the map, keep value position consistent)
						} else if (meta.format == 'bam') {
							if (meta.mapping) {
								def files = [
									file(row.duplex_1, checkIfExists: true),
									file(row.normal_1, checkIfExists: true)
								]
								return [meta, files]
							} else {
								def files = [
									file(row.duplex_1, checkIfExists: true),
									file(row.normal_1, checkIfExists: true),
									row.duplex_2 ? file(row.duplex_2, checkIfExists: true) : error ("ERROR: Sample '${meta.id}' needs an index file when mapping is not required."),
									row.normal_2 ? file(row.normal_2, checkIfExists: true) : error ("ERROR: Sample '${meta.id}' needs an index file when mapping is not required.")
								]
								return [meta, files]
							}

						// CRAM (if remapping is required, don't add indexes to the map, keep value position consistent)
						} else if (meta.format == 'cram') {
							if (meta.mapping) {
								def files = [
									file(row.duplex_1, checkIfExists: true),
									file(row.normal_1, checkIfExists: true)
								]
								return [meta, files]
							} else {
								def files = [
									file(row.duplex_1, checkIfExists: true),
									file(row.normal_1, checkIfExists: true),
									row.duplex_2 ? file(row.duplex_2, checkIfExists: true) : error ("ERROR: Sample '${meta.id}' needs an index file when mapping is not required."),
									row.normal_2 ? file(row.normal_2, checkIfExists: true) : error ("ERROR: Sample '${meta.id}' needs an index file when mapping is not required.")
								]
								return [meta, files]
							}
						}
					}
					.set { ch_samplesheet }

		// Parse reference genome & optional indexes (workflow creates them if missing)

			Channel
				.fromPath("${params.fasta}", checkIfExists: true)
				.map { file ->
					def extension = file.extension
					if (['fa', 'fna', 'fasta'].contains(extension)) {
						return file
					} else {
						error ("ERROR: Reference file '${file}' does not have a '.fa', '.fna', or '.fasta' extension, please provide a FASTA file.")
					}
				}
				.set { ch_reference }

	emit:

		// Emit channels to main workflow

			ch_samplesheet = ch_samplesheet
			ch_reference = ch_reference

}
