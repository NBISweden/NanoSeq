/* 
----------------------------------------------------------------------------------------
Verify base dependencies & check required inputs
----------------------------------------------------------------------------------------
*/

// Verification subworkflow

workflow VERIFY {

	// Help message

		if (params.help) {
			println ("A very helpful message")
			error ("--help							Print this message.")
		}

	// Nextflow version

		if (!nextflow.version.matches('>=24.02.0')) {
			error ("ERROR: This workflow requires Nextflow version '24.02.0-edge' or later. You are running '${nextflow.version}'. Update with 'nextflow self-update'.")
		}

	// Apptainer executable

		if (!"apptainer".execute().text.trim()) {
			error ("ERROR: This workflow requires the Apptainer executable available in PATH.")
		}

	// Check for input samplesheet

		if (!file("${params.samplesheet}").exists()) {
			error ("ERROR: Input file '${params.samplesheet}' was not found (please see docs).")
		}

	// Process samplesheet

		// Initialise empty set to store unique sample ids
		def uniqueIds = new HashSet<String>()

		// Load samplesheet
		def ch_samplesheet = Channel
			.fromPath("${params.samplesheet}")
			.splitCsv(header: true)
			.map { row ->
				// Populate id, format, and mapping metadata fields from samplesheet
				def meta = [id: row.id, format: row.input_format.toLowerCase(), mapping:row.mapping_required.toLowerCase()]

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

				// If checks pass, add metadata fields appropriate to the input format, check files exist, and check the matched normal sequencing method is recognised

				// FASTQ
				if (meta.format == 'fastq') {
					meta.duplex1 = file(row.duplex_1, checkIfExists: true)
					meta.duplex2 = file(row.duplex_2, checkIfExists: true)
					meta.normal1 = file(row.normal_1, checkIfExists: true)
					meta.normal2 = file(row.normal_2, checkIfExists: true)
					meta.normalMethod = row.normal_method.toLowerCase()
					if (!['duplex', 'standard'].contains(meta.normalMethod)) {
						error ("ERROR: Sample '${meta.id}' has an invalid 'normal_method' value '${meta.normalMethod}'. Please only supply 'duplex' or 'standard', depending on how your matched normal samples were sequenced (case insensitive).")
					}
					return meta

				// BAM (if remapping is required, allow the index to be absent)
				} else if (meta.format == 'bam') {
					meta.duplexBam = file(row.duplex_1, checkIfExists: true)
					meta.normalBam = file(row.normal_1, checkIfExists: true)
					if (!meta.mapping) {
						meta.duplexIndex = row.duplex_2 ? file(row.duplex_2, checkIfExists: true) : error ("ERROR: Sample '${meta.id}' needs a duplex index file when mapping is required.")
						meta.normalIndex = row.normal_2 ? file(row.normal_2, checkIfExists: true) : error ("ERROR: Sample '${meta.id}' needs a matched normal index file when mapping is required.")
					} else {
						meta.duplexIndex = row.duplex_2 ? file(row.duplex_2, checkIfExists: false) : null
						meta.normalIndex = row.normal_2 ? file(row.normal_2, checkIfExists: false) : null
					}
					meta.normalMethod = row.normal_method.toLowerCase()
					if (!['duplex', 'standard'].contains(meta.normalMethod)) {
						error ("ERROR: Sample '${meta.id}' has an invalid 'normal_method' value '${meta.normalMethod}'. Please only supply 'duplex' or 'standard', depending on how your matched normal samples were sequenced (case insensitive).")
					}
					return meta

				// CRAM (if remapping is required, allow the index to be absent)
				} else if (meta.format == 'cram') {
					meta.duplexCram = file(row.duplex_1, checkIfExists: true)
					meta.normalCram = file(row.normal_1, checkIfExists: true)
					if (!meta.mapping) {
						meta.duplexIndex = row.duplex_2 ? file(row.duplex_2, checkIfExists: true) : error ("ERROR: Sample '${meta.id}' needs a duplex index file when mapping is required.")
						meta.normalIndex = row.normal_2 ? file(row.normal_2, checkIfExists: true) : error ("ERROR: Sample '${meta.id}' needs a matched normal index file when mapping is required.")
					} else {
						meta.duplexIndex = row.duplex_2 ? file(row.duplex_2, checkIfExists: false) : null
						meta.normalIndex = row.normal_2 ? file(row.normal_2, checkIfExists: false) : null
					}
					meta.normalMethod = row.normal_method.toLowerCase()
					if (!['duplex', 'standard'].contains(meta.normalMethod)) {
						error ("ERROR: Sample '${meta.id}' has an invalid 'normal_method' value '${meta.normalMethod}'. Please only supply 'duplex' or 'standard', depending on how your matched normal samples were sequenced (case insensitive).")
					}
					return meta
				}
			}

	// 

}
