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

			println ("WARN: package versions are not currently reported for the indexing module, due to a Nextflow bug: https://github.com/nextflow-io/nextflow/issues/5785")

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

		// Optional file existence checks

			if (params.bed_to_exclude) {
				if (!file(params.bed_to_exclude).exists()) {
					error ("ERROR: The BED file to exclude does not exist ('${params.bed_to_exclude}').")
				}
			}

			if (params.snp_bed) {
				if (!file(params.snp_bed).exists()) {
					error ("ERROR: The SNP BED file does not exist ('${params.snp_bed}').")
				}
				def snp_bed_index = params.snp_bed + '.tbi'
				if (!file(snp_bed_index).exists()) {
					error ("ERROR: The SNP BED index file does not exist ('${snp_bed_index}'). Ensure you have indexed the sorted BED with tabix.")
				}
			}

			if (params.noise_bed) {
				if (!file(params.noise_bed).exists()) {
					error ("ERROR: The noise BED file does not exist ('${params.noise_bed}').")
				}
				def noise_bed_index = params.noise_bed + '.tbi'
				if (!file(noise_bed_index).exists()) {
					error ("ERROR: The noise BED index file does not exist ('${noise_bed_index}'). Ensure you have indexed the sorted BED with tabix.")
				}
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
							normalMethod: row.normal_method.toLowerCase()
						]
						// Check id names are unique
						if (!uniqueIds.add(meta.id)) {
							error ("ERROR: Duplicate sample ID found: ${meta.id}. All sample IDs must be unique.")
						}
						// Check that matched normal sequencing method is either duplex or standard
						if (!['duplex', 'standard'].contains(meta.normalMethod)) {
							error ("ERROR: Sample '${meta.id}' has an invalid 'normal_method' value '${meta.normalMethod}'. Please only supply 'duplex' or 'standard', depending on how your matched normal samples were sequenced (case insensitive).")
						}
						// Add FASTQ file paths
						def files = [
							file(row.duplex_1, checkIfExists: true),
							file(row.duplex_2, checkIfExists: true),
							file(row.normal_1, checkIfExists: true),
							file(row.normal_2, checkIfExists: true)
						]
						return [meta, files]
					}
					.set { ch_samplesheet }

		// Parse reference genome

			Channel
				.fromPath("${params.fasta}", checkIfExists: true)
				.map { file ->
					// Allow range of FASTA extensions
					def extension = file.extension
					if (!['fa', 'fna', 'fasta'].contains(extension)) {
						error ("ERROR: Reference file '${file}' does not have a '.fa', '.fna', or '.fasta' extension, please provide a valid FASTA file.")
					}
					// Ensure FASTA headers don't contain invalid characters
					file.withReader { reader ->
						reader.eachLine { line ->
							if (line.startsWith('>')) {
								def header = line.substring(1).trim()
								if (header.contains(':') || header.contains(',')) {
									error ("ERROR: Reference FASTA header '${header}' contains invalid characters (':' or ','). These are incompatible with a downstream script.")
								}
							}
						}
					}
					return file
				}
				.set { ch_reference }

	emit:

		// Emit channels to main workflow

			ch_samplesheet
			ch_reference

}
