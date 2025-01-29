/* 
----------------------------------------------------------------------------------------
Main workflow definition
----------------------------------------------------------------------------------------
*/

// Feature flags

	nextflow.preview.topic = true

// Import modules

	include { BWA_MEM2_INDEX } from '../modules/bwa-index.nf'
	include { SAMTOOLS_INDEX } from '../modules/samtools-index.nf'
	include { TEST } from '../modules/test.nf'

// Parse reference genome & indexes (create latter if missing)

	// Define FASTA file
	def fasta_extensions = ['fa', 'fna', 'fasta']
	def reference_fasta = file("${params.referencePath}/${params.fasta}")

	if (reference_fasta.exists()) {
		def file_extension = reference_fasta.extension
		if (fasta_extensions.contains(file_extension)) {
			ch_reference_fasta = Channel.fromPath(reference_fasta)
		} else {
			error ("ERROR: Reference file '${reference_fasta}' does not have a '.fa', '.fna', or '.fasta' extension, please provide a FASTA file.")
		}
	} else {
		error ("ERROR: Reference file not found at '${reference_fasta}'. Please check you've set the parameters '--reference_path' and '--fasta' correctly.")
	}

	// Function to check if BWA indexes exist, returns true if all are found
	indexExtensions = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
	def bwaIndexesExist(fasta) {
		return indexExtensions.every { ext -> file("${fasta}${ext}").exists() }
	}

	// Run function. If indexes exist, add to a channel. If no, set channel variable to null for later creation.
	ch_bwa_indexes = bwaIndexesExist(reference_fasta) ? Channel.value(indexExtensions.collect { ext -> file("${reference_fasta}${ext}") }) : null

	// Function to check if samtools dict exists, returns true if found
	def samtoolsIndexExists(fasta) {
		return file("${fasta}.dict").exists()
	}

	// Run function. If index exists, add to a channel. If no, set channel variable to null for later creation.
	ch_samtools_index = samtoolsIndexExists(reference_fasta) ? Channel.fromPath("${reference_fasta}.dict") : null








//TODO: next the DICT index. check how many outputs are created here also, may need to do similar
	//ch_dict = file("${reference_fasta}.dict").exists() ? Channel.fromPath("${reference_fasta}.dict") : null




















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

// Workflow definition

	workflow NANOSEQ_FORK {

		// Create BWA index files if they don't exist
			if (ch_bwa_indexes == null) {
				BWA_MEM2_INDEX(ch_reference_fasta)
				ch_bwa_indexes = BWA_MEM2_INDEX.out.ch_bwa_indexes
			}

		// Create dict index if it doesn't exist
			if (ch_samtools_index == null) {
				SAMTOOLS_INDEX(ch_reference_fasta)
				ch_samtools_index = SAMTOOLS_INDEX.out.ch_samtools_index
			}

		// TEST PROCESS
			TEST(ch_samtools_index)

		// Report package versions

			Channel.topic('versions')
				.map { process, tool, version ->
					return [process: process, tool: tool, version: version]
				}
				.unique()
				.collect()
				.map { it.join('\n') } // Back to a single string
				.collectFile(name: 'package_versions.txt', newLine: true, storeDir: 'output/package_versions')

	}

