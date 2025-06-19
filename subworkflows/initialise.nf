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
				println (" ")
				println ("NanoSeq (NBIS version)")
				println (" ")
				println ("Command line parameters [default option within square brackets]:")
				println (" ")
				println ("General settings:")
				println (" ")
				println ("--account [${params.account}] 						Used only for the SLURM executor, provide an allocation name for hours billing.")
				println ("--email_report [${params.email_report}] 						Receive an email on pipeline completion.")
				println ("--email [${params.email}] 							Email address for report.")
				println (" ")
				println ("Input settings:")
				println (" ")
				println ("--samplesheet [${params.samplesheet}] 	Path to the samplesheet CSV file.")
				println ("--fasta [${params.fasta}] 	Path to the reference genome FASTA.")
				println (" ")
				println ("Extract tags:")
				println (" ")
				println ("--fastq_tags_m [${params.fastq_tags_m}] 						Sonicated/HpyCH4V/AluI library = 3.")
				println ("--fastq_tags_s [${params.fastq_tags_s}] 						Sonicated library = 2. HpyCH4V/AluI library = 4.")
				println (" ")
				println ("Deduplication:")
				println (" ")
				println ("--min_reads_in_bundle [${params.min_reads_in_bundle}] 					Minimum number of reads in a bundle.")
				println (" ")
				println ("Coverage:")
				println (" ")
				println ("--min_duplex_mapq [${params.min_duplex_mapq}] 						Minimum mapping quality to include a duplex read.")
				println ("--contig_exclude [${params.contig_exclude}] 						List of contigs to exclude, comma separated, %% acts as a wild card (e.g., MT,GL%%,NC_%%,hs37d5).")
				println ("--contig_include [${params.contig_include}] 						List of contigs to include, comma separated, %% acts as a wild card.")
				println ("--contig_larger_than [${params.contig_larger_than}] 					Only include contigs larger than this size.")
				println (" ")
				println ("Partition jobs:")
				println (" ")
				println ("--jobs [${params.jobs}] 							Number of jobs to partition each sample into.")
				println ("--coverage_exclude [${params.coverage_exclude}] 						Exclude regions with coverage values higher than this (0 = no exclusion).")
				println ("--bed_to_exclude [${params.bed_to_exclude}] 						BED file with regions to exclude, empty string if no exclusion.")
				println (" ")
				println ("DSA:")
				println (" ")
				println ("--snp_bed [${params.snp_bed}] 							Full path to BED of sorted SNPs, empty string if none.")
				println ("--noise_bed [${params.noise_bed}] 							Full path to BED of sorted noise regions to be masked, empty string if none.")
				println ("--min_duplex_depth [${params.min_duplex_depth}] 						Minimum duplex depth.")
				println ("--min_normal_base_quality [${params.min_normal_base_quality}] 					Minimum base quality for normal reads.")
				println (" ")
				println ("Variant calling:")
				println (" ")
				println ("--min_as_xs [${params.min_as_xs}] 						Minimum AS-XS value for variant calling.")
				println ("--min_matched_normal_reads [${params.min_matched_normal_reads}] 					Minimum matched normal reads per strand.")
				println ("--clips_fraction [${params.clips_fraction}] 					Fraction of clips.")
				println ("--consensus_fraction [${params.consensus_fraction}] 					Minimum fraction of reads for consensus.")
				println ("--indel_fraction [${params.indel_fraction}] 						Maximum fraction of reads with an indel.")
				println ("--min_cycle_number [${params.min_cycle_number}] 						Minimum cycle number.")
				println ("--max_cycle_number [${params.max_cycle_number}] 						Maximum cycle number.")
				println ("--max_mismatches [${params.max_mismatches}] 						Maximum number of mismatches.")
				println ("--proper_pairs_fraction [${params.proper_pairs_fraction}] 					Minimum fraction of reads that are proper-pairs.")
				println ("--min_consensus_base_quality [${params.min_consensus_base_quality}]				Minimum consensus base quality.")
				println ("--post_trim_read_length [${params.post_trim_read_length}] 					Read length after 5' trimming.")
				println ("--max_normal_vaf [${params.max_normal_vaf}] 					Maximum normal VAF.")
				println ("--min_normal_coverage [${params.min_normal_coverage}] 					Minimum normal coverage.")
				println (" ")
				println ("Indel calling:")
				println (" ")
				println ("--indel_min_bundle_reads [${params.indel_min_bundle_reads}] 					Minimum reads in a bundle for indel calling.")
				println ("--trim_3 [${params.trim_3}] 							Excess bases above this value are trimmed from 3' reads.")
				println ("--trim_5 [${params.trim_5}] 							Bases to trim from 5' reads.")
				println (" ")
				println ("Post processing:")
				println (" ")
				println ("--tri_nuc [${params.tri_nuc}] 							Full path to trinucleotide file, empty string if none.")
				println (" ")


				error ("--help								Print this message.")
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

			if (params.tri_nuc) {
				if (!file(params.tri_nuc).exists()) {
					error ("ERROR: The trinucleotide file does not exist ('${params.tri_nuc}').")
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
