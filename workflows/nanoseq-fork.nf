/* 
----------------------------------------------------------------------------------------
Main workflow definition
----------------------------------------------------------------------------------------
*/

// Feature flags

	nextflow.preview.topic = true

// Import modules



// Process samplesheet

	def ch_samplesheet = Channel
		.fromPath("${params.sampleSheet}")
		.splitCsv(header: true)

// Workflow definition

workflow NANOSEQ_FORK {

	// TEST

	println ("Main workflow")

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
