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

}
