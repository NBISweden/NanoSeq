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

	println ("Verify subworkflow")

}
