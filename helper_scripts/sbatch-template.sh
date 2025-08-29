#!/bin/bash -l

# EDIT THESE VALUES

REPO_PATH=~/nanoseq
TMUX_SESSION_NAME=nanoseq
PIXI_COMMAND="pixi run -e apptainer nextflow main.nf -profile apptainer,baseHPC --samplesheet data/test/samplesheet.csv --fasta data/test/genome.fa --account HOURS_BILLING_PROJECT"
PIPELINE_NAME=NanoSeq


### DO NOT EDIT BELOW THIS LINE ###

if [ -d "$REPO_PATH" ]; then

	cd $REPO_PATH

	TMUX_VERSION_OUTPUT=$(tmux -V 2>&1)

	if [ $? -ne 0 ]; then
		echo "Error: tmux not found"
	else
		if echo "$TMUX_VERSION_OUTPUT" | grep -q "^tmux "; then
			MAJORVERSION=$(echo "$TMUX_VERSION_OUTPUT" | sed 's/^tmux //;s/\..*$//')
			if [ "$MAJORVERSION" -eq 1 ]; then
				tmux new-session -s "$TMUX_SESSION_NAME" -d
				tmux send-keys -t "$TMUX_SESSION_NAME" "$PIXI_COMMAND" C-m
				echo "Found tmux version $MAJORVERSION"
				echo "$PIPELINE_NAME pipeline started in tmux session: '$TMUX_SESSION_NAME', with command:"
				echo "$PIXI_COMMAND"
				echo "Monitor workflow progress in '.nextflow.log' or with 'squeue -u $USER'"
			elif [ "$MAJORVERSION" -ge 2 ]; then
				tmux new -s "$TMUX_SESSION_NAME" -d /bin/bash -c "$PIXI_COMMAND"
				echo "Found tmux version $MAJORVERSION"
				echo "$PIPELINE_NAME pipeline started in tmux session: '$TMUX_SESSION_NAME', with command:"
				echo "$PIXI_COMMAND"
				echo "Monitor workflow progress in '.nextflow.log' or with 'squeue -u $USER'"
			fi
		else
			echo "Error: Unable to parse tmux version"
		fi
	fi

else

	echo "Error: '$REPO_PATH' does not exist, set to the $PIPELINE_NAME repository path."

fi
