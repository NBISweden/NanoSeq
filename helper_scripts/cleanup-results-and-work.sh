# Work and results cleanup

shopt -s extglob

### Edit this line to your cloned repository
cd ~/work/nanoseq

# Clean indexes
rm data/genome/genome.fa.*
rm data/test/genome.fa.*

# Clean results
rm -rf results

# Clean work except pulled images
if [ -d "work" ]; then
	cd work
	rm -rf  !(apptainer)
else
  echo "Work directory does not exist."
fi
