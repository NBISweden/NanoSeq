# Work and results cleanup

shopt -s extglob

rm ./data/genome/genome.fa.*

rm -rf results

if [ -d "work" ]; then
	cd work
	rm -rf  !(apptainer)
else
  echo "Work directory does not exist."
fi
