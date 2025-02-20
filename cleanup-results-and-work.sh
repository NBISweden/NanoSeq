# Work and results cleanup

shopt -s extglob

rm ./data/genome/genome.fa.*

rm -rf results

cd work

rm -rf  !(apptainer)
