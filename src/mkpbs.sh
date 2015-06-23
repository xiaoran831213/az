#!/bin/bash

dst="test.qs"

echo "#!/bin/bash -login" > "$dst"
echo "#PBS -l nodes=1:ppn=1" >> "$dst"
echo "#PBS -l walltime=03:00:00"
echo "#PBS -l mem=4096M"

echo "cd $PBS_O_WORKDIR"

