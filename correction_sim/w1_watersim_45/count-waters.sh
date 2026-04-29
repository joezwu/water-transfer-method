jobname=$1

echo "python count-waters.py ${jobname} 4.5 'resname VS1' 'resname VS2'  > .tmp"
python count-waters.py ${jobname} 4.5 'resname VS1' 'resname VS2'   > ${jobname}.outw

echo "R CMD BATCH -${jobname} count-waters.R"
R CMD BATCH -${jobname} count-waters.R


