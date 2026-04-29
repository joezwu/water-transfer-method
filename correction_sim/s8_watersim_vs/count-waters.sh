jobname=$1

echo "python count-waters.py ${jobname} 'resname VS1' 3.5  > .tmp"
python count-waters.py ${jobname} 'resname VS1' 3.5  > ${jobname}.outw

echo "R CMD BATCH -${jobname} count-waters.R"
R CMD BATCH -${jobname} count-waters.R


