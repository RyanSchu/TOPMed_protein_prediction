dbs=/home/rschubert1/scratch/TOPMed_Proteome/07Valid_CV/valid_cv_filtered_dbs/*db

for d in $dbs
do
	qsub -v db=${d} -N ${d:74:(-3)}_prediction 01impute_db.sh
done
