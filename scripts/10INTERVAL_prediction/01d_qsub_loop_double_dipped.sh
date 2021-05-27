dbs=/data/rschubert1/TOPMED_Proteome/*/PCAIR_modeling/06Elastic_net/dbs_out/*rho0.1_zpval0.05.db

for d in $dbs
do
#	echo ${d:74:(-3)}
	qsub -v db=${d} -N ${d:74:(-3)}_dd_prediction 01c_impute_dd_double_dipped_db.sh
done
