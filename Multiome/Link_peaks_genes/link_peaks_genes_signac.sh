JOB_DIR="correlation_jobs"

DONORS=$(cat donors.txt)

for donor in ${DONORS[*]} ; do

    output_job=${JOB_DIR}"/correlations_"${donor}"_SBATCH.sh"

    echo -e "#!/usr/bin/env bash
#SBATCH --job-name=correlations_${donor}
#SBATCH --output=correlations_${donor}.out
#SBATCH --error=correlations_${donor}.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

        set -e
        ~/start_Rscript.sh ../correlation_one_donor.R ${donor}
        " > ${output_job}
done