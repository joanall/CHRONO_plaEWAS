#! /bin/bash
#SBATCH --job-name=mt
#SBATCH --partition=short
#SBATCH --mail-type=begin       
#SBATCH --mail-type=end         
#SBATCH --mail-user=joana.llaurado@isglobal.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30gb
#SBATCH --output=/PROJECTES/BISC_OMICS/analyses/BiSC_23/PROJECTS_CL_JLLP/EWAS_MT_JLLP/EWAS_MT_%A_%a.txt

module purge > /dev/null 2>&1
module load lang/R

#R CMD BATCH 3_Build_Run_Toy_Model.R $model

# $i = number model = from 1 to 26 

i=$(($SLURM_ARRAY_TASK_ID))

Rscript script/1_1_EWAS_DMRs_mealtime.R $i
#R CMD script/ewas_mealtime/1_EWAS_mealtime_v1.R $i

# To run it : sbatch --array=1-6 your_script.sh