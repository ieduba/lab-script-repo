#! /bin/bash
#SBATCH --partition=cao
#SBATCH --nodelist=node176
#SBATCH --nodes=1 #Choose the node number 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --job-name="211207_demultiplex_Hera_211206_VL00112_70_AAANVKTM5" # annotate job name
#SBATCH -o /rugpfs/fs0/cao_lab/store/zlu/Report_slurm/211207_demultiplex_Hera_211206_VL00112_70_AAANVKTM5.out # STDOUT
#SBATCH --mail-user=zlu@rockefeller.edu # when job is finished, send an email to your address  
#SBATCH --mail-type=ALL 

#accept a run folder, a sample sheet and a output folder; then do demultiplex for the run folder and input the raw data to the ouput foler

run_folder="/rugpfs/fs0/cao_lab/scratch/jcao/bcl2_files/211206_VL00112_70_AAANVKTM5"
sample_sheet="/ru-auth/local/home/zlu/Sample_sheet/2021/211207_Hera.csv"
output_folder="/rugpfs/fs0/cao_lab/scratch/zlu/211207_Hera_Data/"

echo "---------------start demultiplex-----------------------"
echo $(date)
echo "run folder is $run_folder"
echo "sample sheet is $sample_sheet"
echo "output_folder is $output_folder"

# Empty the output folder 
# rm -rf $output_folder

#create the output folder:
echo 
echo start making the output_folder
mkdir -p $output_folder/report

# do the demultiplex
/rugpfs/fs0/cao_lab/scratch/jcao/software/anaconda3/bin/bcl2fastq --runfolder-dir $run_folder -o $output_folder --sample-sheet $sample_sheet --reports-dir $output_folder/report --barcode-mismatches 1 --create-fastq-for-index-reads --no-lane-splitting --use-bases-mask Y*,I8n*,I8n*,Y* --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0

remove un-demultiplexed reads
echo remove undetermined reads
rm $output_folder/Undetermined*.fastq.gz

echo "------------------demultiplex done -----------------"
