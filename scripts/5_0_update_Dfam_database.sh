#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --partition=parallel
#SBATCH --time=2-0:00:00
#SBATCH --job-name=build_dfam_databse
#SBATCH -o /nfs/scratch/oostinto/stdout/build_dfam_databse.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/build_dfam_databse.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.10.2

SINGULARITY_DIR=/nfs/scratch/oostinto/singularity
LIBRARY_DIR=/nfs/scratch/oostinto/databases/repeatmasker

#downlad singularity image
#cd $SINGULARITY_DIR
#singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
#singularity pull dfam-tetools-latest.sif docker://dfam/tetools:1.92


singularity shell --bind /nfs:/nfs /nfs/scratch/oostinto/singularity/dfam-tetools-latest.sif
cp -r /opt/RepeatMasker/Libraries/ /nfs/scratch/oostinto/databases/repeatmasker/
chown -R $USER /nfs/scratch/oostinto/databases/repeatmasker/Libraries

#copy library from inside container to a local path (should be done from command line)


#download dfam library and check if not corrupted
#cd $LIBRARY_DIR
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.0.h5.gz
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.0.h5.gz.md5
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.4.h5.gz
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.4.h5.gz.md5
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.10.h5.gz
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.10.h5.gz.md5
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.11.h5.gz
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.11.h5.gz.md5
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.12.h5.gz
#wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.12.h5.gz.md5

#verify download
#md5sum -c *.md5

#unpack zipped files
#gzip -d dfam39_full.0.h5.gz
#gzip -d dfam39_full.4.h5.gz
#gzip -d dfam39_full.10.h5.gz
#gzip -d dfam39_full.11.h5.gz
#gzip -d dfam39_full.12.h5.gz

#perform update 
singularity shell --bind /nfs/scratch/oostinto/databases/repeatmasker/Libraries:/opt/RepeatMasker/Libraries /nfs/scratch/oostinto/singularity/dfam-tetools-latest.sif


#download repbase libraries 
#cd $LIBRARY_DIR
#wget "https://github.com/yjx1217/RMRB/raw/1644f60bbd93a435a12fe4f3d3ed6916db4df91b/RepBaseRepeatMaskerEdition-20181026.tar.gz"
wget https://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20181026.tar.gz


#remove 
singularity shell --bind /nfs/scratch/oostinto/databases/repeatmasker/Libraries:/opt/RepeatMasker/Libraries /nfs/scratch/oostinto/singularity/dfam-tetools-latest.sif

#cp -r /opt/RepeatMasker/Libraries/ $LIBRARY_DIR/

