
### copy these to the shell ###


# conda install mamba -n base -c conda-forge

# mamba env create --name kripi_rnaseq --file kripi_rnaseq.yaml
# mamba env update --name kripi_rnaseq --file kripi_rnaseq.yaml


conda activate kripi_rnaseq


snakemake -np > snakemake_log.txt
snakemake --rulegraph | dot -Tpng > snakemake_rulegraph.png

export TMPDIR="/localscratch"

sbatch --partition=cpu_p \
--qos=cpu \
--mem=8G \
--time=72:00:00  \
--constraint=Lustre_File_System \
--wrap="snakemake -j 30 --cluster-config cluster.json --cluster 'sbatch --constraint=Lustre_File_System -p {cluster.partition} -q {cluster.q} -c {cluster.c} --mem {cluster.mem} -t {cluster.time}'"


################################
