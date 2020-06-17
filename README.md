# cutAndQC

Snakemake QC pipeline for Cut&Tag or Cut&Run projects. 

# SETUP

Clone this repository into your project directory:

```
# clone to your local called 'my-project'
git clone git@github.com:maxsonBraunLab/cutAndQC.git my-project

# create directory for your fastq files
cd my-project
mkdir -p data/raw

# link your fastqs to here
ln -s /path/to/fastq/files/* data/raw

```

The configuration happens in config.yml

 - Specify the bowtie2 index for the genome you are aligning to.
 * Specify path to bowtie2 indices for genomes in config.yml

# Execution

Setup snakemake profile to run on compute cluster:

SLURM: follow [these instructions](https://github.com/Snakemake-Profiles/slurm)

Example of how to run pipeline

```
snakemake --use-conda --profile slurm -j 60 --latency-wait 60
```

To change runtime parameters for indvidual rules, you can provide a cluster configuration file via the `--cluster-config` flag, for example:

```
snakemake --use-conda --profile slurm --cluster-config src/cluster.yml -j 60 --latency-wait 60
```

Here is an example cluster.yml that sets defaults, and rule specific resources:

```
cat src/cluster.yml
__default__:
    memory: "20G"
    threads: 2
    partition: "exacloud"
    time: 14400 #time in seconds
bowtie2:
    threads: 8
    time: 86400
```


The final QC document will be generated under `data/multiqc/multiqc_report.html`, where alignment statistics, duplication rates, and fastqc results can be compared across all samples. Estimated library complexity is plotted for each sample, as well as the fingerpint plot generated using [deepTools plotFingerprint](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html). 
