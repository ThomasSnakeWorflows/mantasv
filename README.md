
#### Running the manta pipeline on HG002

##### genologin execution

Setting the environment
```bash
module load bioinfo/samtools-1.9
module load bioinfo/bcftools-1.9
module load bioinfo/snakemake-4.8.0
export MANTA_ROOT="/home/faraut/dynawork/SeqOccin/giab/softwares/manta-1.6.0.centos6_x86_64"
```

Snakemake command

```bash
snakemake --configfile config.yaml \
          --cluster-config cluster.yaml \
          --drmaa " --mem={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" \
          --jobs 4 -p -n
```

##### tatum execution

Setting the environment
```bash
conda activate snakemake
export MANTA_ROOT="/home/tfaraut/Travail/Projets/SeqOccin/giab/softwares/manta/manta-1.6.0.centos6_x86_64"
```
Snakemake command
```bash
snakemake --configfile config.yaml --jobs 4 -p -n
```

##### Comparing with the truthset
```bash
truthset_dir="../../../truthset"
truvari -b $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10_DEL.vcf.gz -c mantasv/manta_DEL.vcf.gz --passonly --includebed $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10.bed -o truvari_del --pctsim 0
```
