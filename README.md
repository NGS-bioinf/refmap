# refmap
Multi purpose reference mapper

When in doubt, write alen.suljic@mf.uni-lj.si

The repository includes the scripts and containers for both short reads (_sr) and long reads (_lr). The example is shown for (_sr), but works exactly the same for (_lr).

Prerequisites:
- Singularity (https://singularity-tutorial.github.io/01-installation/)

Usage:
1. Build the Singularity container:
```
sudo singularity build refmap_sr.sif 20240902_refmap_sr.def
#if no sudo privileges:
singularity build --fakeroot refmap_sr.sif 20240902_refmap_sr.def
```
2. Enter refmap_sr.sif:
```
singularity shell /path/to/refmap_sr.sif
```
     
3. Run the refmap_sr script:
```
bash refmap_sr1.sh reference /path/to/data.fastq.gz
```
4. The procedure is the same for refmap_lr.sh

Best use example:
1. To avoid issues with data naming for (_sr), follow the convention in example (sample name should include only "-".
Read orientation is delimited by _.):
   - sample-1_R1.fastq.gz
   - sample-1_R2.fastq.gz
   - ...
  
2. Create a directory where desired.
3. Copy the refmap_sr.sh script and "reference" directory to created directory.
4. Enter refmap_sr.sif container.
5. Optional: If required, open the hicov.sh with your favourite text editor and adjust pipeline parameters (thread count, quality settings, etc...).
6. Run the refmap_sr.sh script. The script takes reference name as positional argument 1 and input data (fastq.gz) directory as positional argument 2.
7. Additional information is provided in header section in refmap_sr.sh script.
