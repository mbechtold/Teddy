Teddy Tool

By Florian Zabel
Department of Environmental Sciences
University of Basel, Switzerland

contact: florian.zabel@unibas.ch

Model Description Paper:
Zabel, Poschlod (2023): The Teddy tool v1.1: temporal disaggregation of daily climate model data for climate impact analysis. Geoscientific Model Development.
DOI: https://doi.org/10.5194/gmd-16-5383-2023

## Large-scale Application of Teddy Tool on a High-Performance Computing Cluster  
**Case Study: VSC-HPC, Belgium**

---

### **Step 1: Compilation of Code**  
Compile the code using a MATLAB license (Tier-2, WICE Cluster).  

**Command to load MATLAB module and compile:**
```
module load MATLAB/2023b
mcc -m TeddyTool_v1_2p.m -a WFDE5Factors.m WFDE5Factors_LST.m wfde5_tiles.m leap_year.m \
local_solar_time.m calc_coordinates_global_land.m merge_tiles.m write_netcdf_tiles.m \
write_netcdf_dailyval.m parsave.m \
-I /apps/leuven/rocky8/icelake/2022b/software/MATLAB/2023b/toolbox/finance \
-I /apps/leuven/rocky8/icelake/2022b/software/MATLAB/2023b/toolbox/fininst \
-I /apps/leuven/rocky8/icelake/2022b/software/MATLAB/2023b/toolbox/parallel
```

### **Step 2: Transfer and Execution**  

Transfer the executable TeddyTool_v1_2p to the Tier-1 cluster.
Create the SLURM job file and configuration file (TeddyTool.ini).

```
#!/bin/bash
#SBATCH -t 72:00:00                 # Job runtime
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Tasks per node
#SBATCH -A 2022_203                 # Project account
#SBATCH --cpus-per-task=124         # CPUs per task
#SBATCH --mem-per-cpu=2012          # Memory per CPU (MB)
#SBATCH -o s85_2015_log.txt         # Standard output log
#SBATCH -e s85_2015_out.txt         # Error log

ulimit -s unlimited                  # Unlimited stack size
cd /dodrio/scratch/projects/2022_200/project_output/rsda/vsc31786/teddy/GFDL-ESM4/s85_2015
module load MCR/R2023b.9             # Load MATLAB Compiler Runtime module
./s85_TeddyTool_v1_2p.sh /readonly/dodrio/apps/RHEL8/zen3-ib/software/MCR/R2023b.9/R2023b
```
