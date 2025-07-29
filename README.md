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
module load MATLAB/2023b-r1
mcc -m TeddyTool_v1_2p.m -a WFDE5Factors.m WFDE5Factors_LST.m wfde5_tiles.m leap_year.m \
local_solar_time.m calc_coordinates_global_land.m write_netcdf_tiles.m \
write_netcdf_dailyval.m parsave.m
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
#SBATCH --cpus-per-task=91         # CPUs per task
#SBATCH --mem-per-cpu=2612          # Memory per CPU (MB)
#SBATCH -o s85_2015_log.txt         # Standard output log
#SBATCH -e s85_2015_out.txt         # Error log

ulimit -s unlimited                  # Unlimited stack size
cd /dodrio/scratch/projects/2022_200/project_output/rsda/vsc31786/teddy/GFDL-ESM4/s85_2015
module load MCR/R2023b.9             # Load MATLAB Compiler Runtime module
./s85_TeddyTool_v1_2p.sh /readonly/dodrio/apps/RHEL8/zen3-ib/software/MCR/R2023b.9/R2023b
```
To stay within the 72 hours walltime, it was necessary to divide the processing into 5 year blocks and a northern and southern hemisphere job that run sequentially and write into the same output netcdf file.
The first line of the TeddyTool.ini file was:
84,9   !latitude (comma seperated, boundaries when using tiles ...
and
9,-56  !latitude (comma seperated, boundaries when using tiles ...

The whole TeddyTool.ini was:

```
9,-56                          !latitude (comma seperated, boundaries when using tiles ->see line 18)
-180,180                        !longitude (comma seperated, boundaries when using tiles ->see line 18)
1                               !select hourly time step for temporal disaggregation (1=1-hourly,2=2-hourly,etc.)
../../OUTPUT/obsclim            !output directory
../../WFDE5factors              !directory containing preprocessed subdaily factors derived from wfde5 data
../../WFDE5_v2p1                        !directory containing hourly WFDE5 reference data
../../ISIMIP3b/InputData/climate/atmosphere/bias-adjusted/global/daily  !directory containing ISIMIP climate model data
GFDL-ESM4 ! GSWP3-W5E5, GFDL-ESM4,IPSL-CM6A-LR,MPI-ESM1-2-HR,MRI-ESM2-0,UKESM1-0-LL      !select ISIMIP climate model(s): GFDL-ESM4,IPSL-CM6A-LR,MPI-ESM1-2-HR,MRI-ESM2-0,UKESM1-0-LL
ssp585 ! ssp585,ssp370,ssp126,historical,picontrol     !select ISIMIP scenario(s): ssp585,ssp370,ssp126,historical,picontrol
2015-2020                               !select timeperiod for ISIMIP climate data: startyear-endyear
tas,hurs,rsds,rlds,ps,sfcwind,pr        !select variable(s): tas,hurs,rsds,rlds,ps,sfcwind,pr
1980-2019  !select time period for hourly wfde5 data that is used to determine the diurnal profile
11   !select time window of +- n days around DOY to search for similar conditions in the past (default=11)
0    !use LST (local solar time) (1=default) or UTC (0) for data processing?
1    !consider precipitation on consecutive days: 1=yes, 0=no; doy window flag (ln 13) must be >= 1
0    !write NaN values for precipitation (0=no=default) in case that no precipitation event can be found in the historical hourly reference dataset (1=yes: in this case mass and energy are not preserved!)
91  !number of parallel workers
5    !value>0:use tiles and set tile size; 0:use coordinates in line 1,2
```
