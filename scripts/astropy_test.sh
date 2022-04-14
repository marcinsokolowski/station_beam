#!/bin/bash -l

#!/bin/bash -l

install_path=~/github/station_beam
export PATH=${install_path}:$PATH
export PYTHON_VERSION=2.7.14

# LOAD MODULES :
# module use /group/mwa/software/modulefiles
# module load MWA_Tools/mwa-sci
# module load astropy
# load modules :
module use /group/mwa/software/modulefiles
# module load MWA_Tools/mwa-sci
# module load pawseytools
# module load python/2.7.14
# module load astropy
module load numpy/1.13.3
module load funcsigs/1.0.2
module load setuptools/38.2.1
module load pytest/3.3.0
module load py/1.5.2
module load astropy/v2.0.12

echo "Modules avail:"
module avail

echo
echo "Modules list:"
module list



python /home/msok/github/station_beam/astropy_test.py > /astro/mwaops/msok/eda2/simulations/EDA2/test.out 2>&1
