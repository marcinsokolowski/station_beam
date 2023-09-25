#!/bin/bash

#!/bin/bash

export PATH=~/github/station_beam/scripts/common/:$PATH

# simulate_total_power_vs_time_all.sh 300 0 9.676768 1 "-q -l -b" - list 
simulate_total_power_vs_time.sh total_power_160_160_MHz.txt - 300 0 9.676768 - - - "--station_name=SKALA4 --use_beam_fits --trcv_type=trcv_aavs2 --size=512"

