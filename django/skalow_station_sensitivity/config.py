from string import Template
import os

# configuration :
station_beam_path_tmp = "$HOME/github/station_beam"
t = Template( station_beam_path_tmp )
station_beam_path = t.substitute(os.environ)


sensitivity_db_path = station_beam_path + "/sql/" 
save_output_path = "/tmp/station_beam/"
