# station_beam
SKA-Low station beam models, sensitivity etc

# Web Interface is currently deployed here:

  - sensitivity.skalow.link ( http://3.25.208.103:8000/ )

# local DJANGO webservice can be easily started by doing:
  - cd django/skalow_station_sensitivity/
  - python manage.py runserver
  - then open http://127.0.0.1:8000/ in the browser on local computer or use IP of the computer
    to open from other machine. Port 8000 might have to be opened to enable access.

# deployment procedure :

  - mkdir ~/github/
  - cd ~/github/
  - git clone https://github.com/marcinsokolowski/station_beam.git
  - Execute ~/github/station_beam/install/deploy_beam_models.sh this scripts assumes that the package is installed in ~/github/station_beam/
  - or do it manually : download beam models of individual dipoles from ... Unpack the file:
      - tar zxvf BeamModels.tar.gz
      - move BeamModels to final location
      - modify file ~/github/station_beam/station_beam_config.py so that parameter beam_model_path is set to full path to BeamModels/"






    