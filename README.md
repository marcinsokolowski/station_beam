# station_beam
SKA-Low station beam models, sensitivity etc

# Web Interface is currently deployed at these locations :

  - sensitivity.skalow.link ( http://3.25.208.103:8000/ )

# local DJANGO webservice can be easily started by doing:
  - cd django/skalow_station_sensitivity/
  - python manage.py runserver
  - then open http://127.0.0.1:8000/ in the browser on local computer or use IP of the computer
    to open from other machine. Port 8000 might have to be opened to enable access.

# deployment procedure :

  - requirements : python3 from anaconda package (https://www.anaconda.com/products/individual-d) usually includes all the requirement packages (astropy, h5py, numpy etc), but when 
    using other python distributions (including standard package included in Linux distributions) the following steps may be required:
    - pip install astropy
    - pip install h5py
    - pip install numpy

  - mkdir ~/github/
  - cd ~/github/
  - git clone https://github.com/marcinsokolowski/station_beam.git
  - Execute ~/github/station_beam/install/deploy_beam_models.sh this scripts assumes that the package is installed in ~/github/station_beam/ , parameters are :
     - first parameter is final location where BeamModel files are copied, default ~/aavs-calibration/BeamModels/
     - parameters 2, 3, 4 are URLs to EDA.tar.gz, SKALA2.tar.gz, SKALA4.tar.gz currently on my Dropbox

    It can also be done manually : download beam models of individual dipoles from ... Unpack the file:
      - tar zxvf BeamModels.tar.gz
      - move BeamModels to final location
      - modify file ~/github/station_beam/station_beam_config.py so that parameter beam_model_path is set to full path to BeamModels/"




  - deployment of databases use script install/deploy_databases.sh which takes 3 parameters:
    - parameter 1 is the location of database file on the target system (default ~/github/station_beam/sql/)
    - parameter 2 is the URL to the EDA2 SQLite database (currently on my Dropbox)
    - parameter 3 is the URL to the AAVS2 SQLite database (currently on my Dropbox)


    





    