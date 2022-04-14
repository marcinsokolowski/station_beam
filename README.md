# station_beam
SKA-Low station beam models, sensitivity etc

# Web Interface is currently deployed at these locations :
  - https://sensitivity.skalow.link/ or https://sensitivity.skalow.link/backup 
    

# Deployment procedure for the python package :

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

    It can also be done manually : download beam models of individual dipoles from https://www.dropbox.com/sh/imb2002vwhmoyxo/AAB36xD5UYGE5PblAiEkWYLma?dl=0 to BeamModels directory, unpack the files:
      - tar zxvf EDA.tar.gz
      - tar zxvf SKALA2.tar.gz
      - tar zxvf SKALA4.tar.gz    
      - move BeamModels to final location, for example : mv BeamModels ~/aavs-calibration/
      - modify file ~/github/station_beam/station_beam_config.py so that parameter beam_model_path is set to full path to BeamModels/"


  - deployment of databases use script install/deploy_databases.sh which takes 3 parameters:
    - parameter 1 is the location of database file on the target system (default ~/github/station_beam/sql/)
    - parameter 2 is the URL to the EDA2 SQLite database (currently on my Dropbox)
    - parameter 3 is the URL to the AAVS2 SQLite database (currently on my Dropbox)


    
# Example commands :

  - Example commands using sensitivity databases for AAVS2 and EDA2 stations are provided in the header of sensitivity_db.py script 

  - Example command to generated sensitivity at a specified time (GPS time), frequency (160 MHz), for AAVS2 station at the position of Hydra-A
    radio-galaxy, using Average Embedded Element (AEE) pattern of SKALA4 antenna in the AAVS2 station:

     - cd python/
     - python ./eda_sensitivity.py --freq=160 -p None -g 1320891467  -m analytic --ra=333.607249950 --dec=-17.02661111 --outsens_file=HydA_aavs2_sensitivity --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err  --nos11 --header=HEADER  --use_beam_fits --station_name=SKALA4 --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_aavs2.txt --projection=aee




    