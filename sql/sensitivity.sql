CREATE TABLE STARS (
   id             serial,      -- star id
   azim_deg       real NOT NULL, -- azimuth from North (0deg) -> East (90deg) -> South (180deg) -> West (270deg)
   za_deg         real NOT NULL, -- zenith angle in degrees
   frequency_mhz  real NOT NULL,
   lst            real NOT NULL, -- local sidreal time [hours decimal]
   unixtime       real NOT NULL, -- unix time stamp of the time for which sensitivity was generated (should correspond to LST)
   gpstime        real,          -- either unixtime or gps time should be filled
   sensitvity     real NOT NULL, -- sensitivity in units [m^2/K]
   t_sys          real NOT NULL, -- system temperature [Kelvin]
   a_eff          real NOT NULL, -- effective area [m^2]
   t_ant          real NOT NULL, -- antenna temperature [Kelvin]
   array_type     int,           -- array type, 1 = EDA1, 2 = EDA2 , 3 = AAVS2 , etc ...
   timestamp      TIMESTAMP WITH TIME ZONE,  -- When this fit was generated
   creator        text,          -- name of creator
   code_version   text           -- What version of the code was used, proposed : eda1 , eda2_beam_fits_trcveda2, eda2_beam_fits_trcveda1 
                                 -- where eda1 - is the code exactly as in the EDA1 paper
);   
   
   
--   magnitude      real,      -- mean magnitude from the measurements
--   sigma_mag      real,      -- rms error of the measured magnitude
--   name           text,      -- star name
--   min_mag        real,      -- minimum magnitude among all the measurements
--   max_mag        real,      -- maximum magnitude among all the measurements
--   no_measurements integer,  -- number of measurements for a star
--   mag_cat        real,       -- do tego pola zapisywana bedzie jasnosc 
                             -- gwiazdy z katalogu odniesienia
