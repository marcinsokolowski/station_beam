-- My cheatshit for sqlite3 : https://www.tutorialspoint.com/sqlite/sqlite_quick_guide.htm

DROP TABLE Sensitivity;

-- List of fields : id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version
CREATE TABLE Sensitivity (
   id             INTEGER PRIMARY KEY AUTOINCREMENT ,      -- star id
   azim_deg       real NOT NULL, -- azimuth from North (0deg) -> East (90deg) -> South (180deg) -> West (270deg)
   za_deg         real NOT NULL, -- zenith angle in degrees
   frequency_mhz  real NOT NULL,
   polarisation   text NOT NULL, -- polarisation X or Y
   lst            real NOT NULL, -- local sidreal time [hours decimal]
   unixtime       real NOT NULL, -- unix time stamp of the time for which sensitivity was generated (should correspond to LST)
   gpstime        real,          -- either unixtime or gps time should be filled
   sensitivity     real NOT NULL, -- sensitivity in units [m^2/K]
   t_sys          real NOT NULL, -- system temperature [Kelvin]
   a_eff          real NOT NULL, -- effective area [m^2]
   t_rcv          real,          -- redudant field should be = t_sys - t_ant
   t_ant          real NOT NULL, -- antenna temperature [Kelvin]
   array_type     int,           -- array type, 1 = EDA1, 2 = EDA2 , 3 = AAVS2 , etc ...
   timestamp      TIMESTAMP WITH TIME ZONE,  -- When this fit was generated
   creator        text,          -- name of creator
   code_version   text           -- What version of the code was used, proposed : eda1 , eda2_beam_fits_trcveda2, eda2_beam_fits_trcveda1 
                                 -- where eda1 - is the code exactly as in the EDA1 paper
);   

CREATE UNIQUE INDEX Sensitivity_unique_index on Sensitivity (array_type,code_version,azim_deg,za_deg,frequency_mhz,unixtime,polarisation);   
   
