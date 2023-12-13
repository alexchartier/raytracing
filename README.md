# simulations


You will need to install PHaRLaP from https://www.dst.defence.gov.au/opportunity/pharlap-provision-high-frequency-raytracing-laboratory-propagation-studies
Install nc_utils from  https://github.com/alexchartier/nc_utils
Install ionogram_classifier from https://github.com/alexchartier/ionogram_classifier

WSPR validation
    download the WSPR data from https://www.wsprnet.org/drupal/wsprnet/spots, then
    run_sami_raytracing.py 

GPS position validation
    Install the Madrigal globalDownload.py, then 
    gps_pos_val.py  [this will fail due to lack of files the first time, but give instructions on how to get them]
    gps_utils.py will generate a sitelist from the Madrigal HDF5 file
    Then run gps_pos_val.py when you have the necessary HDF5 and Madrigal files
    

Ionosonde validation
    ionosonde_val.py

Note you will need to install various Python dependencies, as well as the PHaRLAP and m_map packages referenced in startup.m 


    

