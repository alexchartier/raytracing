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


Ionosonde sites used to generate the classifier:


Site    Lat  Lon

AL945   45.1   276.4
AT135   38.0   23.5
BC840   40.0   254.7
CAJ25   -22.7   315.0
RL055   51.5   359.4
CO765   64.9   212.0
DB049   50.1   4.6
EA655   37.1   353.3
EI764   64.7   212.9
TR169   69.6   19.2
MO155   55.5   37.3
FF055   51.7   358.5
FZA05   -3.9   321.6
GA762   62.4   215.0
GR135   -33.3   26.5
IC435   37.1   127.5
IF843   43.8   247.3
JI915   -12.0   283.2
JR055   54.6   13.4
MH453   52.0   122.5
MHJ45   52.0   122.5
NDA81   81.4   342.5
PF765   65.1   212.6
PSJ55   -51.6   302.1
PQ055   50.0   14.6
PRJ18   18.5   292.9
RO045   41.9   12.5
EB045   40.8   0.5
SMJ65   67.0   309.1
SMK25   67.0   309.1
THJ75   77.5   290.8
THJ75   77.5   290.8
WP937   37.9   284.5

    

