#!/bin/bash 

export chimes_revision=6fa8ee748086c9b94988bd5f814087f879ddf6ca

mkdir chimes-data
cd chimes-data 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/chimes_main_data.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/cross_sections_B87.hdf5 

mkdir HM12_cross_sections
cd HM12_cross_sections
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/redshifts.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z0.000_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z0.200_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z0.400_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z0.600_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z0.800_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z1.000_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z1.200_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z1.400_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z1.600_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z1.800_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z2.000_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z3.000_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z4.000_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z5.000_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z6.000_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z7.000_cross_sections.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/HM12_cross_sections/z8.000_cross_sections.hdf5 

cd .. 
mkdir EqAbundancesTables 
cd EqAbundancesTables 
mkdir colibre_HHe 
cd colibre_HHe 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z0.000_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z0.200_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z0.400_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z0.600_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z0.800_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z1.000_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z1.200_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z1.400_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z1.600_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z1.800_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z2.000_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z3.000_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z4.000_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z5.000_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z6.000_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z7.000_eqm.hdf5 
wget https://bitbucket.org/richings/chimes-data/raw/$chimes_revision/EqAbundancesTables/colibre_HHe/z8.000_eqm.hdf5 
