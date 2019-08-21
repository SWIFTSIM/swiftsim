#!/bin/sh

redshift="30.0"
logdensity="-1.00"
./cooling_rates -z $redshift -d $logdensity
python3 plot_result_cooling.py cooling_output_lognH$logdensity.dat $redshift $logdensity
python3 plot_result_heating.py heating_output_lognH$logdensity.dat $redshift $logdensity

