#!/bin/bash

# Download the tables of the publicly available equations of state used in SWIFT
wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/EoS/planetary_HM80_HHe.txt
wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/EoS/planetary_HM80_ice.txt
wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/EoS/planetary_HM80_rock.txt

mv planetary_HM80_HHe.txt HM80_HHe.txt
mv planetary_HM80_ice.txt HM80_ice.txt
mv planetary_HM80_rock.txt HM80_rock.txt
