#!/bin/bash
cp control.lead12 control.in

cp geometry.lead1 geometry.in
<your aims binary here>
mv lead_self_energy lead_1

cp geometry.lead2 geometry.in
<your aims binary here>
mv lead_self_energy lead_2
