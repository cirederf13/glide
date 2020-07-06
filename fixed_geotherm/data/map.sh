#!/bin/sh

coor="cut"
cpt="topo_fred.cpt"
cpt2="age_fred.cpt"
colr="gray"
region="6.108/7.593/45.149/46.409"

psout=setB
datafile=setB.txt

#grdinfo $coor.grd
#GMT grd2cpt $coor.grd -CGMT_$cpt -E10 > $cpt
GMT makecpt -C$colr -T-10000/2000/1 > $cpt
GMT makecpt -Crainbow -I -Dred/purple -T0/30/0.1 > $cpt2
GMT grdgradient $coor.grd -A0 -Nt0.75 -G$coor_i.grd 
GMT grdimage $coor.grd \
    -I$coor_i.grd \
    -E50 \
    -C$colr.cpt \
    -JM14c -R$region \
    -Bf0.5a.5g.5:."":WeSn \
    -P -V -K > $psout.ps
GMT pscoast -J$proj -R -Di -S -O -K >> $psout.ps
sort -g -r -k 4 $datafile | awk '$6==6 {print $1,$2,$4}' | psxy -R -J -O -Sc0.6 -C$cpt2 -P -K -Wblack >> $psout.ps
sort -g -r -k 4 $datafile | awk '$6==2 {print $1,$2,$4}' | psxy -R -J -O -Sc0.5 -C$cpt2 -P -K -Wblack >> $psout.ps
sort -g -r -k 4 $datafile | awk '$6==4 {print $1,$2,$4}' | psxy -R -J -O -Sc0.4 -C$cpt2 -P -K -Wblack >> $psout.ps
sort -g -r -k 4 $datafile | awk '$6==1 {print $1,$2,$4}' | psxy -R -J -O -Sc0.3 -C$cpt2 -P -K -Wblack >> $psout.ps
sort -g -r -k 4 $datafile | awk '$6==3 {print $1,$2,$4}' | psxy -R -J -O -Sc0.2 -C$cpt2 -P -K -Wblack >> $psout.ps
GMT psscale -D15/5/10/0.5 -C$cpt2 -B5 -O >> $psout.ps

open -a Preview $psout.ps
