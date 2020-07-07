#!/bin/sh

rm map.ps

coor="cut"
cpt="topo_fred.cpt"
cpt2="age_fred.cpt"
colr="gray"
psout="map.ps"
region="6.108/7.593/45.149/46.409"
datafile="setE.txt"

grdinfo $coor.grd
#GMT grd2cpt $coor.grd -CGMT_$cpt -E10 > $cpt
GMT makecpt -C$colr -T-10000/2000/1 -V > $cpt
GMT makecpt -Crainbow -I -Dred/purple -T0/14/0.1 > $cpt2
GMT grdgradient $coor.grd -A0 -Nt0.75 -G$coor_i.grd -V
GMT grdimage $coor.grd \
    -I$coor_i.grd \
    -E50 \
    -C$colr.cpt \
    -JM14c -R$region \
    -Bf0.2a.2g.2:."":WeSn \
    -P -V -K > $psout
#GMT psbasemap -R -J -B2S:"X":/2W:"Y":eWnS -O -K >> $psout
GMT pscoast -J$proj -R -Di -S -O -K >> $psout
psxy ../data/fault.txt -R -J -O -P -K -Wthickest -Sf0.25  >> $psout
psxy true_coords_faults.txt -R -J -O -P -K -Wthick -Sf0.5  >> $psout
psxy true_coords_main_units.txt -R -J -O -P -K -Wthick -Sf0.5  >> $psout
sort -g -r -k 4 $datafile | awk '$6==6 {print $1,$2,$4}' | psxy -R -J -O -Sc0.6 -C$cpt2 -P -K -Wblack >> $psout
sort -g -r -k 4 $datafile | awk '$6==2 {print $1,$2,$4}' | psxy -R -J -O -Sc0.5 -C$cpt2 -P -K -Wblack >> $psout
sort -g -r -k 4 $datafile | awk '$6==4 {print $1,$2,$4}' | psxy -R -J -O -Sc0.4 -C$cpt2 -P -K -Wblack >> $psout
sort -g -r -k 4 $datafile | awk '$6==1 {print $1,$2,$4}' | psxy -R -J -O -Sc0.3 -C$cpt2 -P -K -Wblack >> $psout
sort -g -r -k 4 $datafile | awk '$6==3 {print $1,$2,$4}' | psxy -R -J -O -Sc0.2 -C$cpt2 -P -K -Wblack >> $psout
GMT psscale -D15/5/10/0.5 -C$cpt2 -B2 -O >> $psout
