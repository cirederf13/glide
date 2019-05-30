#!/bin/sh

rm plots/*

#script to make figures and movie of output maps

#Setting some initial things
region=-R6.108/7.593/45.149/46.409
res="0.005"
proj="M8c"
path="../data"
makecpt -Cpolar -T0/0.01/0.0001 -Mgray > apt.cpt
makecpt -Cdrywet -Dwhite/gray -T0./1.8/0.01 > apt2.cpt
makecpt -Cseis -I -Dgray/gray -T0./1./0.001 > apt3.cpt
makecpt -Cdrywet -Dgray/gray -T0.2/.35/0.001 > apt4.cpt
makecpt -Cdrywet -Dgray/gray -T0./1./0.001 > apt4.cpt

gmtset HEADER_FONT_SIZE=28p
gmtset PAPER_MEDIA=A3
gmtset TICK_LENGTH = 0.c
gmtset FRAME_WIDTH = 0.1c

# prepare the topgraphy
blockmedian $path/topo.txt $region -I$res -Vl > $path/cut_block.xyz
surface $path/cut_block.xyz -G$path/cut.grd $region -I$res -V -T0.27
grdgradient $path/cut.grd -A0 -Nt0.75 -G$path/ill.grd -V

#The loop is for every file in maps
for i in 0020 0060

do

if [ ${i: -4} == ".txt" ]
then
  echo ""$i" file ends in .txt"
  continue
fi

ps2=plots/"$i.ps"

#Interpolate surfaces
awk '{print $1,$2,$3}' $i | nearneighbor -G$i.grd $region -I$res -N4/0 -S15k -E-9999 -V
awk '{print $1,$2,$3}' unc/$i | nearneighbor -Gunc/x.grd $region -I$res -N4/0 -S0.1 -E1 -V
awk '{print $1,$2,$4}' unc/$i | nearneighbor -Gunc/$i.grd $region -I$res -N4/0 -S0.1 -E1 -V

#Mask out regions, colour them grey, where the uncertainty is greater than....
grdclip unc/x.grd -Gclip.grd -Sa.99/NaN -V
gmtset COLOR_FOREGROUND=114/114/114

grdmath $i.grd clip.grd OR = new.grd
grdmath unc/$i.grd clip.grd OR = unc/new.grd
grdmath unc/x.grd clip.grd OR = unc/new_redV.grd
#grdmath unc/x.grd clip.grd OR 0.1225 MUL SQRT = unc/new_redV.grd

t1=`awk '{print $7}' $i.txt`
t2=`awk '{print $9}' $i.txt`

#create basemap,grdimage, sample locations and scale
psbasemap $region -J$proj -Y3 -B0.5S:"X":/0.5W:"Y":WeSn -K  > $ps2
grdimage new.grd -J$proj $region -Capt2.cpt -I$path/ill.grd -K -O >> $ps2

#pscoast -J$proj -R -Di -I1 -S -O -K >> $ps2 
awk '$4<4 && $2>'$t2' {print $1,$2}' ../data/ages.txt | psxy -Gwhite -J$proj -R -Wblack -Sc0.4 -V -O -K >> $ps2
awk '$4<'$t1'&& $4>'$t2' {print $1,$2}' ../data/ages.txt | psxy -Gblack -J$proj -R -Wblack -Sc0.4 -V -O -K >> $ps2
awk '{print $1,$2,$3}' ../data/solution.txt | psxy -Capt2.cpt -J$proj -R -Wblack -Ss0.8 -V -O -K >> $ps2
psxy ../mask -L -R -J -O -P -K -W5 -V >> $ps2
psxy ../data/fault.txt -R -J -O -P -K -Wthickest -Sf0.25  >> $ps2
psscale -D4./-1/6/0.5h -B0.4:"mm/yr": -I -Capt2.cpt -O -K >> $ps2

pstext $i.txt  $region -J$proj -X2.5c -Y1c -Wblack -Gwhite -N -V -O -K >> $ps2

#Same as above but for the resolution
grdimage unc/new.grd -J$proj -R -B0.5S:"X":/0.5W:"Y":weSn -Capt3.cpt -X6c -Y-1c -I$path/ill.grd -O -K >> $ps2
awk '$4<'$t1'&& $4>'$t2' {print $1,$2}' ../data/ages.txt | psxy -Gblack -J$proj -R -Wblack -Sc0.4 -V -O -K >> $ps2
psxy ../data/fault.txt -R -J -O -P -K -Wthickest -Sf0.25  >> $ps2
psscale -D4./-1/6/0.5h -B0.2:"resolution": -I -Capt3.cpt -O -K >> $ps2
pstext $i.txt $region -J$proj -X2.5c -Y1c -Wblack -Gwhite -N -O -V -K >> $ps2

#Same as above for reduced variance
grdimage unc/new_redV.grd -J$proj -R -B1S:"X":/1W:"Y":weSn -Capt4.cpt -X6c -Y-1c -I$path/ill.grd -O -K >> $ps2
awk '$4<'$t1'&& $4>'$t2' {print $1,$2}' ../data/ages.txt | psxy -Gblack -J$proj -R -Wblack -Sc0.4 -V -O -K >> $ps2
psxy ../data/fault.txt -R -J -O -P -K -Wthickest -Sf0.25  >> $ps2
psscale -D4./-1/6/0.5h -B0.25:"reduced variance": -I -Capt4.cpt -O -K >> $ps2
pstext $i.txt $region -J$proj -X2.5c -Y1c -Wblack -Gwhite -N -O -V >> $ps2

rm $i.grd
rm new.grd
rm unc/*.grd
rm clip.grd
done

#cleaning up the folders
rm apt*

#Converting the ps files to jpg for movie 
ps2raster -Tj plots/*ps -P -V -A

