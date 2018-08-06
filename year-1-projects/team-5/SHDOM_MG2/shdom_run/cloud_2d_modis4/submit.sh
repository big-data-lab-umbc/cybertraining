# Author: Meng Gao
# DATE: 2018

#wave=0.86 #2.1
#reff=10
deg=60 #60 #0,10,20,30,40,50,60,70
#deg=30 #0,10,20,30,40,50,60,70
pi=`echo "4*a(1)" | bc -l`
rad=`echo "$deg*($pi/180)" | bc -l`
mu0=`echo "c($rad)" | bc -l`
#flag3d=3 #0 #0: 3d, 3: 1d
flag3d=0 #0 #0: 3d, 3: 1d
nx=50
ny=50
#nz=16
#nz=32
nz=50
#nz=100 #50 #32
memory=5000

for wave in 0.86 2.1
do
    for index in 22
    do
	echo "deg wave reff"
	echo $deg $wave $reff
	name=output_${deg}_${wave}_index${index}_nz${nz}_${flag3d}	
	prp=atmos${wave}_index${index}.prp
	
	sed 's/<NX>/'$nx'/g' template_run_lookup > tmp1
	sed 's/<NY>/'$ny'/g' tmp1 > tmp2
	sed 's/<NZ>/'$nz'/g' tmp2 > tmp3
	sed 's/<MEMORY>/'$memory'/g' tmp3 > tmp1
	sed 's/<MU0>/'$mu0'/g' tmp1 > tmp2
	sed 's/<WAVE>/'$wave'/g' tmp2 > tmp1
	sed 's/<3d>/'$flag3d'/g' tmp1 > tmp2
	sed 's/<PRP>/'$prp'/g' tmp2 > tmp1
	sed 's/<NAME>/'$name'/g' tmp1 > run_${deg}_${wave}_index${index}_${flag3d}
	
	bash run_${deg}_${wave}_index${index}_${flag3d}
    done
done