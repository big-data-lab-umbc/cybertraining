# Author: Meng Gao
# DATE: 2018

tempfile=template_run_mie

for wave in 0.86 2.1
do
    for((reff=8;reff<=17;reff=reff+1))
    do
	echo $reff
	sed 's/<WAVE>/'$wave'/g' $tempfile > tmp1
	sed 's/<REFF>/'$reff'/g' tmp1 > tmp2
	sed 's/<NAME>/scat_'$wave'_'$reff'/g' tmp2 > run_mie_${wave}_${reff}
	csh run_mie_${wave}_${reff}
    done
done

