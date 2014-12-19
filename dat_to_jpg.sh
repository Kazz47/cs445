# GNUPLOT script to generate images
echo 'set term jpeg' >> temp_plot.gps

for i in *.dat
do
    filename=`echo $i | sed 's/dat/jpg/'`

    echo "set output '$filename'" >> temp_plot.gps
    echo "plot for [col=2:3] '$i' using 1:col with lines" >> temp_plot.gps
    echo "reset" >> temp_plot.gps
done

gnuplot temp_plot.gps
rm temp_plot.gps

