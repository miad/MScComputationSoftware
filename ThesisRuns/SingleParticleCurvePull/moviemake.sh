 #!/bin/bash
tempdir="MovieTemp"
rm -rf $tempdir
mkdir $tempdir
max=5000
workers=200
framesAfter=1000
framesBefore=1000
movieName="movie.mkv"

seq -f "%05g" 0 ${max} | parallel -j$workers "gnuplot -e \"set terminal pngcairo size 1600,1000 enhanced font 'Veranda,20'; set xra [2.5:2.85]; set yra [-5E-3:0]; set xlabel 'Re[k] /(µm^{-1})'; set ylabel 'Im[k] /(µm^{-1})';  set key right bottom; plot 'KCurves/0_{}.dat' u 1:2 ps 5 lt 1 pt 7 title 'Eigenstates, kinetic operator', 'KFounds/0_{}.dat' u 1:2 ps 5 lt 3 pt 10 lw 4 title 'Eigenstates, full Hamiltonian'\" > $tempdir/{}.png"

seq $framesBefore `expr $framesAfter + $framesBefore` | parallel -j$workers "A=\`expr {} + $max\`;cp $tempdir/\`printf \"%05g.png\" $max\` $tempdir/\`printf \"%05g.png\" \$A \`"
seq $max -1 0 | parallel -j1 "A=\`expr {} + $framesBefore \`; mv $tempdir/\`printf \"%05g.png\" {}\` $tempdir/\`printf \"%05g.png\" \$A \`"

B=`expr $framesBefore + 1`
seq 1 $framesBefore | parallel -j$workers "cp $tempdir/\`printf \"%05g.png\" $B \` $tempdir/\`printf \"%05g.png\" {}\`"



rm -f $movieName
#ffmpeg -r 50 -f image2 -i $tempdir/%05d.png -f mp4 -q:v 0 -vcodec mpeg4 -r 50 $movieName
ffmpeg -r 75 -i $tempdir/%05d.png -c:v libx264 -preset veryslow -qp 0 -r 75 $movieName