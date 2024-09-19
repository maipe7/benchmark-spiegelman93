#!/bin/bash

k0=1e-9
phi0=0.01
etasb=1e19
etaf=1
drho=500.0
g=10.0

rphi1=0.5

delta2=$(python -c "print($k0*$phi0**3*7.0/3.0*$etasb/$etaf)")
delta=$(python -c "print(pow($delta2,0.5))")
w0=$(python -c "print($k0*$drho*$g*$phi0*$phi0/$etaf)")
t0=$(python -c "print($delta/$w0)")
t0yr=$(python -c "print($delta/$w0/(60*60*24*365))")
yres=1000
ndelta=200
outputfolder=output-sp93-ias-ny${yres}-rphi1${rphi1}
echo $outputfolder

echo "# Part added through a script -- overwrites some previous settings:" > 2add.prm

echo "set Output directory                       = ${outputfolder}" >> 2add.prm

echo "set End time                               = $(python -c "print($t0yr*82)")" >> 2add.prm
echo "set Maximum time step                      = $(python -c "print($t0yr*0.1)")" >> 2add.prm
echo "set Maximum first time step                = $(python -c "print($t0yr*0.01)")" >> 2add.prm

echo "subsection Geometry model" >> 2add.prm
echo "  subsection Box" >> 2add.prm
echo "    set X extent = $(python -c "print($delta*$ndelta/$yres)")" >> 2add.prm
echo "    set Y extent = $(python -c "print($delta*$ndelta)")" >> 2add.prm
echo "    set X repetitions = 1" >> 2add.prm
echo "    set Y repetitions = $yres" >> 2add.prm
echo "  end" >> 2add.prm
echo "end" >> 2add.prm

echo "subsection Solver parameters" >> 2add.prm
echo "  subsection Operator splitting parameters" >> 2add.prm
echo "    set Reaction time step                     =  $(python -c "print($t0yr*0.01)")" >> 2add.prm
echo "  end" >> 2add.prm
echo "end" >> 2add.prm

echo "subsection Initial composition model" >> 2add.prm
echo " set Model name = function" >> 2add.prm
echo " subsection Function" >> 2add.prm
echo "   set Function constants = phi0=${phi0}, rphi1=${rphi1}, delta=${delta}, y0=$(python -c "print($delta*50)")" >> 2add.prm
echo "   set Function expression = phi0*(y<y0 ? 1.0 : rphi1 + (1-rphi1)*2.0/( exp( (y-y0)/2.5/delta ) + exp(-(y-y0)/2.5/delta) )  ) ; 0" >> 2add.prm
echo " end" >> 2add.prm
echo "end" >> 2add.prm

echo "subsection Material model" >> 2add.prm
echo "  set Model name = melt global" >> 2add.prm
echo "  subsection Melt global" >> 2add.prm
echo "     set Reference permeability            = ${k0}" >> 2add.prm
echo "     set Reference melt viscosity          = ${etaf}" >> 2add.prm
echo "     set Reference shear viscosity         = ${etasb}" >> 2add.prm
echo "     set Reference bulk viscosity          = ${etasb}" >> 2add.prm
echo "  end" >> 2add.prm
echo "end" >> 2add.prm

dtout=$(python -c "print($t0yr/10.0)")
echo "subsection Postprocess" >> 2add.prm
echo "  subsection Visualization" >> 2add.prm
echo "    set Time between graphical output = ${dtout}" >> 2add.prm
echo "  end" >> 2add.prm
echo " subsection Depth average" >> 2add.prm
echo "    set Number of zones = $(python -c "print($yres*2)")" >> 2add.prm
echo "    set Time between graphical output = ${t0yr}" >> 2add.prm
echo " end" >> 2add.prm
echo "end" >> 2add.prm

cat spiegelman93.prm 2add.prm > spiegelman93_added.prm

# RUN ASPECT SIMULATION
aspect-release spiegelman93_added.prm

awk -f addblanks.awk < $outputfolder/depth_average.txt > $outputfolder/d.txt

gnuplot -persist <<-EOFMarker
    set title "Sp 93 benchmark"
    set xrange [0:${ndelta}]
    set yrange [0:85]
    set nokey
    set out "${outputfolder}/d.png"
    set terminal png
    plot "${outputfolder}/d.txt" using ($ndelta-(\$2)/$delta):((\$3)/$phi0*3+(\$1)/$t0yr-3) every :4 with lines
EOFMarker

gnuplot -persist <<-EOFMarker
    set title "Sp 93 benchmark"
    set xrange [0:${ndelta}]
    set yrange [-2:2]
    set nokey
    set out "${outputfolder}/d80.png"
    set terminal png
    plot "${outputfolder}/d.txt" using ($ndelta-(\$2)/$delta):((\$3)/$phi0) every :80::0::80 with lines
EOFMarker
# rest of script, after gnuplot exits

#done