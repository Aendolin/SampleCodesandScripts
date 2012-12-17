!/bin/bash

. /share/apps/gromacs-4.0/bin/GMXRC.bash

sys='ARA'
TA=../Trajectory_Archive

while [ 1 ]; do

COUNTER=0
i=250
while [ $i -le 200000 ]
  do
  b=`echo $i"pico"`
  #echo $b
  if [ -d $TA/$b -a $i > COUNTER ]; then
      COUNTER=$i
  fi
  let i+=250
done

a=`echo $COUNTER"pico"`

#Area per lipid
#--------------
AV=Areas_Volumes/$a

if [ ! -d $AV -a -d $TA/$a ]; then

#how to echo to stderr?
mkdir $AV

echo "0" | g_traj_mpi -f $TA/$a/start.xtc -s $TA/$a/start.tpr -ob $AV/box.xvg -nojump
cat $AV/box.xvg | grep -v "[@#&]" | awk '{print $1,$2*$3/400}' > $AV/box_$a

DT=`date`

echo -e "set term png\n set grid\n set output '$AV/box.png'\n plot '$AV/box_$a' w l title '$sys Area $DT $a'" > area.plt

gnuplot "area.plt"

scp $AV/box.png aendolin@csmlabfs1.cas.usf.edu:/Volumes/Scratch_1/Wiki/POPE/$sys

rm *#

fi

#Order Parameters
#----------------

OP=Order_Parameters/$a

if [ ! -d $OP -a -d $TA/$a ]; then

mkdir $OP

g_order_mpi -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara1.ndx -od $OP/ara1_$a.xvg
cat $OP/ara1_$a.xvg | grep -v "[@#&]" > $OP/ara1_$a

g_order_mpi -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara2.ndx -od $OP/ara2_$a.xvg
cat $OP/ara2_$a.xvg | grep -v "[@#&]" > $OP/ara2_$a

./g_order_sagar -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara_19-22.ndx -od $OP/ara_19-22_$a.xvg -unsat
cat $OP/ara_19-22_$a.xvg | grep -v "[@#&]" > $OP/ara_19-22_$a

./g_order_sagar -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara_22-25.ndx -od $OP/ara_22-25_$a.xvg -unsat
cat $OP/ara_22-25_$a.xvg | grep -v "[@#&]" > $OP/ara_22-25_$a

./g_order_sagar -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara_25-28.ndx -od $OP/ara_25-28_$a.xvg -unsat
cat $OP/ara_25-28_$a.xvg | grep -v "[@#&]" > $OP/ara_25-28_$a

./g_order_sagar -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara_28-31.ndx -od $OP/ara_28-31_$a.xvg -unsat
cat $OP/ara_28-31_$a.xvg | grep -v "[@#&]" > $OP/ara_28-31_$a


./g_order_sagar -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara_40-43.ndx -od $OP/ara_40-43_$a.xvg -unsat
cat $OP/ara_40-43_$a.xvg | grep -v "[@#&]" > $OP/ara_40-43_$a

./g_order_sagar -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara_43-46.ndx -od $OP/ara_43-46_$a.xvg -unsat
cat $OP/ara_43-46_$a.xvg | grep -v "[@#&]" > $OP/ara_43-46_$a

./g_order_sagar -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara_46-49.ndx -od $OP/ara_46-49_$a.xvg -unsat
cat $OP/ara_46-49_$a.xvg | grep -v "[@#&]" > $OP/ara_46-49_$a

./g_order_sagar -f $TA/$a/start.xtc -s $TA/$a/start.tpr -n Order_Parameters/ara_49-52.ndx -od $OP/ara_49-52_$a.xvg -unsat
cat $OP/ara_49-52_$a.xvg | grep -v "[@#&]" > $OP/ara_49-52_$a


cat $OP/ara_19-22_$a | awk '{print $1+3,$2}' > $OP/unsat19-22
cat $OP/ara_22-25_$a | awk '{print $1+6,$2}' > $OP/unsat22-25
cat $OP/ara_25-28_$a | awk '{print $1+9,$2}' > $OP/unsat25-28
cat $OP/ara_28-31_$a | awk '{print $1+12,$2}' > $OP/unsat28-31

cat $OP/ara_40-43_$a | awk '{print $1+3,$2}' > $OP/unsat40-43
cat $OP/ara_43-46_$a | awk '{print $1+6,$2}' > $OP/unsat43-46
cat $OP/ara_46-49_$a | awk '{print $1+9,$2}' > $OP/unsat46-49
cat $OP/ara_49-52_$a | awk '{print $1+12,$2}' > $OP/unsat49-52

cat $OP/ara1_$a | awk '{if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 6 || $1 == 9 || $1 == 12 || $1 == 15 || $1 == 16 || $1 == 17 || $i == 18) print $1,$2}' > $OP/sat1
cat $OP/sat1 $OP/unsat19-22 $OP/unsat22-25 $OP/unsat25-28 $OP/unsat28-31 | sort -g > $OP/ara1_comb_$a

cat $OP/ara2_$a | awk '{if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 6 || $1 == 9 || $1 == 12 || $1 == 15 || $1 == 16 || $1 == 17 || $i == 18) print $1,$2}' > $OP/sat2
cat $OP/sat2 $OP/unsat40-43 $OP/unsat43-46 $OP/unsat46-49 $OP/unsat49-52 | sort -g > $OP/ara2_comb_$a

DT=`date`

echo -e "set term png\n set grid\n set yrange [-0.1:0.4]\n set output \"$OP/ara1.png\"\n plot '$OP/ara1_comb_$a' w l title '$sys ara1 Order $DT $a'" > order.plt

gnuplot "order.plt"

echo -e "set term png\n set grid\n set yrange [-0.1:0.4]\n set output \"$OP/ara2.png\"\n plot '$OP/ara2_comb_$a' w l title '$sys ara2 Order $DT $a'" > order.plt

gnuplot "order.plt"

scp $OP/ara1.png aendolin@csmlabfs1.cas.usf.edu:/Volumes/Scratch_1/Wiki/POPE/$sys
scp $OP/ara2.png aendolin@csmlabfs1.cas.usf.edu:/Volumes/Scratch_1/Wiki/POPE/$sys

rm *#

fi


#Electron Density
#----------------

ED=Electron_Density

if [ ! -d $ED/$a -a -d $TA/$a ]; then

mkdir $ED/$a

echo "2" | g_density_mpi -f $TA/$a/start.xtc -s $TA/$a/start.tpr -ei $ED/hoh.dat -o $ED/$a/hoh_density_$a.xvg -dens electron -sl 200  
cat $ED/$a/hoh_density_$a.xvg | grep -v "[@#&]" > $ED/$a/hoh_density_$a

echo "1" | g_density_mpi -f $TA/$a/start.xtc -s $TA/$a/start.tpr -ei $ED/ara.dat -o $ED/$a/ara_density_$a.xvg -dens electron -sl 200 
cat $ED/$a/ara_density_$a.xvg | grep -v "[@#&]" > $ED/$a/ara_density_$a

echo -e "set term png\n set grid\n set output \"$ED/$a/density.png\"\n plot '$ED/$a/hoh_density_$a' w l, '$ED/$a/ara_density_$a' w l title '$sys Electron Density $DT $a'" > ed.plt

gnuplot "ed.plt"

scp $ED/$a/density.png aendolin@csmlabfs1.cas.usf.edu:/Volumes/Scratch_1/Wiki/POPE/$sys

rm *#

fi

sleep 4h

done
