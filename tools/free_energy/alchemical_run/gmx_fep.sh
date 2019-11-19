#!/bin/bash

# _________ read inputs from the galaxy wrapper and define some variables ____________ 

lam=$1
iter=$((lam+1))

mkdir MDP
mkdir data
mkdir traj
 
FREE_ENERGY=`pwd`
MDP=$FREE_ENERGY/MDP

set -e

for i in `seq 0 $lam`
 do
   cp em_steep.mdp em_steep_$i.mdp
   sed -i "s/%L%/$i/" em_steep_$i.mdp
   cp nvt.mdp nvt_$i.mdp
   sed -i "s/%L%/$i/" nvt_$i.mdp 
   cp npt.mdp npt_$i.mdp
   sed -i "s/%L%/$i/" npt_$i.mdp 
   cp md.mdp md_$i.mdp
   sed -i "s/%L%/$i/" md_$i.mdp  
 done
mv *.mdp $MDP

for (( i=0; i<$iter; i++ ))
do
    LAMBDA=$i

    # A new directory will be created for each value of lambda 

    mkdir Lambda_$LAMBDA
    cd Lambda_$LAMBDA

# _______ ENERGY MINIMIZATION STEEP _______  

    echo "Starting minimization for lambda = $LAMBDA..." 

    mkdir EM
    cd EM

    # Iterative calls to grompp and mdrun to run the simulations

    gmx grompp -f $MDP/em_steep_$LAMBDA.mdp -c $FREE_ENERGY/morph.gro -p $FREE_ENERGY/morph.top -o min$LAMBDA.tpr

    gmx mdrun -deffnm min$LAMBDA

    sleep 10


# _______ NVT EQUILIBRATION _______ 

    echo "Starting constant volume equilibration..."

    cd ../
    mkdir NVT
    cd NVT

    gmx grompp -f $MDP/nvt_$LAMBDA.mdp -c ../EM/min$LAMBDA.gro -p $FREE_ENERGY/morph.top -o nvt$LAMBDA.tpr

    gmx mdrun -deffnm nvt$LAMBDA

    echo "Constant volume equilibration complete."

    sleep 10

# _______ NPT EQUILIBRATION _______ 

    echo "Starting constant pressure equilibration..."

    cd ../
    mkdir NPT
    cd NPT

    gmx grompp -f $MDP/npt_$LAMBDA.mdp -c ../NVT/nvt$LAMBDA.gro -p $FREE_ENERGY/morph.top -t ../NVT/nvt$LAMBDA.cpt -o npt$LAMBDA.tpr

    gmx mdrun -deffnm npt$LAMBDA

    echo "Constant pressure equilibration complete."

    sleep 10

# ________ PRODUCTION MD ___________ 

    echo "Starting production MD simulation..."

    cd ../
    mkdir Production_MD
    cd Production_MD

    gmx grompp -f $MDP/md_$LAMBDA.mdp -c ../NPT/npt$LAMBDA.gro -p $FREE_ENERGY/morph.top -t ../NPT/npt$LAMBDA.cpt -o md$LAMBDA.tpr

    gmx mdrun -deffnm md$LAMBDA

    echo "Production MD complete."

    # End
    echo "Ending. Job completed for lambda = $LAMBDA"

    cd $FREE_ENERGY
done

cp Lambda_*/Production_MD/*.xvg data/
tar cf data.tar data/

cp Lambda_*/Production_MD/*.trr traj/
tar cf traj.tar traj/

exit;


