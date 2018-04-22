#!/bin/bash
# Aufruf:
# ./start.sh
# Das Programmverzeichnis von NMRSimulation muss in PATH enthalten sein!
# z.B. in .bashrc: PATH=$PATH:$HOME/bin/

# VARS:

# j: Type of simulation ID. Options are:
# 1 -- 1d spectra and omega10/omega30 correlation maps (obsolete, will be
#      changed in future version)
# 2 -- FID
# 3 -- 2D Spectrum
# 4 -- Stimulated Echo Signal
# 5 -- Stimulated Echo Signal 3D (obsolete, will be removed in future version!)
# 6 -- Hahnecho
# j=0

# tm: mixing time
tm=0
# tau: evolution time
tau="0.1171875e-6"
# lifetime: mean lifetime of a state of the motional model. Calculated with
# -log(rand) * lifetime
# lifetime="750.0"
lifetime="1.0"
# dwelltime: time between two samples
dwelltime="0.5e-6"
# trigger: trigger time
trigger=0

# numsamples: number of samples (2^n preferred)
numsamples=4096
# numtrajectories: number of trajectories (to calculate an average)
# numtrajectories=10000000
# numtrajectories=100000 == 50s bei lt = "4.41927079e-05"
# numtrajectories=10000 == 6min bei lt = "3.45255530e-07"
# numtrajectories=4000000
numtrajectories=400000
# storeinterval=100000
storeinterval=20000

# model: motional model ID. Options are:
# 0 -- Isotrop
# 1 -- n side jump (euler angles)
# 2 -- random two side jump
# 3 -- hexagonal ice jump model
# 4 -- cone jump model
# 5 -- Anderson model
# 6 -- n side jump double (euler angles, with second orientation for chemical
#      shift)
model=0

# interaction: interaction IDs, which interaction should be used. Options are:
# 0 -- quadupolar first order
# 1 -- quadrupolar second order
# 2 -- gaussian distributed random frequencies
# 3 -- quadrupolar second order MAS
# 4 -- quadrupolar second order and chemical shift
interaction=1
# eta: asymmetry parameter of the quadrupole
eta=0.0
# deltaq: first order quadrupole coupling constant (anisotropy parameter).
# Also used for:
# deltaqq0: factor for random interactions, mostly equal to -deltaq but
# without $ e q $, which is determined randomly.
deltaq="0.0"
# deltaqq: second order quadrupole coupling constant. Also used for:
# deltaqq0: factor for random interactions, mostly equal to -deltaqq but
# without $ e q $, which is determined randomly.
deltaqq="1000000.0"
# sigma: width of the random $ C_Q $ for the quadrupolar second order
# interaction
sigma="1.0"
# paramdistribution: How eta and deltaq/deltaqq are distributed in the isotrope
# model
paramdistribution=1
# componentselection: choose which component of the frequency is used when the
# frequency is a composition. 0 means all components. For quadrupole second
# order 1 is $\omega_{10} $ and 2 $ \omega_{30} $.
componentselection=0

# basename: the base file name for all output files. If a loop is used to start
# multiple simulations in succession, this name should be changed in the loop
# so the outputs go to different files.
basename="lineshape_lifetime!${lifetime}!"


# sides=9
# jumprot=0
# jumpangle=10

# A loop from generated values
# for i in $(LANG=en seq 0.01 0.005 0.05)

# A loop from a list of values
#  "0.132176" "13.2845" "7863.38" "4.34644e14"

declare -a arr=("3.86980325e-08" "8.01593435e-09")
# declare -a arr=("3.86980325e-08")
for i in "${arr[@]}"
	do
		lifetime=$i
		basename="lineshape_lifetime!${lifetime}!"

		echo "lifetime: ${lifetime}"
		echo "deltaqq: ${deltaqq}"
		echo "numtrajectories: ${numtrajectories}"
		srun ../bin/NMRSimulation -j 2 -model $model -interaction $interaction -storeinterval $storeinterval -numtraj $numtrajectories -numsamples $numsamples -fname $basename -eta $eta -dwelltime $dwelltime -lifetime $lifetime -deltaqq $deltaqq -deltaq0 $deltaq -deltaqq0 $deltaqq -sigma $sigma -trigger $trigger -componentselection $componentselection -paramdistribution $paramdistribution -tau $tau -tm $tm &

		sleep 2s
	done
			
# ../bin/NMRSimulation -j 2 -model $model -interaction $interaction -storeinterval $storeinterval -numtraj $numtrajectories -numsamples $numsamples -fname $basename -eta $eta -dwelltime $dwelltime -lifetime $lifetime -deltaqq $deltaqq -deltaq0 $deltaq -deltaqq0 $deltaqq -sigma $sigma -trigger $trigger -componentselection $componentselection -paramdistribution $paramdistribution
	
