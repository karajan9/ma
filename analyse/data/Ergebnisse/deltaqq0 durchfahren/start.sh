#!/bin/bash
#Aufruf:
#./start.sh
#./parsefiles.sh
#./plotdat.gnuplot
#Das Programmverzeichnis von NMRSimulation muss in PATH enthalten sein!
#PATH=$PATH:$HOME/bin/
#vars:
eta=0.0
tauoffset=0
tm=0
lifetime="8192.0"
dwelltime="2.0"
deltaq="1.0"
# deltaqq="0.010983996"
deltaqq="0.015983996"
sigma="1.0"
paramdistribution=1
componentselection=0

model=0
numsamples=4096
numtrajectories=1000000
storeinterval=$numtrajectories
basename="deltaqq${deltaqq}"
trigger=0

interaction=1
#sides=9
#jumprot=0
#jumpangle=10


for i in $(LANG=en seq 0.01 0.01 0.1)
	do
		deltaqq=$i
		basename="deltaqq${deltaqq}"

		echo "deltaqq=$deltaqq"
		../bin/NMRSimulation -j 2 -model $model -interaction $interaction -storeinterval $storeinterval -numtraj $numtrajectories -numsamples $numsamples -fname $basename -eta $eta -dwelltime $dwelltime -lifetime $lifetime -deltaqq $deltaqq -deltaq0 $deltaq -deltaqq0 $deltaqq -sigma $sigma -trigger $trigger -componentselection $componentselection -paramdistribution $paramdistribution &
	
	done
			
# ../bin/NMRSimulation -j 2 -model $model -interaction $interaction -storeinterval $storeinterval -numtraj $numtrajectories -numsamples $numsamples -fname $basename -eta $eta -dwelltime $dwelltime -lifetime $lifetime -deltaqq $deltaqq -deltaq0 $deltaq -deltaqq0 $deltaqq -sigma $sigma -trigger $trigger -componentselection $componentselection -paramdistribution $paramdistribution
	
