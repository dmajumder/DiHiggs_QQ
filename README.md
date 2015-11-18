# Vuoto

Example of shower analysis

Repository with one analysis code to this file linked to fastjet. 
- output histograms, you can adapt to whatever we want

=====================

To make the pythia customized output, modify the configuration file (main12.cc).
118 - 119 make the txt file with a list of momenta (and pid if you want to),
l55 you set the precision, making the file ligther.

l48 you turn hadronization off (easier to analyse, not big impact in physics).
    - if you turn off also parton level you will have the same file in parton level (decaying the higgses only)
    - usually I use the ".decayed" and ".shower" to indicate what is there, a ".decayed" example is uploaded there
    - BB> hhbb files in this format can be found here
    https://cernbox.cern.ch/index.php/s/jGBiMNUrcYDJC68

l50 you say to not decay b-quarks - this will be our "btag"

===================

The file to run is the Vuoto.cc , instructions to run in the reader.
You need to adapt the Makefile to your Fastjet and root paths
It calls Functions.cc and return the number of jets and leptons
There is a separated roothistos.h to declare the histos

In Functions.cc, in special, there is 

- "recojets" function that applies Mass Drop + Filtering by the moment. It fill the histos there.
Prunned mass is out of box (see fastjet manual) n-subjetiness needs FastJetContrib and adaptation in Makefile.

- "isbtagged", to inspect that there is a b-quark inside a jet and check if it is tagable


