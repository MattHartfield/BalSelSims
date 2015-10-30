BalSelSims:

Simulation of balancing selection simulation with sampling (finite population effects)
Selection; reproduction; mutation based on deterministic recursions, then sampling
Repeats for 'Length' generations

Run by executing:
./BalSelSims N s rec sex self gc reps
Where:
- N is the population size
- s is the fitness disadvantage of homozygotes
- rec is recombination rate
- sex is rate of sex (a value between 0 = obligate asex, and 1 = obligate sex)
- self is selfing rate
- gc is gene conversion
- reps is how many times to introduce linked neutral allele

Note that haplotypes are defined as:
x1 = ab
x2 = Ab
x3 = aB
x4 = AB

Genotypes defined as:
g11 = g1 = ab/ab
g12 = g2 = Ab/ab
g13 = g3 = aB/ab
g14 = g4 = AB/ab
g22 = g5 = Ab/Ab
g23 = g6 = Ab/aB
g24 = g7 = Ab/AB
g33 = g8 = aB/aB
g34 = g9 = aB/AB
g44 = g10 = AB/AB

Simulation uses routines found with the GNU Scientific Library (GSL)
(http://www.gnu.org/software/gsl/)
Since GSL is distributed under the GNU General Public License 
(http://www.gnu.org/copyleft/gpl.html), you must download it 
separately from this file.