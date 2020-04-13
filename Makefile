


fdg: fdg.cc fdg.h MCMC.o
	g++ -Wfatal-errors -o fdg fdg.cc MCMC.o

MCMC.o: MCMC.cc MCMC.h
	g++ -Wfatal-errors -c MCMC.cc
