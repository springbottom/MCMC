


fdg: fdg.cc fdg.h MCMC.o
	g++ -g -Wfatal-errors -o fdg fdg.cc MCMC.o

MCMC.o: MCMC.cc MCMC.h
	g++ -g -Wfatal-errors -c MCMC.cc
