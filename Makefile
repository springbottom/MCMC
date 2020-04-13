
fdg: fdg.cc fdg.h MCMC.o
	g++ -g -Wfatal-errors -stdlib=libc++ -o fdg fdg.cc MCMC.o

MCMC.o: MCMC.cc MCMC.h
	g++ -g -Wfatal-errors -stdlib=libc++ -c MCMC.cc
