
#all: compile execute

#compile:
	#g++-10 -Wall -o main.out $(wildcard *.cpp) -larmadillo -fopenmp -std=c++11 -Ofast

#execute:
	#./main.out

all: compile execute

compile:
	g++-10 -o main.out $(wildcard *.cpp) -larmadillo -std=c++11 -fopenmp -Ofast

execute:
	./main.out

plot:
	python3 plot_sol.py
