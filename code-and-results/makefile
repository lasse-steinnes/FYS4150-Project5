#all: compile execute

#compile:
	#g++-10 -o main.out $(wildcard *.cpp) -larmadillo -std=c++11 -fopenmp -Ofast

#execute:
	#./main.out

#plot:
	#python3 plot_sol.py

#plot2:
	#python3 plot_amplification.py


all: compile execute

compile:
	g++ -o main.out $(wildcard *.cpp) -larmadillo -std=c++11 -Ofast -lomp -Xpreprocessor -fopenmp

execute:
	./main.out

plot:
	python3 plot_sol.py

plot2:
	python3 plot_amplification.py
