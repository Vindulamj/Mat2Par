##You can compile the all and run by just running 
	./run.sh

#Manually compile programs
	##Sequencial Program
	g++ -std=c++11  sequential.cpp -o seq

	##Parallel Program
	g++ -std=c++11  parallel.cpp -o par -fopenmp

	##Optimized Program
	g++ -std=c++11  optimized.cpp -o opt -fopenmp


	#Exicute them
	-First paramete - Number of Iterations

	Ex:- ./seq 10

