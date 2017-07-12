echo Compile Programms
g++ -std=c++11  sequential.cpp -o seq
g++ -std=c++11  parallel.cpp -o par -fopenmp
g++ -std=c++11  optimized.cpp -o opt -fopenmp

echo Run Programs
echo Sequencial Program
./seq 20

echo Parallel Program
./par 20

echo Optimized Parallel Program
./opt 20
