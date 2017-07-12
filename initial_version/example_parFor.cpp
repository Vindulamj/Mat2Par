#include <cmath>
#include <iostream>

using namespace std;

int main()
{
	const int size = 256;
	double sinTable[size];

	#pragma omp parallel for
	for(int n=0; n<size; ++n)
	sinTable[n] = std::sin(2 * M_PI * n / size);

// the table is now initialized
	for(int i; i<size; i++)
		cout << sinTable[i] << endl;

	return 0;
}
