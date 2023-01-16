#include <iostream>
#include <math.h>
#include <fstream>
#include "System_of_linear_equations.h"
using namespace std;

int main()
{
	srand(time(NULL));
	//test 3
	//int N = 10, L = 3;
	//double aver_relative_error = 0;
	//int bad_matrices = 0;
	//int number_of_tests = 1;
	//int i = 1;
	//for (i; i <= number_of_tests; )
	//{
	//	
	//	SLE sys(N, L);
	//	sys.GenerateSLE();
	//	double rel_err = 0;
	//	if (sys.Solve())
	//	{
	//		
	//		//if ((rel_err = sys.RelativeError()) && rel_err < 0.001)
	//		{
	//			i++;
	//			aver_relative_error += sys.RelativeError();
	//			cout << "Precise matrice Q:\n";
	//			sys.GetMatriceQ_precise();
	//			cout << "\nMatrice P: \n";
	//			sys.GetMatriceP();
	//			cout << "\nF: \n";
	//			sys.GetF();
	//			cout << "\nX_precise:\n";
	//			sys.GetPreciseSolution();
	//			cout << "\n";
	//		}
	//	}
	//}
	//aver_relative_error /= i;
	//cout << "Relative error is " << aver_relative_error << endl;
	//cout << "Number of bad matrices: " << bad_matrices << "\n";

	///////////////////////////////////////////////////////////////

	//גגמה טח פאיכא
	fstream file("Generated_system1.txt");
	if (!file.is_open())
	{
		cerr << "File not found.\n";
		return false;
	}
	SLE system(10, 3);
	system.ModifiedTapeSystemInput(file);
	system.PreciseSolutionInput(file);
	if (system.Solve())
	{
		system.GetMatriceQ();
		cout << endl;
		system.GetSolution();
		cout << "Relative error is " << system.RelativeError();
	}
	else cout << "Cannot solve";

	return 0;
}

