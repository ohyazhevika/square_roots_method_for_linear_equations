#pragma once
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
using namespace std;

double Rand(int range)
{
	double ans = (double)(rand() % (2 * abs(range) + 1) - abs(range));
	ans += ans*1.0/rand();
	return ans;
}
class SLE
{
	bool GetY()
	{
		for (int i = 0; i < N; i++)
		{
			double S = f[i];
			for (int k = k0(i); k < i; k++)
				S -= Q[k][i - k] * y[k];
			if (Q[i][0] == 0)
				return false;
			y[i] = S / Q[i][0];
		}
		for (int i = N - 1; i >= 0; i--)
		{
			double S = y[i];
			for (int k = i + 1; k <= kN(i); k++)
				S -= Q[i][k - i] * X[k];
			if (Q[i][0] == 0)
				return false;
			X[i] = S / Q[i][0];
		}
		return true;
	}

	bool GetT()
	{
		for (int i = 0; i < N; i++)
		{
			double S = P[i][L - 1];
			for (int k = k0(i); k < i; k++)
				S -= pow(Q[k][i - k], 2);
			if (S < 0)
			{
				return false;
			}
			S = sqrt(S);
			Q[i][0] = S;

			for (int j = i + 1; j <= kN(i); j++)
			{
				S = P[j][i - j + L - 1];
				for (int k = k0(j); k < i; k++)
					S -= Q[k][i - k] * Q[k][j - k];
				Q[i][j - i] = S / Q[i][0];
			}
		}
		return true;
	}
	int k0(int i)//первый элемент ленты в строке i
	{
		if (i < L)
			return 0;
		else
			return i - L + 1;
	}
	int kN(int i)
	{
		if (i >= N - L)
			return N - 1;
		else
			return i + L - 1;
	}
	bool FullTapeMatriceInput(fstream& file)
	{
		if (!file.is_open())
		{
			cerr << "File not found.\n";
			return false;
		}
		double elem;
		for (int i = 0; i < N; i++)
		{
			int j = 0;
			for (j; j < k0(i); j++)
				file >> elem;
			for (j; j <= i; j++)
				file >> P[i][j - i + L - 1];
			for (j; j < N; j++)
				file >> elem;
		}
		return true;
	}
	bool fInput(fstream& file)
	{
		for (int i = 0; i < N; i++)
			file >> f[i];
		return true;
	}
	bool ModifiedTapeMatriceInput(fstream& file)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < L; j++)
				file >> P[i][j];
		return true;
	}
public:
	double
		** P,//нижняя половина ленты исходной симметричной ленточной матрицы
		** Q, // лента матрицы T искомая
		** Q_precise; //лента матрицы Т заданная
	double
		* f, //правая часть системы
		* y,// решение системы T^(T)y=f
		* y_precise, // точное решение системы T^(T)y=f
		* X, //вектор-столбец Х, получаемый в ходе решения системы Ax=f, где А - исходная симметричная ленточная матрица
		* X_precise; //точное решение системы
	int
		N, //размерность
		L; //ширина ленты
	SLE(int N, int L) : N(N), L(L)
	{
		Q_precise = new double* [N];
		Q = new double* [N];
		P = new double* [N];
		f = new double[N];
		y = new double[N];
		X = new double[N];
		X_precise = new double[N];
		y_precise = new double[N];
		for (int i = 0; i < N; i++)
		{
			X[i] = 0;
			X_precise[i] = 0;
			Q[i] = new double[L];
			P[i] = new double[L];
			Q_precise[i] = new double[L];
			for (int j = 0; j < L; j++)
			{
				Q[i][j] = Q_precise[i][j] = P[i][j] = 0;
			}
		}
	}
	~SLE()
	{
		for (int i = 0; i < N; i++)
		{
			delete[] P[i];
			delete[] Q[i];
			delete[] Q_precise[i];
		}
		delete[] P;
		delete[] f;
		delete[] X;
		delete[] X_precise;
		delete[] y;
		delete[] y_precise;
	}

	bool FullTapeSystemInput(fstream& file)
	{
		if (!file.is_open())
		{
			cerr << "File not found.\n";
			return false;
		}
		return FullTapeMatriceInput(file)&&fInput(file);
	}
	bool ModifiedTapeSystemInput(fstream& file)
	{
		if (!file.is_open())
		{
			return false;
		}
		return ModifiedTapeMatriceInput(file) && fInput(file);
	}
	bool MatricesQandFInput(fstream& file)
	{
		if (!file.is_open())
		{
			return false;
		}

	}
	void GenerateQ()
	{
		for (int i = 0; i < N; i++)
		{
			X_precise[i] = Rand(10);
			for (int j = 0; j < min(L, N - i); j++)
				Q_precise[i][j] = Rand(10); //точная верхнедиагональная матрица сгенерирована
		}
		for (int i = 0; i < N; i++)
			Q_precise[i][0] = fabs(Q_precise[i][0]);
	}
	void GenerateSLE1()//генерация матрицы A (то есть Р) непосредственно
	{//НЕ ГАРАНТИРУЕТСЯ, что будет получена положительно определённая матрица!!!
		for (int i = 0; i < N; i++)
		{
			X_precise[i] = Rand(10);
			for (int j = max(0, L-i-1) ; j < L; j++)
				P[i][j] = Rand(10);
		}
		//получаем вектор f=Ax
		for (int i = 0; i < N; i++)
		{
			double S = 0;
			for (int j = k0(i); j <= i; j++)
				S += P[i][j - i + L - 1] * X_precise[j];
			for (int j = i + 1; j <= kN(i); j++)
				S += P[j][i - j + L - 1] * X_precise[j];
			f[i] = S;
		}
	}
	void GenerateSLE(int kk = 0)
	{
		GenerateQ();
		if (kk != 0)
			for (int i = 0; i < N; i++)
				Q_precise[i][0] /= pow(10, kk);//сгенерировали плохо обусловленную матрицу

		for (int i=0; i<N; i++)
			for (int j = k0(i); j <= i; j++)
			{
				double S = 0;
				for (int k = k0(i); k <= j; k++)
					S += Q_precise[k][i-k] * Q_precise[k][j-k];
				P[i][j-i+L-1] = S; //A=T^t*T - получили коэффициенты системы
			}
		// находим значение вектора y_precise=Tx
		for (int i = 0; i < N; i++)
		{
			double S = 0;
			for (int j = i; j <= kN(i); j++)
				S += Q_precise[i][j - i] * X_precise[j];
			y_precise[i] = S;
		}
		// находим точное значение вектора f=(T^t)y
		for (int i = 0; i < N; i++)
		{
			double S = 0;
			for (int j = k0(i); j <= i; j++)
				S += Q_precise[j][i - j] * y_precise[j];
			f[i] = S;
		}

	}


	bool Solve()
	{
		if (!GetT())
		{
		//	cout << "The matrice is not positively defined, thus the system cannot be solved this way.\n";
			return false;
		}
		if (!GetY())
		{
		//	cout << "The system cannot be solved.\n";
			return false;
		}
		return true;
	}
	void PreciseSolutionInput(fstream& file)
	{
		for (int i = 0; i < N; i++)
			file >> X_precise[i];
	}
	bool GetSolution() 
	{
		for (int i = 0; i < N; i++)
			cout << "X[" << i << "] = " << X[i] << '\n';
		return true;
	}
	bool GetPreciseSolution()
	{
		for (int i = 0; i < N; i++)
			cout << X_precise[i] << '\n';
		return true;
	}
	bool ShowSolutions()
	{
		for (int i = 0; i < N; i++)
			cout << "X_precise[" << i << "] = " << X_precise[i] << "   " << "X[" << i << "] = " << X[i] << '\n';
		return true;
	}
	//////////////////////
	double RelativeError()
	{
		double ans = (X_precise[0] == 0) ? fabs(X_precise[0] - X[0]) : fabs((X_precise[0] - X[0]) / X_precise[0]);
		for (int i = 1; i < N; i++)
		{
			double error = (X_precise[i] == 0) ? fabs(X_precise[i] - X[i]) : fabs((X_precise[i] - X[i]) / X_precise[i]);
			if (error > ans)
				ans = error;
		}
		return ans;
	}

	void GetMatriceP()
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < L; j++)
				cout << P[i][j] << "  ";
			cout << '\n';
		}
	}

	void GetF()
	{
		for (int i = 0; i < N; i++)
			cout << f[i] << '\n';
	}
	void GetMatriceQ()
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < L; j++)
				cout << Q[i][j] << "  ";
			cout << '\n';
		}
	}
	void GetMatriceQ_precise()
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < L; j++)
				cout << Q_precise[i][j] << "  ";
			cout << '\n';
		}
	}
};
/*
0 0 0 16
0 0 4 17
0 -8 -6 14
0 8 1 6
*/