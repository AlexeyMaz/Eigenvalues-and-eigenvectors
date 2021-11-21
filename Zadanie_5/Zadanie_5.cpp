#include <iostream>
#include <math.h>
#include<fstream>
#include <vector>
#include <string>
using namespace std;

void print_matr(int n, vector<vector<double>> A)
{
	cout << "Матрица: " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << A[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

bool check_symmetry(int size, vector<vector<double>> A)
{
	int i, j;
	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++)
			if ((i != j) && (A[i][j] != A[j][i]))
				return 0;

	return 1;
}

vector<double> product(vector<vector<double>> A, vector<double> v, vector<double> res)
{
	for (int i = 0; i < A.size(); i++)
	{
		res[i] = 0;
		for (int j = 0; j < v.size(); j++)
			res[i] += A[i][j] * v[j];
	}
	return res;
}

void to_norm(vector<double>& solution)
{
	double norm = 0;
	for (int i = 0; i < solution.size(); i++)
		norm += solution[i] * solution[i];
	norm = sqrt(norm);

	for (int i = 0; i < solution.size(); i++)
		solution[i] /= norm;
}

int rotation(int size, vector<vector<double>>& A, vector<vector<double>>& solution, double eps)
{
	int result = 1;
	int maxI, maxJ;
	double max, fi;
	vector<vector<double>> matrRotation;
	matrRotation.assign(size, vector<double>(size));
	vector<vector<double>> temp;
	temp.assign(size, vector<double>(size));

	double fault = 0.0;
	for (int i = 0; i < size; i++)
		for (int j = i + 1; j < size; j++)
			fault = fault + A[i][j] * A[i][j];
	fault = sqrt(2 * fault);

	while (fault > eps)
	{
		max = 0.0;
		for (int i = 0; i < size; i++)
			for (int j = i + 1; j < size; j++)
				if (A[i][j] > 0 && A[i][j] > max)
				{
					max = A[i][j];
					maxI = i;
					maxJ = j;
				}
				else if (A[i][j] < 0 && -A[i][j] > max)
				{
					max = -A[i][j];
					maxI = i;
					maxJ = j;
				}

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
				matrRotation[i][j] = 0;
			matrRotation[i][i] = 1;
		}

		if (A[maxI][maxI] == A[maxJ][maxJ])
		{
			matrRotation[maxI][maxI] = matrRotation[maxJ][maxJ] =
				matrRotation[maxJ][maxI] = sqrt(2.0) / 2.0;
			matrRotation[maxI][maxJ] = -sqrt(2.0) / 2.0;
		}
		else
		{
			fi = 0.5 * atan((2.0 * A[maxI][maxJ]) /
				(A[maxI][maxI] - A[maxJ][maxJ]));
			matrRotation[maxI][maxI] = matrRotation[maxJ][maxJ] = cos(fi);
			matrRotation[maxI][maxJ] = -sin(fi);
			matrRotation[maxJ][maxI] = sin(fi);
		}

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				temp[i][j] = 0.0;
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				for (int k = 0; k < size; k++)
					temp[i][j] += matrRotation[k][i] * A[k][j];

		//вычисляем собственные значения
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				A[i][j] = 0.0;
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				for (int k = 0; k < size; k++)
					A[i][j] += temp[i][k] * matrRotation[k][j];

		//снова выделяем максимальный элемент
		fault = 0.0;
		for (int i = 0; i < size; i++)
			for (int j = i + 1; j < size; j++)
				fault = fault + A[i][j] * A[i][j];
		fault = sqrt(2 * fault);

		//вычисляем собственные вектора
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				temp[i][j] = 0.0;
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				for (int k = 0; k < size; k++)
					temp[i][j] = temp[i][j] + solution[i][k] * matrRotation[k][j];

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				solution[i][j] = temp[i][j];
		result++;
	}

	return result;
}

int str_iter(int size, vector<vector<double>> A, vector<double>& solution, double& l, double eps)
{
	double l_prev, x_cur;
	int result = 2;
	x_cur = solution[0];
	solution = product(A, solution, solution);
	l_prev = solution[0] / x_cur; x_cur = solution[0];

	solution = product(A, solution, solution);
	l = solution[0] / x_cur; x_cur = solution[0];

	while (fabs(l - l_prev) > eps)
	{
		solution = product(A, solution, solution);
		l_prev = l;
		l = solution[0] / x_cur; x_cur = solution[0];
		result++;
	}
	
	to_norm(solution);
	return result;
}

void main()
{
	setlocale(LC_ALL, "RUS");
	string tmp;
	vector<vector<double>> matr1;
	vector<double> elem;

	//считывание матрицы из файла
	ifstream matr_f("f.txt");
	if (matr_f.is_open())
	{
		string line;
		while (getline(matr_f, line))
		{
			line += " ";
			tmp = "";
			for (int i = 0; i < line.size(); i++)
			{
				if (line[i] == ' ')
				{
					elem.push_back(stod(tmp));
					tmp = "";
				}
				else tmp += line[i];
			}
			matr1.push_back(elem);
			elem.clear();
		}
	}
	matr_f.close();

	int size = matr1.size();
	print_matr(size, matr1);
	vector<vector<double>> solution1;
	solution1.assign(size, vector<double>(size));
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			solution1[i][j] = 0;
		}
		solution1[i][i] = 1;
	}

	cout << "Метод вращений\n\n";
	cout << "Введите точность: ";
	double eps; cin >> eps;
	if (!check_symmetry(size, matr1))
	{
		cout << endl << "Матрица не симметричная" << endl;
	}
	else
	{
		int iters_count = rotation(size, matr1, solution1, eps);
		cout << endl << endl;
		for (int i = 0; i < size; i++)
		{
			cout << "Собственный вектор h" << i + 1 << ": " << endl;
			for (int j = 0; j < size; j++)
			{
				cout << solution1[j][i] << endl;
			}
			cout << endl;
		}
		cout << "Собственные значения: " << endl;
		for (int i = 0; i < size; i++)
		{
			cout << matr1[i][i] << endl;
		}
		cout << endl << "Число итераций: " << iters_count << endl;
	}

	//МЕТОД ПРЯМОЙ ИТЕРАЦИИ

	vector<vector<double>> matr2;
	//считывание матрицы из файла
	ifstream mar_f("f.txt");
	if (mar_f.is_open())
	{
		string line;
		while (getline(mar_f, line))
		{
			line += " ";
			tmp = "";
			for (int i = 0; i < line.size(); i++)
			{
				if (line[i] == ' ')
				{
					elem.push_back(stod(tmp));
					tmp = "";
				}
				else tmp += line[i];
			}
			matr2.push_back(elem);
			elem.clear();
		}
	}
	mar_f.close();

	vector<double> solution2(size);
	double l = 0;
	for (int i = 0; i < size; i++)
		solution2[i] = 1;

	cout << "Метод вращений\n\n";
	cout << "Введите точность: ";
	eps; cin >> eps;
	if (!check_symmetry(size, matr2))
	{
		cout << endl << "\n\nМатрица не симметричная" << endl;
	}
	else
	{
		int iters_count = str_iter(size, matr2, solution2, l, eps);
		cout << "\n\nСпектральный радиус: " << l << "\n\n";
		cout << "Собственный вектор h:\n";
		for (int i = 0; i < size; i++)
			cout << solution2[i] << endl;
		cout << endl << "Число итераций: " << iters_count << endl;
	}
}