#include <stdio.h>
#include <mpi.h>
#include <math.h> 
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int** initMatrix(int n)
{
	int **a = new int*[n];
	for (int i = 0; i < n; i++)
		a[i] = new int[n];

	return a;
}

int** fillInMatrix(int** a, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i][j] = i;// rand() % 100;
		}
	}
	return a;
}

int** fillInZeroMatrix(int** a, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i][j] = 0;
		}
	}
	return a;
}

int main(int argc, char* argv[])
{
	int rank, size;
	double t1, t2;
	int matrixSize;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0)
	{
		cout << "Enter matrix size\n";
		cin >> matrixSize;
	}
	MPI_Bcast(&matrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int **A, **B, **C;
	if (rank == 0)
	{
		A = initMatrix(matrixSize);
		B = initMatrix(matrixSize);
		C = initMatrix(matrixSize);
		A = fillInMatrix(A, matrixSize);
		B = fillInMatrix(B, matrixSize);
		C = fillInZeroMatrix(C, matrixSize);
		if (matrixSize <= 20)
		{
			cout << "\nMatrix A:\n";
			for (int i = 0; i < matrixSize; i++)
			{
				for (int j = 0; j < matrixSize; j++)
				{
					cout << A[i][j] << " ";
				}
				cout << "\n";
			}
			cout << "\nMatrix B:\n";
			for (int i = 0; i < matrixSize; i++)
			{
				for (int j = 0; j < matrixSize; j++)
				{
					cout << B[i][j] << " ";
				}
				cout << "\n";
			}
			cout << "\n";
		}
	}

	int residue = matrixSize % size;
	int rowAmount = matrixSize / size;
	int columnAmount = rowAmount;

	//Если остаток > 0 и это один из последних потоков,
	//то он должен принять на 1 больше столбцов
	if (residue > 0 && size - rank <= residue)
	{
		columnAmount++;
	}
	int *column = new int[matrixSize*columnAmount];

	//0 поток рассылает столбцы матрицы остальным потокам
	if (rank == 0)
	{
		t1 = MPI_Wtime();
		int r = 0;
		for (int m = 1; m < size; m++)
		{
			int amount = columnAmount;
			//Если есть остаток, то последним потокам посылаем на 1 больше столбцов
			if (residue > 0 && size - m <= residue)
			{
				amount++;
			}
			int *sendingColumn = new int[matrixSize*amount];

			for (int l = 0; l < amount; l++)
			{
				for (int j = 0; j < matrixSize; j++)
				{
					sendingColumn[j + l*matrixSize] = B[j][m*columnAmount + l + r];
				}
			}

			/*for (int e = 0; e < matrixSize*amount; e++)
			{
				cout << "rank: " << m << " elem: " << e << " - " << sendingColumn[e] << "\n";
			}
			cout << "\n";*/
			MPI_Send(sendingColumn, matrixSize*amount, MPI_INT, m, m, MPI_COMM_WORLD);
			//если m потоку послали на 1 больше столбцов, 
			//то m+1 поток должен получить столбцы из исходной матрицы со смещением += 1
			if (amount > columnAmount)
				r++;
		}
	}
	else
	{
		//Все ненулевые потоки принимают столбцы матрицы
		MPI_Status status;
		MPI_Recv(column, matrixSize*columnAmount, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);
		/*for (int e = 0; e < matrixSize*columnAmount; e++)
		{
		cout << "rank: " << rank << " elem: " << e << " - " << column[e] << "\n";
		}
		cout << "\n";*/
	}

	//0 поток рассылает строки матрицы остальным потокам
	if (rank == 0)
	{
		int r = 0;
		for (int k = 0; k < size; k++)
		{
			int amount = rowAmount;
			if (residue > 0 && size - k <= residue)
			{
				amount++;
			}
			int *row = new int[matrixSize*amount];
			for (int l = 0; l < amount; l++)
			{
				for (int j = 0; j < matrixSize; j++)
				{
					row[j + l*matrixSize] = A[k*rowAmount + l + r][j];
				}
			}
			/*cout << "k: " << k<<"\n";
			for (int e = 0; e < matrixSize*amount; e++)
			{
			cout << row[e] << " ";
			}
			cout << "\n";*/

			//Посылаем потокам строки матрицы
			for (int m = 1; m < size; m++)
			{
				MPI_Send(row, matrixSize*amount, MPI_INT, m, k, MPI_COMM_WORLD);
			}
			if (amount > rowAmount)
				r++;
		}
	}
	else
	{
		int * elemsC = new int[columnAmount*matrixSize];
		int r = 0;
		//ненулевые потоки принимают строки матрицы
		//считают свои квадраты матрицы
		//и отправляют их нулевому потоку
		for (int k = 0; k < size; k++)
		{
			int n = rowAmount;
			if (residue > 0 && size - k <= residue)
			{
				n++;
			}
			MPI_Status rowStatus;
			int * row = new int[matrixSize*n];
			MPI_Recv(row, matrixSize*n, MPI_INT, 0, k, MPI_COMM_WORLD, &rowStatus);

			/*cout << "rank: " << rank << "\n";
			for (int e = 0; e < matrixSize*n; e++)
			{
				cout << row[e] << " ";
			}
			cout << "\n";*/

			for (int q = 0; q < columnAmount; q++)
			{
				int offset = matrixSize - rowAmount;
				for (int l = 0; l < n; l++)
				{
					int c = 0;
					for (int t = 0; t < matrixSize; t++)
					{
						c += row[t + l*matrixSize] * column[t + q*matrixSize];

					}
					//cout <<"rank: "<<rank<<" "<< c<<"\n";
					int index = l + q*rowAmount + q * offset + k*rowAmount + r;
					elemsC[index] = c;// k*rowAmount*columnAmount] = c;
				}
			}
			if (n > rowAmount)
				r++;
		}
		MPI_Send(elemsC, columnAmount*matrixSize, MPI_INT, 0, rank, MPI_COMM_WORLD);
	}


	if (rank == 0)
	{
		//0 поток считает свою часть матрицы
		for (int i = 0; i < columnAmount; i++)
		{
			for (int j = 0; j < matrixSize; j++)
			{
				int c = 0;
				for (int k = 0; k < matrixSize; k++)
				{
					c += A[j][k] * B[k][i];
				}
				C[j][i] = c;
			}
		}

		int r = 0;
		for (int m = 1; m < size; m++)
		{
			int amount = columnAmount;
			//Если есть остаток, то последним потокам посылаем на 1 больше столбцов
			if (residue > 0 && size - m <= residue)
			{
				amount++;
			}
			MPI_Status status;
			int * receivedColumns = new int[matrixSize*amount];
			MPI_Recv(receivedColumns, matrixSize*amount, MPI_INT, m, m, MPI_COMM_WORLD, &status);

			for (int l = 0; l < amount; l++)
			{
				for (int j = 0; j < matrixSize; j++)
				{
					C[j][m*columnAmount + l + r] = receivedColumns[j + l*matrixSize];
				}
			}
			if (amount > columnAmount)
				r++;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0)
	{
		t2 = MPI_Wtime();
		if (matrixSize <= 20)
		{
			cout << "\nResult matrix:\n";
			for (int i = 0; i < matrixSize; i++)
			{
				for (int j = 0; j < matrixSize; j++)
				{
					cout << C[i][j] << " ";
				}
				cout << "\n";
			}
		}
		cout << "\nTime: " << (t2 - t1);
	}

	MPI_Finalize();
}