#include <stdio.h>
#include <mpi.h>
#include <math.h> 
#include <iostream>

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
			a[i][j] = i;
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

	int **A = initMatrix(matrixSize);
	int **B = initMatrix(matrixSize);
	int **C = initMatrix(matrixSize);

	if (rank == 0)
	{
		A = fillInMatrix(A, matrixSize);
		B = fillInMatrix(B, matrixSize);
		C = fillInZeroMatrix(C, matrixSize);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	int residue = matrixSize % size;
	int rowAmount = matrixSize / size;
	int columnAmount = rowAmount;
	int tempSize = size*rowAmount;
	
	if (rank == 0)
	{
		t1 = MPI_Wtime();
	}
	bool first = true;
	for (int i = 0; i < tempSize; i += rowAmount)
	{
		if (rank == 0)
		{
			//0 поток вычисляет сам первую часть матрицы
			for (int q = 0; q < rowAmount; q++)
			{
				for (int l = 0; l < columnAmount; l++)
				{
					int c = 0;
					for (int j = 0; j < matrixSize; j++)
					{
						c += A[i + q][j] * B[j][l];
					}
					C[i + q][l] = c;
				}
			}

			int r = 0;
			//Посылаем потокам части матрицы
			for (int m = 1; m < size; m++)
			{				
				int amount = columnAmount;
				if (residue > 0 && size - m <= residue)
				{
					amount++;
				}
				int *row = new int[matrixSize*rowAmount];
				int *column = new int[matrixSize*amount];
				int n;
				if (rowAmount > amount)
					n = rowAmount;
				else
					n = amount;
				for (int l = 0; l < n; l++)
				{
					for (int j = 0; j < matrixSize; j++)
					{
						if (l < rowAmount)
							row[j + l*matrixSize] = A[i + l][j];
						if (l < amount)
							column[j + l*matrixSize] = B[j][m*columnAmount + l + r];
					}
				}
				MPI_Send(row, matrixSize*rowAmount, MPI_INT, m, m, MPI_COMM_WORLD);
				MPI_Send(column, matrixSize*amount, MPI_INT, m, m, MPI_COMM_WORLD);


				//Принимаем подсчитанные элементы матрицы от других потоков
				//и заполняем ими конечную матрицу
				MPI_Status status;
				int * elemsC = new int[rowAmount*amount];
				MPI_Recv(elemsC, rowAmount*amount, MPI_INT, m, m, MPI_COMM_WORLD, &status);
				for (int q = 0; q < rowAmount; q++)
				{
					for (int l = 0; l < amount; l++)
					{
						C[i + q][m*columnAmount + l + r] = elemsC[l + amount*q];
					}
				}
				if (amount > columnAmount)
					r++;
			}
		}
		else
		{
			MPI_Status rowStatus, columnStatus;
			int amount = columnAmount;
			if (residue > 0 && size - rank <= residue)
			{
				amount++;
			}
			int * row = new int[matrixSize*rowAmount];
			int *column = new int[matrixSize*amount];
			MPI_Recv(row, matrixSize*rowAmount, MPI_INT, 0, rank, MPI_COMM_WORLD, &rowStatus);
			MPI_Recv(column, matrixSize*amount, MPI_INT, 0, rank, MPI_COMM_WORLD, &columnStatus);

			int * elemsC = new int[rowAmount*amount];
			for (int q = 0; q < rowAmount; q++)
			{
				for (int l = 0; l < amount; l++)
				{
					int c = 0;
					for (int k = 0; k < matrixSize; k++)
					{
						c += row[k + q*matrixSize] * column[k + l*matrixSize];
					}
					elemsC[l + q*amount] = c;
				}
			}
			MPI_Send(elemsC, rowAmount*amount, MPI_INT, 0, rank, MPI_COMM_WORLD);
		}

		if (first && residue > 0 && i + rowAmount >= tempSize)
		{
			i = i + rowAmount- residue;
			tempSize = matrixSize;
			rowAmount = residue;
			first = false;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0)
	{
		t2 = MPI_Wtime();
		if (matrixSize < 20)
		{
			for (int i = 0; i < matrixSize; i++)
			{
				for (int j = 0; j < matrixSize; j++)
				{
					cout << C[i][j] << " ";
				}
				cout << "\n";
			}
		}
		cout << "Time: " << (t2 - t1);
	}

	MPI_Finalize();
}