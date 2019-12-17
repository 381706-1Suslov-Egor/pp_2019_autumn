// Copyright 2019 Suslov Egor
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iostream>
#include "../../../modules/task_3/suslov_e_matrix_crs/matrix_crs.h"

bool isSrandCalled = false;

void InitializeMatrix(int N, int NZ, crsMatrix& mtx)
{
    mtx.N = N;
    mtx.NZ = NZ;
    mtx.Value = new double[NZ];
    mtx.Col = new int[NZ];
    mtx.RowIndex = new int[N + 1];
}

void FreeMatrix(crsMatrix& mtx)
{
    delete[] mtx.Value;
    delete[] mtx.Col;
    delete[] mtx.RowIndex;
}

void GenerateRegularCRS(int N, int cntInRow, crsMatrix& mtx)
{
    std::mt19937 mersenne;
    mersenne.seed(static_cast<unsigned int>(time(0)));
    int i, j, k, f, tmp, notNull, c;
    notNull = cntInRow * N;
    InitializeMatrix(N, notNull, mtx);
    for (i = 0; i < N; i++) {
        for (j = 0; j < cntInRow; j++) { // Формируем номера столбцов в строке i
            do {
                mtx.Col[i * cntInRow + j] = mersenne() % N;
                f = 0;
                for (k = 0; k < j; k++)
                    if (mtx.Col[i * cntInRow + j] ==
                        mtx.Col[i * cntInRow + k])
                        f = 1;
            } while (f == 1);
        }
        for (j = 0; j < cntInRow - 1; j++) // Сортируем номера столцов в строке i
            for (k = 0; k < cntInRow - 1; k++)
                if (mtx.Col[i * cntInRow + k] > mtx.Col[i * cntInRow + k + 1]) {
                    tmp = mtx.Col[i * cntInRow + k];
                    mtx.Col[i * cntInRow + k] =
                        mtx.Col[i * cntInRow + k + 1];
                    mtx.Col[i * cntInRow + k + 1] = tmp;
                }
    }
    for (i = 0; i < cntInRow * N; i++) // Заполняем массив значений
        mtx.Value[i] = (double(mersenne() % 1000)/100);
    c = 0; // Заполняем массив индексов строк
    for (i = 0; i <= N; i++) {
        mtx.RowIndex[i] = c;
        c += cntInRow;
    }
}

int Multiplicatew(crsMatrix A, crsMatrix B, crsMatrix& C, double& time)
{
    if (A.N != B.N)
        return 1;
    int N = A.N;
    int i, j, k;
    std::vector<int> columns;
    std::vector<double> values;
    std::vector<int> row_index;
    clock_t start = clock();
    int NZ = 0;
    int* temp = new int[N];
    row_index.push_back(0);
    for (i = 0; i < N; i++)
    {
        memset(temp, -1, N * sizeof(int)); // i-я строка матрицы A, Обнуляем массив указателей на элементы
        int ind1 = A.RowIndex[i], ind2 = A.RowIndex[i + 1]; // массив указателей Идем по ненулевым элементам строки и заполняем
        for (j = ind1; j < ind2; j++) {
            int col = A.Col[j];
            temp[col] = j; 
        }
        for (j = 0; j < N; j++) { // j-я строка матрицы B
            double sum = 0;   
            int ind3 = B.RowIndex[j], ind4 = B.RowIndex[j + 1]; // Все ненулевые элементы строки j матрицы B
            for (k = ind3; k < ind4; k++) {
                int bcol = B.Col[k];
                int aind = temp[bcol];
                    if (aind != -1)
                        sum += A.Value[aind] * B.Value[k];
            }
            if (fabs(sum) > 0.000001) {
                columns.push_back(j);
                values.push_back(sum);
                NZ++;
            }
        }
        row_index.push_back(NZ);
    }
    InitializeMatrix(N, NZ, C);
    for (j = 0; j < NZ; j++) {
        C.Col[j] = columns[j];
        C.Value[j] = values[j];
    }
    for (i = 0; i <= N; i++)
        C.RowIndex[i] = row_index[i];
    delete[] temp;
    clock_t finish = clock();
    time = (double)(finish - start) / CLOCKS_PER_SEC;
    return 0;
}

int Create_part_crs(double* polnaya_stroka_C, int &N, double* &crs_C_value, int* &crs_C_col)
{
    int push_elem = 0;
    std::vector<double> vec_value;
    std::vector<int> vec_col;
    for (int k = 0; k < N; k++)
    {
        if (polnaya_stroka_C[k] != 0)
            vec_value.push_back(polnaya_stroka_C[k]);
            vec_col.push_back(k);
    }
    push_elem = int(vec_col.size());
    crs_C_value = new double[push_elem];
    crs_C_col = new int[push_elem];
    for (int i = 0; i < push_elem; i++) {
        crs_C_value[i] = vec_value.at(i); //vec_value[i];
        crs_C_col[i] = vec_col.at(i); //vec_col[i];
    }
    return push_elem;
}

void MultiplicateGustafson(crsMatrix A, crsMatrix B, crsMatrix& C)
{
    int add_elem_vsego = 0, add_elem = 0;
    int strok_peredali_is_mtrA = A.N;
    std::vector<double> polnaya_stroka_C(A.N, 0);
    std::vector<double> crs_C_value_full;
    std::vector<int> crs_C_col_full;
    std::vector<int> crs_C_row_index_full;
    crs_C_row_index_full.push_back(0);
    for (int i = 0; i < strok_peredali_is_mtrA; i++) { //i-ая строка матрицы А
        for (int j = A.RowIndex[i]; j < A.RowIndex[i + 1]; j++) //Элементы i-ой строки матр А
        {
            int Col_elem_A = A.Col[j];
            for (int k = B.RowIndex[Col_elem_A]; k < B.RowIndex[Col_elem_A + 1]; k++) //Элементы B.Col[j]-ой строки матр В (Соответсвующей столбцу j эл-та из А)
            {
                polnaya_stroka_C[B.Col[k]] += A.Value[j] * B.Value[k];
            }
        }
        for (int k = 0; k < A.N; k++)
        {
            if (polnaya_stroka_C[k] != 0) {
                crs_C_value_full.push_back(polnaya_stroka_C[k]);
                crs_C_col_full.push_back(k);
                polnaya_stroka_C[k] = 0;
            }
        }
        crs_C_row_index_full.push_back(crs_C_value_full.size()); //заполняем rowindex
    }
    int row_index_size = crs_C_row_index_full.size(); //row_index_size = количество заполненных строк
    C.N = A.N; // кол-во обр строк A, С
    add_elem_vsego = int(crs_C_value_full.size());
    C.NZ = add_elem_vsego;
    C.Col = new int[add_elem_vsego];
    C.Value = new double[add_elem_vsego];
    C.RowIndex = new int[row_index_size]; //начинается с нуля
    for (int i = 0; i < add_elem_vsego; i++)
    {
        C.Col[i] = crs_C_col_full.at(i);
        C.Value[i] = crs_C_value_full.at(i);
    }
    for (int i = 0; i < row_index_size; i++)
    {
        C.RowIndex[i] = crs_C_row_index_full.at(i);
    }
}

int MultiplicateMPI(crsMatrix& A, crsMatrix& B, crsMatrix& C)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size == 1) {
        MultiplicateGustafson(A, B, C);
        return 0;
    }
    int* upakovka;
    upakovka = new int[3];
    int row_on_proc, ostatok_last_proc, vsego_elem_mtr_B;
    if (rank == 0) {
        C.N = A.N;
        row_on_proc = A.N / size;
        ostatok_last_proc = A.N % size;
        vsego_elem_mtr_B = B.RowIndex[B.N];
        upakovka[0] = A.N;
        upakovka[1] = vsego_elem_mtr_B;
        upakovka[2] = row_on_proc;
    } else {

    }
    MPI_Bcast(&upakovka[0], 3, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        row_on_proc = upakovka[2];
        vsego_elem_mtr_B = upakovka[1];
        int N = upakovka[0];
        A.N = N;
        B.N = N;
        int size_row_index = N+1;
        A.RowIndex = new int[row_on_proc+1];
        B.Col = new int[vsego_elem_mtr_B];
        B.Value = new double[vsego_elem_mtr_B];
        B.RowIndex = new int[size_row_index];
    }
    MPI_Bcast(B.Col, vsego_elem_mtr_B, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(B.Value, vsego_elem_mtr_B, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B.RowIndex, B.N+1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            MPI_Send(&A.RowIndex[row_on_proc * i + ostatok_last_proc], row_on_proc + 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&A.Col[A.RowIndex[row_on_proc * i + ostatok_last_proc]], A.RowIndex[row_on_proc * (i + 1) + ostatok_last_proc] - A.RowIndex[row_on_proc * i + ostatok_last_proc], MPI_INT, i, 1, MPI_COMM_WORLD); // отправляем количество элементов содержащееся в данных строках
            MPI_Send(&A.Value[A.RowIndex[row_on_proc * i + ostatok_last_proc]], A.RowIndex[row_on_proc * (i+1) + ostatok_last_proc] - A.RowIndex[row_on_proc * i + ostatok_last_proc], MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
        }
        crsMatrix C_temp;
        create_part_crs_C(row_on_proc + ostatok_last_proc, A, B, C_temp);
        int* massiv_elem_proc_helperov;
        massiv_elem_proc_helperov = new int[size];
        massiv_elem_proc_helperov[0] = C_temp.RowIndex[ostatok_last_proc + row_on_proc - 1]; // первый процесс добавил элементов
        MPI_Status st;
        for (int i = 1; i < size; i++) {
            MPI_Recv(&massiv_elem_proc_helperov[i], 1, MPI_INT, i, 5, MPI_COMM_WORLD, &st); // принимаем количество элементов созданных другими процессами
        }
        int vsego_elem_so_vseh_proc = 0;
        for (int i = 0; i < size; i++) {
            vsego_elem_so_vseh_proc += massiv_elem_proc_helperov[i];
        }
        C.RowIndex = new int[A.N + 1];
        C.Col = new int[vsego_elem_so_vseh_proc];
        C.Value = new double[vsego_elem_so_vseh_proc];
        
        C.RowIndex[0] = 0;
        for (int i = 0; i < massiv_elem_proc_helperov[0]; i++)
        {
            C.Col[i] = C_temp.Col[i];
            C.Value[i] = C_temp.Value[i];
        }
        for (int i = 1; i <= ostatok_last_proc + row_on_proc; i++)
        {
            C.RowIndex[i] = C_temp.RowIndex[i-1];
        }
        MPI_Status status1, status2, status3;
        for (int i = 1; i < size; i++) {
            MPI_Recv(&C.RowIndex[row_on_proc * i + ostatok_last_proc + 1], row_on_proc, MPI_INT, i, 1, MPI_COMM_WORLD, &status1);
            for (int j = 0; j < row_on_proc; j++) {
                C.RowIndex[row_on_proc * i + ostatok_last_proc+1+j] += C.RowIndex[row_on_proc * i + ostatok_last_proc];
            }
            MPI_Recv(&C.Col[C.RowIndex[row_on_proc * i + ostatok_last_proc]], massiv_elem_proc_helperov[i], MPI_INT, i, 2, MPI_COMM_WORLD, &status2);
            MPI_Recv(&C.Value[C.RowIndex[row_on_proc * i + ostatok_last_proc]], massiv_elem_proc_helperov[i], MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status3);
        }
        FreeMatrix(C_temp);
        delete[] massiv_elem_proc_helperov;
    } else {
        MPI_Recv(&A.RowIndex[0], row_on_proc + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        int sdvig_strok_k_starty = A.RowIndex[0];
        for (int i = 0; i < row_on_proc + 1; i++) { // вычитаем начальный элемент чтобы row_index характеризовал value и col
            A.RowIndex[i] -= sdvig_strok_k_starty;
        }
        A.Col = new int[A.RowIndex[row_on_proc]];
        A.Value = new double[A.RowIndex[row_on_proc]];
        MPI_Status status1, status2;
        MPI_Recv(&A.Col[0], A.RowIndex[row_on_proc], MPI_INT, 0, 1, MPI_COMM_WORLD, &status1);
        MPI_Recv(&A.Value[0], A.RowIndex[row_on_proc], MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status2);
        create_part_crs_C(row_on_proc, A, B, C);
        int sozdali_elem = C.RowIndex[row_on_proc - 1]; //-1 т.к. функция удаляет нулевой элемент
        MPI_Send(&sozdali_elem, 1, MPI_INT, 0, 5, MPI_COMM_WORLD); //Посылаем кол-во элементов созданное данным процессом
       
        MPI_Send(&C.RowIndex[0], row_on_proc, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&C.Col[0], sozdali_elem, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&C.Value[0], sozdali_elem, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        FreeMatrix(A);
        FreeMatrix(B);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank != 0) {
        FreeMatrix(C);
    }
    return 0;
}
void create_part_crs_C(int row_peredali, crsMatrix& A, crsMatrix& B, crsMatrix& C)
{
    int add_elem_vsego = 0, add_elem = 0;
    int strok_peredali_is_mtrA = row_peredali;
    std::vector<double> polnaya_stroka_C(A.N, 0);
    std::vector<double> crs_C_value_full;
    std::vector<int> crs_C_col_full;
    std::vector<int> crs_C_row_index_full;
    for (int i = 0; i < strok_peredali_is_mtrA; i++) { //i-ая строка матрицы А
        for (int j = A.RowIndex[i]; j < A.RowIndex[i + 1]; j++) { //Элементы i-ой строки матр А
            int Col_elem_A = A.Col[j];
            for (int k = B.RowIndex[Col_elem_A]; k < B.RowIndex[Col_elem_A + 1]; k++) { //Элементы B.Col[j]-ой строки матр В (Соответсвующей столбцу j эл-та из А)
                polnaya_stroka_C[B.Col[k]] += A.Value[j] * B.Value[k];
            }
        }
        for (int k = 0; k < A.N; k++) {
            if (polnaya_stroka_C[k] != 0) {
                crs_C_value_full.push_back(polnaya_stroka_C[k]);
                crs_C_col_full.push_back(k);
                polnaya_stroka_C[k] = 0;
            }
        }
        crs_C_row_index_full.push_back(crs_C_value_full.size()); //заполняем rowindex
    }
    int row_index_size = crs_C_row_index_full.size(); //row_index_size = количество заполненных строк
    C.N = A.N; // кол-во обр строк A, С
    add_elem_vsego = crs_C_value_full.size();
    //C.NZ = add_elem_vsego;
    C.Col = new int[add_elem_vsego];
    C.Value = new double[add_elem_vsego];
    C.RowIndex = new int[row_index_size]; //не начинается с нуля
    for (int i = 0; i < add_elem_vsego; i++) {
        C.Col[i] = crs_C_col_full.at(i);
        C.Value[i] = crs_C_value_full.at(i);
    }
    for (int i = 0; i < row_index_size; i++) {
        C.RowIndex[i] = crs_C_row_index_full.at(i);
    }
}

double** create_norm_mtr(crsMatrix A)
{
    double** norm_mtr;
    norm_mtr = new double*[A.N];
    for (int j = 0; j < A.N; j++) {
        norm_mtr[j] = new double[A.N];
    }
    for (int j = 0; j < A.N; j++) {
        for (int i = 0; i < A.N; i++) {
            norm_mtr[i][j] = 0;
        }
    }
    for (int i = 0; i < A.N; i++) {
        for (int j = A.RowIndex[i]; j < A.RowIndex[i + 1]; j++) {
        
            norm_mtr[i][A.Col[j]] = A.Value[j];
        }
    }
    return norm_mtr;
}
void print_norm_mtr(double** norm_mtr, int N) {
    std::cout << std::endl;
    std::cout << std::endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << norm_mtr[i][j]<<"   ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
}
double** mult_norm_matr(double** A, double** B, int N) {
    double** C;
    C = new double* [N];
    for (int j = 0; j < N; j++) {
        C[j] = new double[N];
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0;
            for (int k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
        std::cout << std::endl;
    }
    return C;
}
