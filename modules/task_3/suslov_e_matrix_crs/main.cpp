// Copyright 2019 Suslov Egor
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./matrix_crs.h"

TEST(Multyplicate_crs_matrix, crs_matrix_row_index_test_1) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    crsMatrix A, B, C, Ctest;
    int NZ = 2; // ��������� � ������ ������
    int N = 30;
    if (rank == 0) {
        GenerateRegularCRS(N, NZ, A);
        GenerateRegularCRS(N, NZ, B);
    }
    MultiplicateMPI(A, B, C);
    if (rank == 0) {
        MultiplicateGustafson(A, B, Ctest);
        for (int i = 0; i < N + 1; i++) {
            ASSERT_EQ(C.RowIndex[i], Ctest.RowIndex[i]);
        }

    }
}

TEST(Multyplicate_crs_matrix, crs_matrix_col_test_2) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    crsMatrix A, B, C, Ctest;
    int NZ = 2;
    int N = 30;
    if (rank == 0) {
        GenerateRegularCRS(N, NZ, A);
        GenerateRegularCRS(N, NZ, B);
    }
    MultiplicateMPI(A, B, C);
    if (rank == 0) {
        MultiplicateGustafson(A, B, Ctest);
        for (int i = 0; i < C.RowIndex[N]; i++) {
            ASSERT_EQ(C.Col[i], Ctest.Col[i]);
        }

    }
}

TEST(Multyplicate_crs_matrix, crs_matrix_value_test_3) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    crsMatrix A, B, C, Ctest;
    int NZ = 2;
    int N = 30;
    if (rank == 0) {
        GenerateRegularCRS(N, NZ, A);
        GenerateRegularCRS(N, NZ, B);
    }
    MultiplicateMPI(A, B, C);
    if (rank == 0) {
        MultiplicateGustafson(A, B, Ctest);
        for (int i = 0; i < C.RowIndex[N]; i++) {
            ASSERT_EQ(C.Value[i], Ctest.Value[i]);
        }

    }
}

TEST(Multyplicate_crs_matrix, crs_matrix_one_row_index_test_4) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    crsMatrix A, B, C, Ctest;
    int NZ = 1; // ��������� � ������ ������
    int N = 40;
    if (rank == 0) {
        GenerateRegularCRS(N, NZ, A);
        GenerateRegularCRS(N, NZ, B);
    }
    MultiplicateMPI(A, B, C);
    if (rank == 0) {
        MultiplicateGustafson(A, B, Ctest);
        for (int i = 0; i < N + 1; i++) {
            ASSERT_EQ(C.RowIndex[i], Ctest.RowIndex[i]);
        }

    }
}

TEST(Multyplicate_crs_matrix, crs_matrix_one_col_test_5) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    crsMatrix A, B, C, Ctest;
    int NZ = 1;
    int N = 40;
    if (rank == 0) {
        GenerateRegularCRS(N, NZ, A);
        GenerateRegularCRS(N, NZ, B);
    }
    MultiplicateMPI(A, B, C);
    if (rank == 0) {
        MultiplicateGustafson(A, B, Ctest);
        for (int i = 0; i < C.RowIndex[N]; i++) {
            ASSERT_EQ(C.Col[i], Ctest.Col[i]);
        }

    }
}

TEST(Multyplicate_crs_matrix, crs_matrix_one_value_test_6) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    crsMatrix A, B, C, Ctest;
    int NZ = 1;
    int N = 40;
    if (rank == 0) {
        GenerateRegularCRS(N, NZ, A);
        GenerateRegularCRS(N, NZ, B);
    }
    MultiplicateMPI(A, B, C);
    if (rank == 0) {
        MultiplicateGustafson(A, B, Ctest);
        for (int i = 0; i < C.RowIndex[N]; i++) {
            ASSERT_EQ(C.Value[i], Ctest.Value[i]);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}



