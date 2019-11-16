// Copyright 2019 Suslov Egor
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./producer_consumer.h"

TEST(Producer_Consumer, Producer_Test1) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    int kol_resursov = 5;
    if(size == 1) {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, rank, 1);
        }
        for (int i = 0; i < kol_resursov; i++) {
            ASSERT_EQ(buffer[i], 1);
        }
    } else {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, kol_resursov %(size-1), 2);
        }
        for (int i = 0; i < kol_resursov; i++) {
            ASSERT_EQ(buffer[i], 2);
        }
    }
}

TEST(Producer_Consumer, Producer_Test2) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    int kol_resursov = 10;
    if (size == 1) {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, rank, 1);
        }
        for (int i = 0; i < kol_resursov; i++) {
            ASSERT_EQ(buffer[i], 1);
        }
    } else {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, kol_resursov % (size - 1), 2);
        }
        for (int i = 0; i < kol_resursov; i++) {
            ASSERT_EQ(buffer[i], 2);
        }
    }
}

TEST(Producer_Consumer, Producer_Test3) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    int kol_resursov = 5;
    if (size == 1) {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, rank, 1);
        }
        for (int i = 0; i < kol_resursov; i++) {
            ASSERT_EQ(buffer[i], 1);
        }
    } else {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, kol_resursov % (size - 1), 2);
        }
        for (int i = 0; i < kol_resursov; i++) {
            ASSERT_EQ(buffer[i], 2);
        }
    }
}

TEST(Producer_Consumer, Consumer_Test1) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    int kol_resursov = 5;
    if (size == 1) {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, rank, 1);
        }
        std::vector<int> resurce_consume(kol_resursov, -1);
        for (int i = 0; i < kol_resursov; i++) {
            Consumer(*A, rank, resurce_consume[i]);
            ASSERT_EQ(1, resurce_consume[i]);
        }
    } else {
        std::vector<int> resurce_consume(kol_resursov, -1);
        for (int i = 0; i < kol_resursov; i++) {
            Consumer(*A, 0, resurce_consume[i]);
            ASSERT_EQ(1, resurce_consume[i]);
        }
        
    }
}

TEST(Producer_Consumer, Consumer_Test2) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    int kol_resursov = 10;
    if (size == 1) {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, rank, 1);
        }
        std::vector<int> resurce_consume(kol_resursov, -1);
        for (int i = 0; i < kol_resursov; i++) {
            Consumer(*A, 0, resurce_consume[i]);
            ASSERT_EQ(1, resurce_consume[i]);
        }
    } else {
        if (rank == 1) {
            for (int i = 0; i < kol_resursov; i++) {
                Producer(*A, rank, 1);
            }
        }
        if (rank == 0) {
            std::vector<int> resurce_consume(kol_resursov, -1);
            for (int i = 0; i < kol_resursov; i++) {
                Consumer(*A, 0, resurce_consume[i]);
                ASSERT_EQ(1, resurce_consume[i]);
            }
        }

    }
}
TEST(Producer_Consumer, Consumer_Test3) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    int kol_resursov = 5;
    if (size == 1) {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, rank, 1);
        }
        std::vector<int> resurce_consume(kol_resursov, -1);
        for (int i = 0; i < kol_resursov; i++) {
            Consumer(*A, 0, resurce_consume[i]);
            ASSERT_EQ(1, resurce_consume[i]);
        }
    } else {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, kol_resursov % (size - 1), 2);
        }
        for (int i = 0; i < kol_resursov; i++) {
            ASSERT_EQ(buffer[i], 2);
        }
    }
}
TEST(Producer_Consumer, Producer_Consumer_Test1) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    int kol_resursov = 5;
    if (size == 1) {
        for (int i = 0; i < kol_resursov; i++) {
            Producer(*A, rank, 3);
        }
        std::vector<int> resurce_consume(kol_resursov, -1);
        for (int i = 0; i < kol_resursov; i++) {
            Consumer(*A, 0, resurce_consume[i]);
            ASSERT_EQ(3, resurce_consume[i]);
        }
    } else {
        Producer(*A, rank, 3);
        for (int i = 0; i < size; i++) {
            ASSERT_EQ(buffer[i], 3);
        }
        std::vector<int> resurce_consume(kol_resursov, -1);
        for (int i = 0; i < size; i++) {
            Consumer(*A, 0, resurce_consume[i]);
            ASSERT_EQ(3, resurce_consume[i]);
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

