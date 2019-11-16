// Copyright 2019 Suslov Egor
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./producer_consumer.h"

TEST(Producer_Consumer, Producer_Test1) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int kol_resursov = 2;
    std::vector<int> buffer(100,-1);
    std::vector<int> *A = &buffer;
    std::vector<int> vec_resursov;
    vec_resursov = getRandomVector(kol_resursov);
    Producer(*A, vec_resursov);
    if (size == 1) {
        if (rank == 0) {
            ASSERT_EQ(getPositive_elem(buffer, buffer.size()), kol_resursov);
        }
    }
    else {
        int vec_resursov_size = vec_resursov.size();
        for (int i = 0; i < vec_resursov_size; i++) {
            ASSERT_EQ(buffer[i], vec_resursov[i]);
        }
    }
}
TEST(Producer_Consumer, Producer_Test2) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int kol_resursov = 2;
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    std::vector<int> vec_resursov;
    vec_resursov = getRandomVector(kol_resursov);
    Producer(*A, vec_resursov);
    if (rank == 0) {
        ASSERT_EQ(getPositive_elem(buffer, buffer.size()), kol_resursov);
    }
}
TEST(Producer_Consumer, Producer_Test3) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int kol_resursov = 2;
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    std::vector<int> vec_resursov;
    if (rank % 2 == 0) {
        vec_resursov = getRandomVector(kol_resursov);
    }
    Producer(*A, vec_resursov);
    if (size == 1) {
        if (rank == 0) {
            int vec_resursov_size = vec_resursov.size();
            for (int i = 0; i < vec_resursov_size; i++) {
                ASSERT_EQ(buffer[i], vec_resursov[i]);
            }
        }
    }
    else {
        if (rank == 0) {
            ASSERT_EQ(getPositive_elem(buffer, buffer.size()), size*2);
        }
    }
}
TEST(Producer_Consumer, Consumer_Test1) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int kol_resursov = 2;
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    std::vector<int> vec_resursov_produce;
    std::vector<int> vec_resursov_consume;
    if (rank % 2 == 0) {
        vec_resursov_produce = getRandomVector(kol_resursov);
    }
    Producer(*A, vec_resursov_produce);
    if (rank % 2 == 1 || rank == 0) {
        std::vector<int> vec_resursov_consume(kol_resursov, -1);
    }
    std::vector<int> *B = &vec_resursov_consume;
    Consumer(*A, *B);
    if (size == 1) {
        if (rank == 0) {
            int vec_resursov_consume_size = vec_resursov_consume.size();
            for (int i = 0; i < vec_resursov_consume_size; i++) {
                ASSERT_EQ(vec_resursov_consume[i], vec_resursov_produce[i]);
            }
        }
    }
    else {
        if (rank == 1) {
            ASSERT_EQ(getPositive_elem(buffer, buffer.size()), size * 2);
        }
    }
}
TEST(Producer_Consumer, Consumer_Test2) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int kol_resursov = 2;
    std::vector<int> buffer(100, -1);
    std::vector<int> *A = &buffer;
    std::vector<int> vec_resursov_produce;
    std::vector<int> vec_resursov_consume;
    if (rank % 2 == 0) {
        vec_resursov_produce = getRandomVector(kol_resursov);
    }
    Producer(*A, vec_resursov_produce);
    if (rank % 2 == 1 || rank == 0) {
        std::vector<int> vec_resursov_consume(kol_resursov,-1);
    }
    std::vector<int> *B = &vec_resursov_consume;
    Consumer(*A, *B);
    if (size == 1) {
        if (rank == 0) {
            int vec_resursov_consume_size = vec_resursov_consume.size();
            for (int i = 0; i < vec_resursov_consume_size; i++) {
                ASSERT_EQ(vec_resursov_produce[i], vec_resursov_consume[i]);
            }
        }
    }
    else {
        if (size % 2 == 1) {

            if (rank == 1) {
                ASSERT_EQ(getPositive_elem(vec_resursov_consume, vec_resursov_consume.size()), kol_resursov);
            }
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

