// Copyright 2019 Suslov Egor
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_2/suslov_e_producer_consumer/producer_consumer.h"


std::vector<int> getRandomVector(int size_v) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vect(size_v);
    for (int i = 0; i < size_v; i++) { vect[i] = gen() % 100; }
    return vect;
}
int getPositive_elem(std::vector<int> vector, int count_size_vector) {
    if (count_size_vector < 2) {
        return 0;
    }
    int negative_elem = 0;
    for (int c = 0; c < count_size_vector; c++) {
        if (vector[c] == -1) {
            negative_elem++;
        }
    }
    return vector.size() - negative_elem;
}
int* Create_dinamic_massiv_from_vector(std::vector<int> vec) {
    int vec_size = vec.size();
    int* mas = new int[vec_size];
    for (int i = 0; i < vec_size; i++) {
        mas[i] = vec[i];
    }
    return mas;
}

int Producer(std::vector<int> &buffer_silka, std::vector<int> resurces) { //Передаем ресурсы которые нужно создать векторами

    MPI_Group MPI_GROUP_WORLD;
    MPI_Group first_row_group;
    MPI_Comm first_row_comm;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size == 1) {
        for (int i = 0; i < resurces.size(); i++) {
            buffer_silka[i] = resurces[i];
        }
        return 0;
    }
    int kol_chet_nomerov_proc;
    int chetniy_rank = 0;
    if (rank % 2 == 0) {
        chetniy_rank = 1;
    }
    MPI_Reduce(&chetniy_rank, &kol_chet_nomerov_proc, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    int* process_ranks = new int[kol_chet_nomerov_proc];
    int count = 2;
    process_ranks[0] = 0;
    for (int i = 1; i < kol_chet_nomerov_proc; i++) {
        process_ranks[i] = count;
        count += 2;
    }
    std::vector<int> buffer;
    buffer = buffer_silka;
    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
    MPI_Group_incl(MPI_GROUP_WORLD, kol_chet_nomerov_proc, process_ranks,&first_row_group);
    MPI_Comm_create(MPI_COMM_WORLD, first_row_group, &first_row_comm);
    if (rank % 2 == 1) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    int buffer_size = buffer.size();
    int one_resurces_size = resurces.size();
    int* buffer_ = Create_dinamic_massiv_from_vector(buffer); //new int[buffer_size]
    int* resurce = Create_dinamic_massiv_from_vector(resurces); //new int[one_resurces_size];
    int all_resurce_size = 0;
    int* rcount = new int[kol_chet_nomerov_proc];
    int* displs = new int[kol_chet_nomerov_proc];
    MPI_Reduce(&one_resurces_size, &all_resurce_size, 1, MPI_INT, MPI_SUM, 0, first_row_comm);
    if (rank == 0) {
        if (buffer_size < all_resurce_size) {
            for (int i = 0; i < resurces.size(); i++) {
                buffer_silka[i] = resurces[i];
            }
            return 0;
        }
    }
    MPI_Allgather(&one_resurces_size, 1, MPI_INT, &rcount, kol_chet_nomerov_proc, MPI_INT, first_row_comm);
    /*for (int i = 0; i < new_group_size; i++) {
        all_resurce_size += rcount[i];
    }*/
    int temp;
    displs[0] = 0;
    for (int i = 1; i < kol_chet_nomerov_proc; i++) {
        temp += rcount[i-1];
        displs[i] = temp;
    }
    MPI_Gatherv(resurce, one_resurces_size, MPI_INT, buffer_,  rcount, displs, MPI_INT, 0, first_row_comm);
    for (int i = 0; i < all_resurce_size; i++) {
        buffer_silka[i]  = buffer_[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

int Consumer(std::vector<int> &buffer, std::vector<int> &resurces) { //,int* buf, int size_buf, int* res, int size_res
    
    MPI_Group MPI_GROUP_WORLD;
    MPI_Group first_row_group;
    MPI_Comm first_row_comm;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size == 1) {
        for (int i = 0; i < resurces.size(); i++) {
            resurces[i] = buffer[i];
        }
        return 0;
    }
    int kol_nechet_nomerov_proc_and_null;
    int nechetniy_rank = 0;
    if (rank == 0 || rank % 2 == 1) {
        nechetniy_rank = 1;
    }
    MPI_Reduce(&nechetniy_rank, &kol_nechet_nomerov_proc_and_null, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (size % 2 == 0) { //Блокируем последний процесс если четное число процессов, т.к. для него не создали ресурс
        kol_nechet_nomerov_proc_and_null -= 1;
        if (rank = size - 1) {
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    int* process_ranks = new int[kol_nechet_nomerov_proc_and_null];
    int count = 1;
    process_ranks[0] = 0;
    for (int i = 1; i < kol_nechet_nomerov_proc_and_null; i++) {
        process_ranks[i] = count;
        count += 2;
    }
    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
    MPI_Group_incl(MPI_GROUP_WORLD, kol_nechet_nomerov_proc_and_null, process_ranks, &first_row_group);
    MPI_Comm_create(MPI_COMM_WORLD, first_row_group, &first_row_comm);
    if (rank % 2 == 0 && rank !=0) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    int buffer_size = buffer.size();
    int one_consume_resurces_size = resurces.size();
    int* buffer_ = Create_dinamic_massiv_from_vector(buffer); //new int[buffer_size]
    int* resurce = Create_dinamic_massiv_from_vector(resurces); //new int[one_resurces_size];
    int resurs_take_elem_is_bufera;
    MPI_Reduce(&one_consume_resurces_size, &resurs_take_elem_is_bufera, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    int all_resurce_size = 0;
    int* rcount = new int[kol_nechet_nomerov_proc_and_null];
    int* displs = new int[kol_nechet_nomerov_proc_and_null];
    if (buffer_size < resurces.size()) {
         for (int i = 0; i < buffer_size; i++) {
              resurces[i] = buffer[i];
              buffer[i] = -1;
         }
         return 0;
    }
    MPI_Allgather(&one_consume_resurces_size, 1, MPI_INT, &rcount, kol_nechet_nomerov_proc_and_null, MPI_INT, first_row_comm);
    /*for (int i = 0; i < new_group_size; i++) {
    all_resurce_size += rcount[i];
    }*/
    int temp;
    displs[0] = 0;
    for (int i = 1; i < kol_nechet_nomerov_proc_and_null; i++) {
        temp += rcount[i - 1];
        displs[i] = temp;
    }
    MPI_Scatterv(buffer_, rcount, displs, MPI_INT, resurce, one_consume_resurces_size, MPI_INT, 0, first_row_comm);
    for (int i = 0; i < resurs_take_elem_is_bufera; i++) {
        buffer[i] = -1;
    }
    for (int i = 0; i < one_consume_resurces_size; i++) {
        resurces[i] = resurce[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}










