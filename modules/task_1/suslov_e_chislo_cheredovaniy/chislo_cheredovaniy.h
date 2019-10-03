// Copyright 2018 Nesterov Alexander
#ifndef MODULES_TEST_TASKS_TEST_MPI_OPS_MPI_H_
#define MODULES_TEST_TASKS_TEST_MPI_OPS_MPI_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(int  sz);
int getChisloCheredovaniy(std::vector<int> global_vec, int count_size_vector);
int getParallelOperations(std::vector<int> global_vec, int count_size_vector);

#endif  // MODULES_TASK_1_CHISLO_CHEREDOVANIY_H_

