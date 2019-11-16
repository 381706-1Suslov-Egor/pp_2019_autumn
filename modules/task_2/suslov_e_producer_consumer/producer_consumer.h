// Copyright 2019 Suslov Egor
#ifndef MODULES_TASK_2_SUSLOV_E_PRODUCER_CONSUMER_H_
#define MODULES_TASK_2_SUSLOV_E_PRODUCER_CONSUMER_H_

#include <vector>


//int Producernew(std::vector<int> buffer, std::vector<int> vec_resursov);
int getPositive_elem(std::vector<int> vector, int count_size_vector);
std::vector<int> getRandomVector(int  sz);
int Consumer(std::vector<int> &buffer, std::vector<int> &resurces);
int Producer(std::vector<int> &buffer, std::vector<int> resurces);
//int Producer1(std::vector<int> buffer, std::vector<int> vec_resursov);
//std::vector<int> Consumers(std::vector<int> buffer, int kol_resursov);

#endif  // MODULES_TASK_2_SUSLOV_E_PRODUCER_CONSUMER_H_
