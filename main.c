#include "median.h"
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef float number_t;

typedef struct {
  number_t *data;
  int size;
} Array;

int world_size;
int world_rank;

float masterPart(int world_size, int world_rank, int size, int partLength,
                 float *numberPart, MPI_Comm comm);

void slavePart(int world_rank, int partLength, float *numberPart, int size,
               MPI_Comm comm);

Array *array_new_random(int size) {
  Array *array = (Array *)malloc(sizeof(*array));
  array->size = size;
  array->data = (number_t *)malloc(sizeof(number_t) * size);
  int cal = 5;
  srand((cal + 1) * time(NULL));

  for (int i = 0; i < size; i++)
    array->data[i] =
        world_rank * size + i; // (number_t)(rand() - rand()) * 0.05f;
  return array;
}

void array_free(Array *array) {
  free(array->data);
  free(array);
}

int group_size(int l) { return world_size / (1 << l); }

int group_number(int l) { return world_rank / group_size(l); }

Array *array_distances_from_vp(Array *dataset, number_t vantage_point) {
  Array *distances = (Array *)malloc(sizeof(Array));
  distances->size = dataset->size;
  distances->data = (number_t *)malloc(sizeof(number_t) * distances->size);
  for (int i = 0; i < dataset->size; i++) {
    distances->data[i] = fabs(
        vantage_point - dataset->data[i]); // we calculate the distance of the
                                           // ith point from the vantage point
  }
  return distances;
}

Array *points_for_transfer(Array *distances, number_t, vantage_point,
                           number_t median, int is_process_left) {
  Array *
}

int less_than_median(Array *distances, number_t vantage_point,
                     number_t median) {
  int count_points =
      0; // counts how many points are closer to vantage point than median.
  for (int i = 0; i < distances->size; i++) {
    if (distances[i]->data < fabs(median - vantage_point)) {
      count_points++;
    }
  }
  return count_points;
}

void test_array_distances_from_vp() {
  Array *dataset = array_new_random(5);
  number_t vp = 4.56f;
  Array *distances = array_distances_from_vp(dataset, vp);
  for (int i = 0; i < dataset->size; i++) {
    printf("PROCESS %d: vantage_point = %f, point[%d] = %f and distance = %f\n",
           world_rank, vp, i, dataset->data[i], distances->data[i]);
  }
  array_free(dataset);
  /** @todo Do some actual testing with assert() */
}

int master_rank(int l) { return world_rank - (world_rank % group_size(l)); }

void test_partitioning() {
  world_size = 16;
  world_rank = 12;
  assert(master_rank(0) == 0);
  assert(master_rank(1) == 8);
  assert(master_rank(2) == 12);
}

int main(int argc, char **argv) {
  test_partitioning();
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  Array *dataset = array_new_random(5);
  int dataset_size = dataset->size * world_size;

  if (world_rank == 0) {
    float median = masterPart(world_size, world_rank, dataset_size,
                              dataset->size, dataset->data, MPI_COMM_WORLD);
    printf("Median: %.2f\n", median);
  } else {
    slavePart(world_rank, dataset->size, dataset->data, dataset_size,
              MPI_COMM_WORLD);
  }

  array_free(dataset);
  MPI_Finalize();
  return 0;
}
