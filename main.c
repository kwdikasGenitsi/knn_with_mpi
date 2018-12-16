#include "median.h"
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef float number_t;

typedef struct
{
    number_t* data;
    int size;
} Array;

int world_size;
int world_rank;

/**
 * To know the depth of our tree, we need to calculate log2(x) accurately.
 * We can't use floating point numbers since we risk getting an off-by-one error
 * due to floating point accuracy.
 * We therefore use this GCC bultin, which returns the number of leading 0s in
 * an integer. On x86 this is done with the BSR instruction.
 */
static inline int
log2i(int x)
{
    assert(x > 0);
    return sizeof(int) * 8 - __builtin_clz(x) - 1;
}

void
test_log2i()
{
    int n = 1;
    for (int i = 0; i < 15; i++) {
        assert(log2i(n) == i);
        n *= 2;
    }
}

float
masterPart(int world_size,
           int world_rank,
           int size,
           int partLength,
           float* numberPart,
           MPI_Comm comm);

void
slavePart(int world_rank,
          int partLength,
          float* numberPart,
          int size,
          MPI_Comm comm);

Array*
array_new_random(int size)
{
    Array* array = (Array*) malloc(sizeof(*array));
    array->size = size;
    array->data = (number_t*) malloc(sizeof(number_t) * size);
    int cal = 5;
    srand((cal + 1) * time(NULL));

    for (int i = 0; i < size; i++)
        array->data[i]
          = world_rank * size + i; // (number_t)(rand() - rand()) * 0.05f;
    return array;
}

void
array_free(Array* array)
{
    free(array->data);
    free(array);
}

int
group_size(int l)
{
    return world_size / (1 << l);
}

int
group_number(int l)
{
    return world_rank / group_size(l);
}

Array*
array_distances_from_vp(Array* dataset, number_t vantage_point)
{
    Array* distances = (Array*) malloc(sizeof(Array));
    distances->size = dataset->size;
    distances->data = (number_t*) malloc(sizeof(number_t) * distances->size);
    for (int i = 0; i < dataset->size; i++) {
        distances->data[i] = fabs(
          vantage_point - dataset->data[i]); // we calculate the distance of the
                                             // ith point from the vantage point
    }
    return distances;
}

void
test_array_distances_from_vp()
{
    Array* dataset = array_new_random(5);
    number_t vp = 4.56f;
    Array* distances = array_distances_from_vp(dataset, vp);
    for (int i = 0; i < dataset->size; i++) {
        printf(
          "PROCESS %d: vantage_point = %f, point[%d] = %f and distance = %f\n",
          world_rank, vp, i, dataset->data[i], distances->data[i]);
    }
    array_free(dataset);
    /** @todo Do some actual testing with assert() */
}

int
less_than_median(Array* array, number_t vantage_point)
{
    int count_points
      = 0; // counts how many points are closer to vantage point than median.
    return count_points;
}

int
master_rank(int l)
{
    return world_rank - (world_rank % group_size(l));
}

void
test_partitioning()
{
    world_size = 16;
    world_rank = 12;
    assert(master_rank(0) == 0);
    assert(master_rank(1) == 8);
    assert(master_rank(2) == 12);
}

MPI_Comm
communicator_for_level(int l)
{
    MPI_Comm c;
    MPI_Comm_split(MPI_COMM_WORLD, group_number(l), 0, &c);
    return c;
}

void
split_parallel(int l)
{
    MPI_Comm comm = communicator_for_level(l);
    MPI_Comm_free(&comm);
}

int
main(int argc, char** argv)
{
    test_partitioning();
    test_log2i();

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    Array* dataset = array_new_random(5);
    int dataset_size = dataset->size * world_size;

#if 0
    for (int l = 0; group_size(l) > 1; l++) {
        split_parallel(l);
    }
#endif

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
