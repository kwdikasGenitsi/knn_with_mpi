#pragma once
#include "dataset.h"
#include <mpi.h>

/**
 * Partitions the given dataset across a communicator, so that entries for which
 * values[i] < threshold are in the left processes and vice versa.
 */
void mpi_partition_by_value (Dataset dataset, Array values, number_t threshold,
                             MPI_Comm comm);

void test_mpi_partition_by_value ();