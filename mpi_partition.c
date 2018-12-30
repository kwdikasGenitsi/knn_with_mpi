#include "mpi_partition.h"
#include "stack.h"
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

void
mpi_partition_by_value (Dataset dataset, Array values, number_t threshold,
                        MPI_Comm comm)
{
  int rank, comm_size;
  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &comm_size);
  assert (comm_size > 1);

  size_t stride
    = dataset.feature_count + 1; /* Useful for navigating the stacks. */

  bool is_left = rank < (comm_size / 2) ? true : false;

  /* Create a buffer and fill it with the values that need to be sent. */
  Stack *stack = stack_new (dataset.data.size);

  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t value = values.data[i];
      if ((is_left && (value > threshold))
          || (!is_left && (value <= threshold)))
        {
          number_t *point = dataset_point (dataset, i);
          stack_push_buffer (stack, point, stride);
          /* Mark the empty slots by setting their first coordinate to
           * the threshold. */
          point[1] = threshold;
        }
    }

  size_t elements_to_swap = stack_get_size (stack);
  /* Stores the number of items each process wants to swap.
   * Only valid for the master.
   */
  size_t elements[comm_size];
  MPI_Gather (&elements_to_swap, 1, MPI_UINT64_T, elements, 1, MPI_UINT64_T, 0,
              comm);

  if (rank == 0)
    {
      printf ("Elements: ");
      for (int i = 0; i < comm_size; i++)
        {
          printf ("[%zu]", elements[i] / (dataset.feature_count + 1));
        }
      printf ("\n");
    }

  /* Print the stacks to send. */
  for (int i = 0; i < comm_size; i++)
    {
      MPI_Barrier (comm);
      if (i == rank)
        {
          printf ("Rank %d: ", i);
          stack_dump (stack);
          printf ("\n");
        }
    }
  MPI_Barrier (comm);

  Array receive_buffer = array_new (elements_to_swap);
  if (rank == 0)
    {
      Stack *buffer_stack = stack_new (dataset.data.size * 4);
      int middle_rank = comm_size / 2;
      int sending_rank = middle_rank;

      for (int i = 0; i < middle_rank; i++)
        {
          /* How many elements does this process need? */
          size_t requested_elements = elements[i];

          while (stack_get_size (buffer_stack) < requested_elements)
            {
              if (sending_rank >= comm_size)
                fprintf (stderr, "Imbalanced request going left.\n");

              size_t provided_elements = elements[sending_rank];
              Array buffer = array_new (provided_elements);

              MPI_Recv (buffer.data, provided_elements, MPI_NUMBER_T,
                        sending_rank, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

              stack_push_array (buffer_stack, buffer);
              sending_rank++;

              array_free (buffer);
            }

          /* Send the elements to the process. */
          Array buffer = array_new (requested_elements);

          stack_pop_array (buffer_stack, buffer, 0, requested_elements);
          if (i == 0)
            memcpy (receive_buffer.data, buffer.data,
                    sizeof (number_t) * requested_elements);
          else
            MPI_Send (buffer.data, requested_elements, MPI_NUMBER_T, i, 0,
                      comm);

          array_free (buffer);
        }

      /* The stack should be empty by this point. */
      assert (stack_get_size (buffer_stack) == 0);

      sending_rank = 0;
      for (int i = middle_rank; i < comm_size; i++)
        {
          /* How many elements does this process need? */
          size_t requested_elements = elements[i];

          while (stack_get_size (buffer_stack) < requested_elements)
            {
              if (sending_rank >= middle_rank)
                fprintf (stderr, "Imbalanced request going right.\n");

              size_t provided_elements = elements[sending_rank];
              Array buffer = array_new (provided_elements);

              if (sending_rank == 0)
                {
                  memcpy (buffer.data, stack->stack,
                          provided_elements * sizeof (number_t));
                  stack->top = 0; /* The data has been sent, empty the stack. */
                }
              else
                MPI_Recv (buffer.data, provided_elements, MPI_NUMBER_T,
                          sending_rank, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

              stack_push_array (buffer_stack, buffer);
              sending_rank++;

              array_free (buffer);
            }

          /* Send the elements to the process. */
          Array buffer = array_new (requested_elements);

          stack_pop_array (buffer_stack, buffer, 0, requested_elements);
          MPI_Send (buffer.data, requested_elements, MPI_NUMBER_T, i, 0, comm);

          array_free (buffer);
        }

      /* The stack should once again be empty. */
      assert (stack_get_size (buffer_stack) == 0);

      stack_free (buffer_stack);
    }
  else
    {
      MPI_Status status;

      /* Do a send/receive so that the master can coordinate the rest.
       */
      MPI_Sendrecv (stack->stack, elements_to_swap, MPI_NUMBER_T, 0, 0,
                    receive_buffer.data, elements_to_swap, MPI_NUMBER_T, 0,
                    MPI_ANY_TAG, comm, &status);

      int actual_elements_received = -1;
      MPI_Get_count (&status, MPI_NUMBER_T, &actual_elements_received);
      assert ((int) elements_to_swap == actual_elements_received);

      /* Empty the stack, since the data has been sent. */
      stack->top = 0;
    }
  stack_push_array (stack, receive_buffer);
  array_free (receive_buffer);

  /* All MPI communication should be finished by this point. */
  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      if (point[1] == threshold)
        {
          stack_pop_buffer (stack, point, stride);
        }
    }
  assert (stack_get_size (stack) == 0);
  stack_free (stack);
}

void
test_mpi_partition_by_value ()
{
  Dataset dataset = dataset_new (1, 1024);
  int rank, comm_size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &comm_size);

  Array values = array_new (dataset.size);

  number_t local_sum = 0.0f;
  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      point[0] = (number_t) (i + rank * dataset.size); /* Random ID. */
      /* Should produce an equal number of 0s and 1s. */
      point[1] = (number_t) (i + (comm_size - rank - 1) * dataset.size);
      values.data[i] = point[1];
      local_sum += point[1];
    }

  /* Sum all the elements. */
  number_t sum = 0.0f;
  MPI_Reduce (&local_sum, &sum, 1, MPI_NUMBER_T, MPI_SUM, 0, MPI_COMM_WORLD);

  number_t median = (dataset.size * comm_size - 1) / 2.0f;

  mpi_partition_by_value (dataset, values, median, MPI_COMM_WORLD);

  /* Check that the sum is the same. */
  number_t new_sum = 0.0f;
  local_sum = 0.0f;
  for (size_t i = 0; i < dataset.size; i++)
    local_sum += dataset_point (dataset, i)[1];
  MPI_Reduce (&local_sum, &new_sum, 1, MPI_NUMBER_T, MPI_SUM, 0,
              MPI_COMM_WORLD);

  assert (new_sum == sum);

  bool is_left = rank < (comm_size / 2) ? true : false;

  for (size_t i = 0; i < dataset.size; i++)
    {
      number_t *point = dataset_point (dataset, i);
      if (is_left)
        {
          assert (point[1] <= median);
        }
      else
        {
          assert (point[1] >= median);
        }
    }

  dataset_free (dataset);
}
