#include "array.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

Array
array_new (size_t size)
{
  Array array;
  array.data = (number_t *) malloc (sizeof (number_t) * size);
  array.size = size;
  return array;
}

Array
array_copy (Array rhs)
{
  Array new_array;
  new_array.data = (number_t *) malloc (sizeof (number_t) * rhs.size);
  new_array.size = rhs.size;
  memcpy (new_array.data, rhs.data, sizeof (number_t) * rhs.size);
  return new_array;
}

void
array_free (Array array)
{
  free (array.data);
}

ArraySlice
array_get_slice (Array array, size_t offset, size_t length)
{
  ArraySlice slice;
  slice.data = array.data + offset;
  slice.size = length - offset;
  return slice;
}

void
array_slice_free (ArraySlice *slice)
{
  free (slice);
}

void
array_dump (Array array)
{
  printf ("[");
  for (size_t i = 0; i < array.size; i++)
    {
      printf ("%s" PRINTF_NUMBER_FMT, i > 0 ? ", " : "", array.data[i]);
    }
  printf ("]\n");
}

void
array_slice_dump (ArraySlice slice)
{
  printf ("[");
  for (size_t i = 0; i < slice.size; i++)
    {
      printf ("%s" PRINTF_NUMBER_FMT, i > 0 ? ", " : "", slice.data[i]);
    }
  printf ("]\n");
}

void
array_fill_random (Array array)
{
  int cal = 5;
  MPI_Comm_rank (MPI_COMM_WORLD, &cal);
  srand ((cal + 1) * time (NULL));
  for (size_t i = 0; i < array.size; i++)
    array.data[i] = (number_t) (rand () % 100000 - rand () % 100000) * 0.0001f;
}

size_t
heap_parent_of (size_t i)
{
  return (i - 1) / 2;
}

size_t
heap_left_child_of (size_t i)
{
  return 2 * i + 1;
}

size_t
heap_right_child_of (size_t i)
{
  return 2 * i + 2;
}