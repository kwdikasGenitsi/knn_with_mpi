#include "vp_master_buffer.h"
#include <stdlib.h>
#include <stdio.h>

void
array_free (Array *array)
{
  free (array->data);
  free (array);
}

Array *
array_new (int size)
{
  Array *array = (Array *) malloc (sizeof (*array));
  array->data = (number_t *) malloc (sizeof (number_t) * size);
  array->size = size;
  return array;
}

MasterBuffer *
master_buffer_new (int size)
{
  MasterBuffer *master_buffer
    = (MasterBuffer *) malloc (sizeof (*master_buffer));
  master_buffer->buffer = array_new (size);
  master_buffer->top = 0;
  return master_buffer;
}
void
master_buffer_free (MasterBuffer *master_buffer)
{
  array_free (master_buffer->buffer);
  free (master_buffer);
}

// adds the points of an array to the master_buffer if there is enough space.
// returns 0 if master_buffer has been filled successfully or -1 if there was no
// space.
int
master_buffer_fill (MasterBuffer *master_buffer, Array *new_points)
{
  int free_space = master_buffer->buffer->size - master_buffer->top;
  if (free_space >= new_points->size)
    {
      int k = 0;
      for (int i = master_buffer->top;
           i < master_buffer->top + new_points->size; i++)
        {
          master_buffer->buffer->data[i] = new_points->data[k];
          k++;
        }
      master_buffer->top += new_points->size;
      return 0;
    }
  else
    {
      return -1;
    }
}

// pops an array of size number_of_points with points from the master_buffer.
Array *
master_buffer_throw (MasterBuffer *master_buffer, int number_of_points)
{

  if (master_buffer->top >= number_of_points)
    {
      Array *output_array = array_new (number_of_points);
      int k = master_buffer->top - 1;
      for (int i = 0; i < number_of_points; i++)
        {
          output_array->data[i] = master_buffer->buffer->data[k];
          k--;
        }
      master_buffer->top -= number_of_points;
      return output_array;
    }
  else
    {
      printf ("master_buffer_throw: not enough points to throw!");
      return NULL;
    }
}
