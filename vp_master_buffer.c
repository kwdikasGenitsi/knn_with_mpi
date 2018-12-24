#include "vp_master_buffer.h"
#include <stdio.h>
#include <stdlib.h>

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

int
free_space (MasterBuffer *master_buffer)
{
  return master_buffer->buffer->size - master_buffer->top;
}

int
filled_space (MasterBuffer *master_buffer)
{
  return master_buffer->top;
}

// adds the points of an array to the master_buffer if there is enough space.
// returns 0 if master_buffer has been filled successfully or -1 if there was no
// space.
int
master_buffer_fill (MasterBuffer *master_buffer, Array *new_points)
{
  if (free_space (master_buffer) >= new_points->size)
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

  if (filled_space (master_buffer) >= number_of_points)
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
