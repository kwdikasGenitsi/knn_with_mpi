#include "array.h"
#include <stdio.h>
#include <stdlib.h>

void
array_free (Array *array)
{
  free (array->data);
  free (array);
}

Array *
array_new (size_t size)
{
  Array *array = (Array *) malloc (sizeof (*array));
  array->data = (number_t *) malloc (sizeof (number_t) * size);
  array->size = size;
  return array;
}

void
array_dump (Array *array)
{
  printf ("[");
  for (size_t i = 0; i < array->size; i++)
    {
      printf ("%s%f", i > 0 ? ", " : "", array->data[i]);
    }
  printf ("]\n");
}