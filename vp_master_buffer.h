#pragma once

#include "array.h"

#if 0
typedef struct
{
  Array *buffer;
  int top;
} MasterBuffer;

int free_space (MasterBuffer *master_buffer);
int filled_space (MasterBuffer *master_buffer);

MasterBuffer *master_buffer_new (int size);
void master_buffer_free (MasterBuffer *master_buffer);

int master_buffer_fill (MasterBuffer *master_buffer, Array *new_points);
Array *master_buffer_throw (MasterBuffer *master_buffer, int number_of_points);
#endif