#pragma once

typedef float number_t;


typedef struct
{
  number_t *data;
  int size;
} Array;

Array* array_new(int size);
void array_free (Array *array);



typedef struct
{
	Array *buffer;
	int top;
} MasterBuffer;



MasterBuffer* master_buffer_new(int size);
void master_buffer_free(MasterBuffer *master_buffer);

int master_buffer_fill(MasterBuffer *master_buffer, Array *new_points);
Array* master_buffer_throw(MasterBuffer *master_buffer, int number_of_points);
