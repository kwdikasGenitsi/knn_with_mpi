#pragma once

#include <stddef.h>

/**
 * The number type used to represent the dataset.
 * Defined here so that it can be changed to a double.
 */
typedef float number_t;
#define MPI_NUMBER_T MPI_FLOAT
#define PRINTF_NUMBER_FMT "%f"

typedef struct
{
  number_t *data; /**< Pointer to a data buffer of size size. */
  size_t size;    /**< The number of elements the array can store. */
} Array;

/**
 * Creates a new array.
 * @param size The number of elements the array must hold.
 * @return A pointer to an Array structure that must be freed
 *         with array_free().
 */
Array *array_new (size_t size);

/**
 * Frees the memory occupied by an Array.
 */
void array_free (Array *array);

/**
 * Prints the contents of the array in human-readable form to stdout.
 * For instance "[1, 2, 3]\n"
 * @param array The array to dump.
 */
void array_dump (Array *array);

/**
 * Fills the array with random values.
 */
void array_fill_random (Array *array);