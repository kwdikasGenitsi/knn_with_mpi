#pragma once

#include <assert.h>
#include <stddef.h>

#define UNUSED(expr)                                                           \
  do                                                                           \
    {                                                                          \
      (void) (expr);                                                           \
    }                                                                          \
  while (0)

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

typedef struct
{
  number_t *data;
  size_t size;
} ArraySlice;

/**
 * Creates a new array.
 * @param size The number of elements the array must hold.
 * @return A pointer to an Array structure that must be freed
 *         with array_free().
 */
Array array_new (size_t size);

/**
 * Creates a new array that points to a subarray within the original Array.
 * We call this a "slice", in accordance with Rust terminology.
 * Slices are allocated statically and need not be freed.
 */
ArraySlice array_get_slice (Array array, size_t offset, size_t length);

/**
 * Frees the memory occupied by an Array.
 */
void array_free (Array array);

/**
 * Prints the contents of the array in human-readable form to stdout.
 * For instance "[1, 2, 3]\n"
 * @param array The array to dump.
 */
void array_dump (Array array);

void array_slice_dump (ArraySlice slice);

/**
 * Fills the array with random values.
 */
void array_fill_random (Array array);

/* Heap indexing helpers. */
size_t heap_parent_of (size_t i);
size_t heap_left_child_of (size_t i);
size_t heap_right_child_of (size_t i);

/**
 * To know the depth of our tree, we need to calculate log2(x) accurately.
 * We can't use floating point numbers since we risk getting an off-by-one error
 * due to floating point accuracy.
 * We therefore use this GCC bultin, which returns the number of leading 0s in
 * an integer. On x86 this is done with the BSR instruction.
 */
static inline int
log2i (int x)
{
  assert (x > 0);
  return sizeof (int) * 8 - __builtin_clz (x) - 1;
}
