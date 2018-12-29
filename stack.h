#pragma once
#include "array.h"
#include <stddef.h>

typedef struct
{
  number_t *stack;
  size_t max_size;
  size_t top;
} Stack;

Stack *stack_new (size_t max_size);
void stack_free (Stack *vp_stack);

void stack_push (Stack *vp_stack, number_t element);
number_t stack_pop (Stack *vp_stack);

/**
 * Pushes every element from an Array to the stack.
 */
void stack_push_array (Stack *stack, Array *array);

/**
 * Pops a number of elements from the stack and writes them, in the same
 * order as they were pushed, in an array.
 * @param stack The stack to pop from.
 * @param array The array to write the data to.
 * @param offset The offset within the array to start the write at.
 * @param length The number of elements to be popped.
 */
void stack_pop_array (Stack *stack, Array *array, size_t offset, size_t length);

size_t stack_get_max_size (Stack *vp_stack);
size_t stack_get_size (Stack *vp_stack);

/**
 * Gets a random item from the stack.
 * The index is checked against the array boundaries with an
 * assertion.
 */
number_t stack_get_at (Stack *vp_stack, size_t index);

/**
 * Minimal unit test.
 */
void test_stack ();

size_t stack_get_free_space (Stack *stack);