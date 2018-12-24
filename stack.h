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
