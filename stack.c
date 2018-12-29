#include "stack.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

Stack *
stack_new (size_t max_size)
{
  Stack *stack = malloc (sizeof (*stack));
  stack->stack = malloc (sizeof (*stack->stack) * max_size);
  stack->max_size = max_size;
  stack->top = 0;
  return stack;
}

void
stack_free (Stack *stack)
{
  free (stack->stack);
  free (stack);
}

void
stack_push (Stack *stack, number_t element)
{
  /* Prevent a stack overflow. */
  assert ((stack->top + 1) < stack->max_size);

  stack->stack[stack->top++] = element;
}

number_t
stack_pop (Stack *stack)
{
  /* Prevent a stack underflow. */
  assert (stack->top > 0);

  return stack->stack[--stack->top];
}

void
stack_push_array (Stack *stack, Array *array)
{
  /* Prevent a stack overflow. */
  assert ((stack->top + array->size) < stack->max_size);

  memcpy (stack->stack, array->data, sizeof (number_t) * array->size);
  stack->top += array->size;
}

void
stack_pop_array (Stack *stack, Array *array, size_t offset, size_t length)
{
  /* Make sure the array can store the result. */
  assert ((array->size - offset) >= length);

  /* Prevent an underflow. */
  assert (stack->top >= length);

  memcpy (array->data + offset, stack->stack + stack->top - length,
          length * sizeof (number_t));

  stack->top -= length;
}

size_t
stack_get_max_max_size (Stack *stack)
{
  return stack->max_size;
}

size_t
stack_get_size (Stack *stack)
{
  return stack->top;
}

/**
 * Gets a random item from the stack.
 * The index is checked against the array boundaries with an
 * assertion.
 */
number_t
stack_get_at (Stack *stack, size_t index)
{
  assert (index < stack->top);
  return stack->stack[index];
}

void
test_stack ()
{
  Stack *stack = stack_new (5);
  stack_push (stack, (number_t) 57.0f);
  stack_push (stack, (number_t) 11.0f);

  assert (stack_get_at (stack, 0) == (number_t) 57.0f);
  assert (stack_get_at (stack, 1) == (number_t) 11.0f);
  assert (stack_get_max_max_size (stack) == 5);
  assert (stack_get_size (stack) == 2);
  assert (stack_pop (stack) == (number_t) 11.0f);
  assert (stack_get_size (stack) == 1);
  assert (stack_pop (stack) == (number_t) 57.0f);
  assert (stack_get_size (stack) == 0);

  stack_free (stack);
}

size_t
stack_get_free_space (Stack *stack)
{
  return stack->max_size - stack->top;
}
