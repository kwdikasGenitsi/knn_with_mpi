#pragma once
#include <stddef.h>

typedef float number_t;

typedef struct
{
  number_t vantage_point;
  number_t median;
} VpNode;

typedef struct
{
  VpNode *stack;
  size_t length;
  size_t top;
} VpStack;

VpStack *vp_stack_new (size_t length);
void vp_stack_free (VpStack *vp_stack);

void vp_stack_push (VpStack *vp_stack, VpNode node);
VpNode vp_stack_pop (VpStack *vp_stack);

size_t vp_stack_get_max_length (VpStack *vp_stack);
size_t vp_stack_get_size (VpStack *vp_stack);

/**
 * Gets a random item from the stack.
 * The index is checked against the array boundaries with an
 * assertion.
 */
VpNode vp_stack_get_at (VpStack *vp_stack, size_t index);

/**
 * Minimal unit test.
 */
void test_vp_stack ();
