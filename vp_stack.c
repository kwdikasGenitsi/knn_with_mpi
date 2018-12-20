#include "vp_stack.h"
#include <assert.h>
#include <stdlib.h>

VpStack *
vp_stack_new (size_t length)
{
  VpStack *vp_stack = malloc (sizeof (*vp_stack));
  vp_stack->stack = malloc (sizeof (*vp_stack->stack) * length);
  vp_stack->length = length;
  vp_stack->top = 0;
  return vp_stack;
}

void
vp_stack_free (VpStack *vp_stack)
{
  free (vp_stack->stack);
  free (vp_stack);
}

void
vp_stack_push (VpStack *vp_stack, VpNode node)
{
  /* Make sure the item fits. */
  assert ((vp_stack->top + 1) < vp_stack->length);
  vp_stack->stack[vp_stack->top++] = node;
}

VpNode
vp_stack_pop (VpStack *vp_stack)
{
  assert (vp_stack->top > 0);
  return vp_stack->stack[--vp_stack->top];
}

size_t
vp_stack_get_max_length (VpStack *vp_stack)
{
  return vp_stack->length;
}

size_t
vp_stack_get_size (VpStack *vp_stack)
{
  return vp_stack->top;
}

/**
 * Gets a random item from the stack.
 * The index is checked against the array boundaries with an
 * assertion.
 */
VpNode
vp_stack_get_at (VpStack *vp_stack, size_t index)
{
  assert (index < vp_stack->top);
  return vp_stack->stack[index];
}

void
test_vp_stack ()
{
  VpStack *stack = vp_stack_new (5);
  vp_stack_push (stack, (VpNode){.vantage_point = 3.0, .median = 1.0});
  vp_stack_push (stack, (VpNode){.vantage_point = 5.0, .median = 2.0});

  assert (vp_stack_get_at (stack, 0).median == 1.0);
  assert (vp_stack_get_at (stack, 1).vantage_point == 5.0);
  assert (vp_stack_get_max_length (stack) == 5);
  assert (vp_stack_get_size (stack) == 2);
  assert (vp_stack_pop (stack).median == 2.0);
  assert (vp_stack_get_size (stack) == 1);
  assert (vp_stack_pop (stack).median == 1.0);
  assert (vp_stack_get_size (stack) == 0);

  vp_stack_free (stack);
}
