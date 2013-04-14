#include <VDM.h>

int choose_pivot(bigrational_matrix& X, const unsigned long diag_num)
{
  assert(X.get_num_row() == X.get_num_col());
  assert(X.get_num_row() > diag_num);

  unsigned long i, j, k;
  unsigned long ord = X.get_num_row();
  int           c;
  bigrational   y;

  //  1. Find a pivot X(i, j).

  c = 1;
  i = j = diag_num;

  while (c && i < ord)
  {
    while (c && j < ord)
      if (X(i, j) != 0)
        c = 0;
      else
        j++;

    if (c)
    {
      i++;
      j = diag_num;
    }
  }

  assert(!c);
  assert(i < ord && j < ord);

  //  2. Interchange column diag_num and column j.

  if (j != diag_num)
    for (k = 0; k < ord; k++)
    {
      y              = X(k, diag_num);
      X(k, diag_num) = X(k, j);
      X(k, j)        = y;
    }

  //  3. Interchange row diag_num and row i.

  if (i != diag_num)
    for (k = diag_num; k < ord; k++)
    {
      y              = X(diag_num, k);
      X(diag_num, k) = X(i, k);
      X(i, k)        = y;
    }

  return 0;
}

bigrational_vector Solve_VDM_GE(const bigrational_vector& Val,
                                const bigrational_vector& RHS)
{
  assert(Val.get_dim() == RHS.get_dim());

  long               i, j;
  unsigned long      k;
  unsigned long      ord = Val.get_dim();
  bigrational_matrix VDM(ord, ord);
  bigrational        f;
  bigrational_vector Sol(ord);

  //  1. Setup Vandermonte matrix VDM and initialize the solution vector Sol.

  for (i = 0; i < ord; i++)
  {
    VDM(i, 0) = 1;
    Sol[i]    = RHS[i];
  }

  for (i = 0; i < ord; i++)
    for (j = 1; j < ord; j++)
      VDM(i, j) = VDM(i, j - 1) * Val[i];

  //  2. Apply Gaussian elimination.

  //  2-1. Diagonalize the system.

  for (k = 0; k < ord - 1; k++)
  {
    choose_pivot(VDM, k);

    for (i = k + 1; i < ord; i++)
    {
      f         = VDM(i, k) / VDM(k, k);
      VDM(i, k) = 0;

      for (j = k + 1; j < ord; j++)
        VDM(i, j) -= f * VDM(k, j);

      Sol[i] -= f * Sol[k];
    }
  }

  assert(VDM(ord - 1, ord - 1) != 0);  //  Assert non-singularity.

  //  2-2. Solve the diagonal system.

  for (i = ord - 1; i >= 0; i--)
  {
    for (j = ord - 1; j > i; j--)
      Sol[i] -= Sol[j] * VDM(i, j);

    Sol[i] /= VDM(i, i);
  }

  return Sol;
}

