#pragma once

#include "Const.h"

extern double tempera[l + 3][m + 3]; //���x
extern double rho[l + 3][m + 3]; //�R���̖��x

void solvePoisson();
void solvePoissonByCG();