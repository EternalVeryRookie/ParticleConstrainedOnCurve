import numpy as np


def tdma(left_mat, right_vector):
  n = left_mat.shape[0]
  P = [None] * n
  Q = [None] * n
  X = [None] * n
  P[0] = -left_mat[0][1] / left_mat[0][0]
  Q[0] = right_vector[0] / left_mat[0][0]
  for i in range(1, n-1):
    P[i] = -left_mat[i][i + 1] / (left_mat[i][i] + left_mat[i][i - 1]*P[i - 1])

  for i in range(1, n):
    Q[i] = (right_vector[i] - left_mat[i][i - 1] * Q[i - 1]) / (left_mat[i][i] + left_mat[i][i - 1] * P[i - 1])

  X[n - 1] = Q[n - 1]
  for i in range(n-2, -1, -1):
    X[i] = P[i] * X[i + 1] + Q[i]

  return X
