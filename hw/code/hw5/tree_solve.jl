# [x] = tree_solve(p, d, w, b)
#
# Solve Ax = b where A is a symmetric positive definite matrix
# with diagonal d and off-diagonal entries A(i,p(i)) = A(p(i),i) = w(i)
# (except where p(i) = 0).  Your solution should run in O(n) time, and
# should *not* use the sparse solvers in MATLAB or Julia directly --
# write your own tree factorization!  You may assume that p is in ascending
# order (i.e. the nodes are ordered from leaves to root).
#
using LinearAlgebra
function tree_solve(p, d, w, b)
  A = zeros(Float64, length(p), length(p))
  L = zeros(Float64, length(p), length(p))
  x = zeros(Float64, length(p))
  y = zeros(Float64, length(p))
  for i = 1:length(p)
    A[i, i] = d[i]
    if p[i] > 0
      A[i, p[i]] = w[i]
      A[p[i], i] = w[i]
    end
  end

  L = A
  for i = 1:length(p)
    L[i, i] = sqrt(L[i,i])
    # L[i+1:end,i] /= L[i,i]
    # L[i+1:end,i+1:end] -= L[i+1:end,i]*L[i+1:end,i]'
    if p[i] > 0
      L[p[i], i] = L[p[i], i]/L[i,i]
      L[p[i], p[i]] -= L[p[i], i]^2
    end
  end
  L = tril(L)
  for i = 1:length(p)
    y[i] = b[i]/L[i, i]
    if p[i] > 0
      b[p[i]] -= y[i]*L[p[i],i]
    end
  end
  for i = length(p):-1:1
    if p[i] > 0
      y[i] -= x[p[i]]*L[p[i],i]
    end
    x[i] = y[i]/L[i,i]
  end

  return x
end

p = [2,3,0]
d = [2,2,2]
w = [-1,-1,0]
b = rand(3)
x = tree_solve(p,d,w,b)
