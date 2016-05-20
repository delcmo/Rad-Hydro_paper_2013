import numpy as np

def b(xq,p):
  p+=1
  dx = np.linspace(-1, 1, num=p)
  v = [1]*p
  """Calculate b_j(x_xi)"""
  for i in xrange(p):
    for k in xrange(p):
      if k != i:
        v[i] *= (xq-dx[k]) / (dx[i]-dx[k])
  return v