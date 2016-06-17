#### import python modules ####
import sys, os, math, csv, scipy
import numpy as np
import itertools as it
from Lagrange_test_function import b
from scipy.interpolate import interp1d
from scipy.interpolate import splrep, splev
import matplotlib.pyplot as plt
#### import python modules ####

def compute_error_norms(x_offset, x_coord, num_value, x_coord_exact, exact_value, quad_order):
  # SET VARIABLES
  x_coord_exact_offset = [float(x)+float(x_offset) for x in x_coord_exact]
  num_value = [float(x) for x in num_value]
  x_coord = [float(x) for x in x_coord]
  exact_value = [float(x) for x in exact_value]
  nb_nodes = len(x_coord)
  nb_cells = nb_nodes-1
  nb_nodes_exact = len(x_coord_exact_offset)
  l2_norm_cell = [0]*nb_cells
  l1_norm_cell = [0]*nb_cells
  # compute quadrature points and weights
  [xq, wq] = np.polynomial.legendre.leggauss(quad_order)
  # compute values of Lagrange functions at quadrature points
  lq = b(xq,1)
  lq = np.asarray(lq)
  jac = 0.5*(x_coord[-1]-x_coord[0])/nb_cells
  # check
  if (nb_nodes_exact < nb_nodes):
    print 'Number of exact nodes is lower than number of numerical nodes.'
    exit()

  # LOOP OVER NUMERICAL MESH
  index_left = 0
  for cell in xrange(0,nb_cells,1):
    # get coordinates of left and right nodes
    node_left = x_coord[cell]
    node_right = x_coord[cell+1]
    diff_node = node_right - node_left
#    jac = 0.5*diff_node
    sum_node = node_right + node_left
    # compute coordinates of xq in cell 'cell'
    xq_cell = 0.5*(xq*diff_node+sum_node)
    # get the indeces of closest nodes to the cell from the exact mesh
    while (x_coord_exact_offset[index_left] < node_left):
      index_left += 1
    index_left = index_left-10 if index_left !=0 else index_left
    index_right = index_left
    while (x_coord_exact_offset[index_right] < node_right):
      index_right += 1
    index_right = index_right+10 if index_right != nb_nodes_exact-1  else nb_nodes_exact
    # reduce the number of nodes from analytical solution down to 100
    interval=1
    if index_right-index_left>100:
      interval=int((index_right-index_left)/100)
    exact_value_cell = exact_value[index_left:index_right:interval]
    exact_value_cell[0] = exact_value[index_left]
    exact_value_cell[-1] = exact_value[index_right]
    x_coord_exact_cell = x_coord_exact_offset[index_left:index_right:interval]
    x_coord_exact_cell[0] = x_coord_exact_offset[index_left]
    x_coord_exact_cell[-1] = x_coord_exact_offset[index_right]
    # interpolate values at the quadrature points 'xq_cell'
#    f = splrep(x_coord_exact_cell, exact_value_cell, s=0)
#    exact_value_cell_xq = splev(xq_cell, f, der=0)
    f = interp1d(x_coord_exact_cell, exact_value_cell, kind='cubic')
    exact_value_cell_xq = f(xq_cell)
    # compute the numerical solution value at the quadrature points
    num_value_cell = num_value[cell:cell+2]
    value_xq = [0]*quad_order
    for qp in xrange(quad_order):
      value_xq[qp] = np.dot(lq[:,qp],num_value_cell)
    # compute mass difference between the exact and numerical solutions
    l1_norm_cell[cell] = np.dot(abs(value_xq-exact_value_cell_xq),wq)*jac
    l2_norm_cell[cell] = np.dot((value_xq-exact_value_cell_xq)**2, wq)*jac

  l2_norm = math.sqrt(sum(l2_norm_cell))
  l1_norm = sum(l1_norm_cell)
  return (l1_norm, l2_norm)

