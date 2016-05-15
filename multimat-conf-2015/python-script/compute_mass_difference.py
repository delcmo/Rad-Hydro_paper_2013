#### import python modules ####
import sys, os, math, csv, scipy
import numpy as np
import itertools as it
from Lagrange_test_function import b
from scipy.interpolate import interp1d
from scipy.interpolate import splrep, splev
import matplotlib.pyplot as plt
from decimal import *
#### import python modules ####
np.set_printoptions(precision=20)

def compute_mass_diff(x_offset, x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order):
  # SET VARIABLES
  x_coord_exact_offset = [float(x)+float(x_offset) for x in x_coord_exact]
  mat_density = [float(x) for x in mat_density]
#  mat_density = [float(x)/float(mat_density[0]) for x in mat_density]
  x_coord = [float(x) for x in x_coord]
  mat_density_exact = [float(x)/float(mat_density_exact[0]) for x in mat_density_exact]
#  mat_density_exact = [float(x)/float(mat_density_exact[0]) for x in mat_density_exact]
  nb_nodes = len(x_coord)
  nb_cells = nb_nodes-1
  nb_nodes_exact = len(x_coord_exact_offset)
  mass_num = [0]*nb_cells
  mass_exact = [0]*nb_cells
  mass_diff = 0
  # compute quadrature points and weights
  [xq, wq] = np.polynomial.legendre.leggauss(quad_order)
  xq = np.asarray(xq)
  # compute values of Lagrange functions at quadrature points
  lq = b(xq,1)
  lq = np.array(lq, dtype=float)
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
    jac = 0.5*diff_node
    sum_node = node_right + node_left
    # compute coordinates of xq in cell 'cell'
    xq_cell = 0.5*(xq*diff_node+sum_node)
    # get the indeces of closest nodes to the cell from the exact mesh
    while (x_coord_exact_offset[index_left] < node_left):
      index_left += 1
    index_left = index_left-2 if index_left !=0 else index_left
    index_right = index_left
    while (x_coord_exact_offset[index_right] < node_right):
      index_right += 1
    index_right = index_right+2 if index_right != nb_nodes_exact-1  else nb_nodes_exact
    # interpolate values at the quadrature points 'xq_cell'
    exact_value_cell = mat_density_exact[index_left:index_right]
    x_coord_exact_cell = x_coord_exact_offset[index_left:index_right]
#    f = splrep(x_coord_exact_cell, exact_value_cell, s=0)
#    exact_value_cell_xq = splev(xq_cell, f, der=0)
    f = interp1d(x_coord_exact_cell, exact_value_cell, kind=3)
    exact_value_cell_xq = f(xq_cell)
    # compute the numerical solution value at the quadrature points
    mat_density_cell = mat_density[cell:cell+2]
    value_xq = [0]*quad_order
    for qp in xrange(quad_order):
      value_xq[qp] = np.dot(lq[:,qp],mat_density_cell)
    # compute mass difference between the exact and numerical solutions
    mass_num[cell] = np.dot(value_xq,wq)*jac
    mass_exact[cell] = np.dot(exact_value_cell_xq, wq)*jac
    mass_diff += np.dot(exact_value_cell_xq-value_xq, wq)*jac

#  plt.plot(x_coord_exact_offset, mat_density_exact)
#  plt.plot(x_coord, mat_density)
#  plt.show()

  return abs(mass_diff)