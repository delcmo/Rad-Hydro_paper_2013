% clear all; close all; clc;
function [normL1, normL2, n_cells] = post_process_norm(numeric_filename,options,...
                                              exact_filename, x_offset)
    
% numeric_filename = file storing the numerical solution
% options = quadrature points, interpolation type, index of the numerical
% solution
% exact_filename = used to pass the text files storing the exact solution (x,rho,eps,v,temp)

% name of the variables
variable_names = {'density'; 'radiation'; 'mach'; 'material temperature'};

% load data from exact and numeric filenames
nb_exact_file = size(exact_filename);
for ifile=1:nb_exact_file(2)
   file_id = fopen(exact_filename{ifile});
   exact_value(:,ifile) = textread(exact_filename{ifile}, '%f');
   fclose(file_id);   
end

% make sure the x-values from the exact solution are linearly increasing 
% (remove values with the same x-coord)
[exact_value_unique(:,1), index_unique] = unique(exact_value(:,1));
for ifile=2:nb_exact_file(2)
    exact_value_unique(:,ifile) = exact_value(index_unique,ifile);
end
exact_value_unique(:,1) = exact_value_unique(:,1) - x_offset;
nb_exact = size(exact_value_unique(:,1));
nb_exact = nb_exact(1);

%index for numerical solution
index_num = options.index;
nb_index_num = size(options.index);
assert(nb_index_num(1)==nb_exact_file(2),'the number of exact files does not match the number of indexes of the numerical solutions.');

% load vectors from 'num_filename'
numeric = csvread(numeric_filename,1,0);
for i=1:nb_index_num(1)
    numeric_value(:,i) = numeric(:,index_num(i));
end

% normalize the material temperature and radiation variables
% for both exact and numerical values.
for i=2:nb_index_num(1)
    exact_value_unique(:,i) = exact_value_unique(:,i) / exact_value_unique(1,i);
    numeric_value(:,i) = numeric_value(:,i) / numeric_value(1,i);    
end

% get number of cells and return the value
n_cells = size(numeric_value(:,1))-1;
n_cells = n_cells(1);

% load 1d quadrature
[xq,wq] = GLNodeWt(options.nquad);

% jacobian term
jac = 0.5 * (numeric_value(end,1)-numeric_value(1,1)) / n_cells;

% load test functions
f = Lagrange_poly(xq,1);

% initialize norms
normL1 = zeros(nb_index_num(1)-1,1);
normL2_sq = zeros(nb_index_num(1)-1,1);

% initialize integer
index_left = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over the cells of the numerical mesh %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_cell=1:1:n_cells
% coordinates of the nodes belonging to cell 'i_cell'
node_cell_coord = numeric_value(i_cell:i_cell+1,1);
diff_node = node_cell_coord(2)-node_cell_coord(1);
sum_node = node_cell_coord(1)+node_cell_coord(2);

% coordinates of the quadrature points in the real mesh
xq_coord = xq*0.5*diff_node+0.5*sum_node;

% find node of the semi-analytical mesh that is the closest from the left
% to 'xq_coord(1)'
while (exact_value_unique(index_left,1)<node_cell_coord(1))
    index_left = index_left+1;
end
index_right = index_left;
while (exact_value_unique(index_right,1)<node_cell_coord(end))
    index_right = index_right+1;
end
if (i_cell ~= 1) 
    index_left = index_left-1; 
end
if(i_cell ~= n_cells)
    index_right = index_right+1;
end
if (i_cell == 1), index_first_node_exact_sol = index_left; end

% % % % % find node of the semi-analytical mesh that is the closest from the left
% % % % % to 'xq_coord(1)'
% % % % while (exact_value_unique(index_left,1)<node_cell_coord(1))
% % % %     index_left = index_left+1;
% % % %     if (index_left>nb_exact), break; end
% % % % end
% % % % index_left = min(index_left, nb_exact);
% % % % index_right = index_left;
% % % % while (exact_value_unique(index_right,1)<node_cell_coord(end))
% % % %     index_right = index_right+1;
% % % %     if (index_right>nb_exact), break; end
% % % % end
% % % % if (i_cell ~= 1), index_left = max(index_left-1,1); end
% % % % if (i_cell == 1), index_first_node_exact_sol = index_left; end
% % % % if(i_cell ~= n_cells), index_right = index_right+1; end
% % % % if (i_cell == n_cells), index_right = nb_exact; index_left = index_right-options.nquad; end
% % % % index_right = min(index_right, nb_exact);

% % % % normalize the exact solution using the node closest to first node of
% % % % numerical mesh
% % % % for i=2:nb_index_num(1)
% % % %     exact_value(:,i) = exact_value(:,i) / exact_value(index_first_node_exact_sol,i);
% % % % end

% interpolate the values of the exact solution at the quadrature points
% (only for rho, eps, mach and temp)
exact_xq = zeros(options.nquad,nb_exact_file(2)-1);
for ifile=2:nb_exact_file(2)
    exact_xq(:,ifile-1) = interp1(exact_value_unique(index_left:index_right,1), exact_value_unique(index_left:index_right,ifile), xq_coord, options.interpolation_type);
end

% interpolate the values of the numerical solution at the quadrature points
numeric_xq = zeros(options.nquad,nb_exact_file(2)-1);
for ifile=2:nb_index_num(1)
    for nq=1:options.nquad
        numeric_xq(nq,ifile-1) = dot(numeric_value(i_cell:i_cell+1,ifile),f(nq,:));
    end
end

% compute the cell contribution to the L1 and L2 norms
for i=1:nb_index_num(1)-1
    normL1(i) = normL1(i) + dot(wq,abs(exact_xq(:,i)-numeric_xq(:,i)))*jac;
    normL2_sq(i) = normL2_sq(i) + dot(wq,(exact_xq(:,i)-numeric_xq(:,i)).^2)*jac;
end

% set index_left equal to index_right-1 to speed up next search
index_left = index_right-1;
end % end of loop over cells

% compute L2 norm
normL2 = sqrt(normL2_sq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print the values of the norms %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (options.output)
    for i=1:nb_index_num(1)-1
        fprintf('L1 norm in domain for the %s is: %12.7e \n', variable_names{i}, normL1(i));
        fprintf('L2 norm in domain for the %s is %12.7e \n \n', variable_names{i}, normL2(i));    
        fprintf('----------------------------------------------------------- \n');    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot data if required %
%%%%%%%%%%%%%%%%%%%%%%%%%
if (options.plot)
   % xmin and xmax
   xmin = min(numeric_value(:,1));
   xmax = max(numeric_value(:,1));
   % legend 
   lgd{1} = 'exact';
   lgd{2} = 'numerical';
   % title
   ttl{1}='density'; ttl{2}='radiation'; 
   ttl{3}='mach number'; ttl{4}='material temperature';
   for i=2:nb_index_num(1)
       figure(n_cells+i)
       plot(exact_value_unique(:,1),exact_value_unique(:,i),'-',numeric_value(:,1),numeric_value(:,i),'-');
       legend(lgd);
       title(ttl{i-1});
       xlim([xmin xmax]);
   end
end
end % end function