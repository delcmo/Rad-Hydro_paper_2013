function [ diff_rhoE ] = compute_shift_with_total_energy_conservation(numeric_filename,options,...
                                              exact_filename, x_offset)
% numeric_filename = file storing the numerical solution
% options = quadrature points, interpolation type, index of the numerical
% solution
% exact_filename = used to pass the text files storing the exact solution (x,rho,eps,Mach or v,temp)

% NOTE: in this function, ALL of the semi-analytical solutions are required
% since the total energy (material + radiation) is computed.

% load data from exact and numeric filenames
nb_exact_file = size(exact_filename);
for ifile=1:nb_exact_file(2)
   file_id = fopen(exact_filename{ifile});
   exact_value(:,ifile) = textread(exact_filename{ifile}, '%f');
   fclose(file_id);   
end
nb_exact = size(exact_value(:,1));
nb_exact = nb_exact(1);

% make sure the x-values from the exact solution are linearly increasing 
% (remove values with the same x-coord)
[exact_value_unique(:,1), index_unique] = unique(exact_value(:,1));
for ifile=2:nb_exact_file(2)
    exact_value_unique(:,ifile) = exact_value(index_unique,ifile);
end

%index for numerical solution
index_num = options.index;
nb_index_num = size(options.index);
% assert(nb_index_num(1)==nb_exact_file(2),'the number of exact files does not match the number of indexes for the numerical solutions.');

% load vectors from 'num_filename'
numeric = csvread(numeric_filename,1,0);
for i=1:nb_index_num(1)
    numeric_value(:,i) = numeric(:,index_num(i));
end

% normalize the variables for both exact and numerical values. 
% for i=2:nb_index_num(1)
%     exact_value(:,i) = exact_value(:,i) / exact_value(1,i);
%     numeric_value(:,i) = numeric_value(:,i) / numeric_value(1,i);
% end
% exact_value_unique(:,3) = exact_value_unique(:,3)/exact_value_unique(1,3);
% numeric_value(:,3) = numeric_value(:,3)/numeric_value(1,3);

% offset the x-vector from the semi-analytical solution
exact_value_unique(:,1) = exact_value_unique(:,1) - x_offset;

% compute the numerical and semi-analytical total energy vectors:
gamma=options.eos(1);
Cv=options.eos(2);
a=options.eos(3);

rhoE_star_num = 0.;
% rhoE_star_num = Cv*numeric_value(:,2).*numeric_value(:,4);
% rhoE_star_num = rhoE_star_num/rhoE_star_num(1,1);
% rhoE_star_num_bis = 0.5*numeric_value(:,3).*numeric_value(:,3)./numeric_value(:,2);
% rhoE_star_num_bis = rhoE_star_num_bis/rhoE_star_num_bis(1,1);
% rhoE_star_num = rhoE_star_num + rhoE_star_num_bis;
% rhoE_star_num = rhoE_star_num/rhoE_star_num(1,1);
% rhoE_star_num_bis = numeric_value(:,5);
% rhoE_star_num_bis = rhoE_star_num_bis/rhoE_star_num_bis(1,1);
% rhoE_star_num = rhoE_star_num + rhoE_star_num_bis;
% rhoE_star_num = rhoE_star_num/rhoE_star_num(1,1);

rhoE_star_num = numeric_value(:,2);
rhoE_star_num = rhoE_star_num/rhoE_star_num(1,1);
rhoE_star_num_bis = numeric_value(:,3);
rhoE_star_num_bis = rhoE_star_num_bis/rhoE_star_num_bis(1,1);
rhoE_star_num = rhoE_star_num + rhoE_star_num_bis;
rhoE_star_num = rhoE_star_num/rhoE_star_num(1,1);

rhoE_star_exact = 0.;
rhoE_star_exact = Cv*exact_value_unique(:,2).*exact_value_unique(:,5);
rhoE_star_exact = rhoE_star_exact/rhoE_star_exact(1,1);
rhoE_star_exact_bis = 0.5*exact_value_unique(:,2).*exact_value_unique(:,4).*exact_value_unique(:,4);
rhoE_star_exact_bis = rhoE_star_exact_bis/rhoE_star_exact_bis(1,1);
rhoE_star_exact = rhoE_star_exact + rhoE_star_exact_bis;
rhoE_star_exact = rhoE_star_exact/rhoE_star_exact(1,1);
rhoE_star_exact_bis = exact_value_unique(:,3);
rhoE_star_exact_bis = rhoE_star_exact_bis/rhoE_star_exact_bis(1,1);
rhoE_star_exact = rhoE_star_exact + rhoE_star_exact_bis;
rhoE_star_exact = rhoE_star_exact/rhoE_star_exact(1,1);

plot(numeric_value(:,1), rhoE_star_num, exact_value_unique(:,1), rhoE_star_exact);

% % % % % sp_num = Cv*gamma*(gamma-1)*numeric_value(:,5);
% % % % rhou_num = numeric_value(:,2); % .*numeric_value(:,4); % .*sp_num;
% % % % rhou_num = rhou_num/rhou_num(1,1);

% % % % % sp_exact = Cv*gamma*(gamma-1)*exact_value_unique(:,5);
% % % % rhou_exact = exact_value_unique(:,2).*exact_value_unique(:,4); % .*sp_exact;
% % % % rhou_exact = rhou_exact/rhou_exact(1,1);

% % % % plot(numeric_value(:,1), rhou_num, exact_value_unique(:,1), rhou_exact);

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
diff_rhoE = 0.;

% initialize integer
index_left = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over the cells of the numerical mesh %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_cell=1:1:n_cells
% coordinates of the nodes belonging to cell 'i_cell'
node_cell_coord = rhoE_star_num(i_cell:i_cell+1,1);
diff_node = node_cell_coord(2)-node_cell_coord(1);
sum_node = node_cell_coord(1)+node_cell_coord(2);

% coordinates of the quadrature points in the real mesh
xq_coord = xq*0.5*diff_node+0.5*sum_node;

% find node of the semi-analytical mesh that is the closest from the left
% to 'xq_coord(1)'

while (rhoE_star_exact(index_left,1)<node_cell_coord(1))
    index_left = index_left+1;
end
index_right = index_left;
while (rhoE_star_exact(index_right,1)<node_cell_coord(end))
    index_right = index_right+1;
end
if (i_cell ~= 1) 
    index_left = index_left-1; 
end
if(i_cell ~= n_cells)
    index_right = index_right+1;
end

% interpolate the values of the exact solution at the quadrature points
% (only for rho, eps, mach and temp)
exact_xq = zeros(options.nquad,1);
exact_xq(:,1) = interp1(exact_value_unique(index_left:index_right,1), rhoE_star_exact(index_left:index_right,1), xq_coord, options.interpolation_type);


% interpolate the values of the numerical solution at the quadrature points
numeric_xq = zeros(options.nquad,1);
for nq=1:options.nquad
    numeric_xq(nq,1) = dot(rhoE_star_num(i_cell:i_cell+1,1),f(nq,:));
end

% compute the cell contribution to the L1 and L2 norms
diff_rhoE = diff_rhoE + dot(wq,exact_xq(:,1)-numeric_xq(:,1))*jac;

% set index_left equal to index_right-1 to speed up next search
index_left = index_right-1;
end % end of loop over cells
diff_rhoE=abs(diff_rhoE);
end % end function

