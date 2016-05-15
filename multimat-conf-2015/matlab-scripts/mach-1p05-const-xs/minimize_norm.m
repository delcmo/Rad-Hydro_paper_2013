% clear all; close all; clc;
function [x_offset] = minimize_norm(numeric_filename,options, exact_filename)
    
% numeric_filename = file storing the numerical solution
% options = quadrature points, interpolation type, index of the numerical
% solution
% exact_filename = used to pass the text files storing the exact solution (x,rho,eps,v,temp)

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
% for both exact and numerical values. No need to do it for mach number
% and density.
for i=2:nb_index_num(1)
    exact_value(:,i) = exact_value(:,i) / exact_value(1,i);
    numeric_value(:,i) = numeric_value(:,i) / numeric_value(1,i);    
end

% get number of cells and return the value
n_cells = size(numeric_value(:,1))-1;
n_cells = n_cells(1);

% get far-field numerical density values
rho_num_left = numeric_value(1,2);
rho_num_right = numeric_value(end,2);
rho_num_middle = 0.5 * (rho_num_left+rho_num_right);
rho_num_third = 1./3. * (rho_num_left+rho_num_right);
rho_num_two_third = 2./3. * (rho_num_left+rho_num_right);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over the cells of the numerical mesh %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable storing the 
index_left = 1; x_offset_third = 0.; x_offset_middle = 0.; x_offset_two_third = 0.;
for i_cell=1:1:n_cells
    % numerical density value belonging to cell 'i_cell'
    node_cell_rho_num = numeric_value(i_cell:i_cell+1,2);
    
    % check if the value 'rho_num_third' belongs to the current cell
    if ( node_cell_rho_num(1,1)<=rho_num_third && node_cell_rho_num(end,1)>rho_num_third )
        rho_num_ref = node_cell_rho_num(1);
        x_num_ref_third = numeric_value(i_cell,1);
        % find density value of the semi-analytical mesh that is the closest from the left
        % to 'rho_num_ref'
        while (exact_value_unique(index_left,2)<rho_num_ref)
            index_left = index_left+1;
        end
        % re-check
        diff_m1 = abs(rho_num_ref-exact_value_unique(index_left-1,2));
        diff = abs(rho_num_ref-exact_value_unique(index_left,2));
        diff_p1 = abs(rho_num_ref-exact_value_unique(index_left+1,2));
        if (diff > diff_m1), index_left=index_left-1; diff=diff_m1; end;
        if (diff > diff_p1), index_left=index_left+1; end;
        % get the corresponing x_exact value:
        x_exact_ref_third = exact_value_unique(index_left,1);
        x_offset_third = abs(x_exact_ref_third-x_num_ref_third);    
    end
    
    % check if the value 'rho_num_middle' belongs to the current cell
    if ( node_cell_rho_num(1,1)<=rho_num_middle && node_cell_rho_num(end,1)>rho_num_middle )
        rho_num_ref = node_cell_rho_num(1);
        x_num_ref_middle = numeric_value(i_cell,1);
        % find density value of the semi-analytical mesh that is the closest from the left
        % to 'rho_num_ref'
        while (exact_value_unique(index_left,2)<rho_num_ref)
            index_left = index_left+1;
        end
        % re-check
        diff_m1 = abs(rho_num_ref-exact_value_unique(index_left-1,2));
        diff = abs(rho_num_ref-exact_value_unique(index_left,2));
        diff_p1 = abs(rho_num_ref-exact_value_unique(index_left+1,2));
        if (diff > diff_m1), index_left=index_left-1; diff=diff_m1; end;
        if (diff > diff_p1), index_left=index_left+1; end;
        % get the corresponing x_exact value:
        x_exact_ref_middle = exact_value_unique(index_left,1);
        x_offset_middle = abs(x_exact_ref_middle-x_num_ref_middle);
    end    
    
    % check if the value 'rho_num_third' belongs to the current cell
    if ( node_cell_rho_num(1,1)<=rho_num_two_third && node_cell_rho_num(end,1)>rho_num_two_third )
        rho_num_ref = node_cell_rho_num(1);
        x_num_ref_two_third = numeric_value(i_cell,1);
        % find density value of the semi-analytical mesh that is the closest from the left
        % to 'rho_num_ref'
        while (exact_value_unique(index_left,2)<rho_num_ref)
            index_left = index_left+1;
        end
        % re-check
        diff_m1 = abs(rho_num_ref-exact_value_unique(index_left-1,2));
        diff = abs(rho_num_ref-exact_value_unique(index_left,2));
        diff_p1 = abs(rho_num_ref-exact_value_unique(index_left+1,2));
        if (diff > diff_m1), index_left=index_left-1; diff=diff_m1; end;
        if (diff > diff_p1), index_left=index_left+1; end;
        % get the corresponing x_exact value:
        x_exact_ref_two_third = exact_value_unique(index_left,1);
        x_offset_two_third = abs(x_exact_ref_two_third-x_num_ref_two_third);  
    end  

end % end of loop over cells
% compute the x_offset value and return it
x_offset = [x_offset_third; x_offset_middle; x_offset_two_third];
end % end function