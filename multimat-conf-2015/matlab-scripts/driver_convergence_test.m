clear; close all; clc; format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb_var = 4; % number of variables to consider among [rho,epsilon,mach,T]
options.nquad = 5; % quadrature rule
options.interpolation_type = 'spline'; % interpolation type
indexes = [11;5;1;2;8];
options.index = indexes(1:nb_var+1,1); %(x,rho,eps,mach,mat_temp)
options.plot = false;
options.output = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0;
%i=i+1; numeric_filename{i}='mach-3-nel-25-points0.csv';
%i=i+1; numeric_filename{i}='mach-3-nel-50-points0.csv';
i=i+1; numeric_filename{i}='mach-3-nel-100-points0.csv';
i=i+1; numeric_filename{i}='mach-3-nel-200-points0.csv';
i=i+1; numeric_filename{i}='mach-3-nel-400-points0.csv';
i=i+1; numeric_filename{i}='mach-3-nel-800-points0.csv';
i=i+1; numeric_filename{i}='mach-3-nel-1600-points0.csv';
i=i+1; numeric_filename{i}='mach-3-nel-3200-points0.csv';
%i=i+1; numeric_filename{i}='mach-3-nel-6400-points0.csv';
i=0; % (x,rho,radiation,mach,material temperature)
i=i+1; exact_filename{i}='data_x.dat';
i=i+1; exact_filename{i}='data_Density.dat';
if(i<=nb_var), i=i+1; exact_filename{i}='data_RED.dat'; end
if(i<=nb_var), i=i+1; exact_filename{i}='data_Mach.dat'; end;
if(i<=nb_var), i=i+1; exact_filename{i}='data_Temp.dat'; end;

options_min.interpolation_type = 'spline'; % interpolation type
options_min.index = [11;5]; %(x,rho,eps,mach,mat_temp)
options_min.output = false;
i=0;
i=i+1; exact_filename_min{i}='data_x.dat';
i=i+1; exact_filename_min{i}='data_RED.dat';

% fprintf('----------------------------------------------------------- \n');
% fprintf('----------------------------------------------------------- \n');
% fprintf('Offset value for the filename %s: \n', numeric_filename{end});
% fprintf('----------------------------------------------------------- \n');    
% x_offset = minimize_norm(numeric_filename{end},...
%                        options_min,...
%                        exact_filename_min);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ifile=1:length(numeric_filename)
    fprintf('----------------------------------------------------------- \n');
    fprintf('----------------------------------------------------------- \n');
    fprintf('Offset value for the filename %s: \n', numeric_filename{ifile});
    fprintf('----------------------------------------------------------- \n');    
    x_offset = minimize_norm(numeric_filename{ifile},...
                           options_min,...
                           exact_filename_min);    
    fprintf('Offset value: x_offset = %8.6e. \n', x_offset);
    fprintf('----------------------------------------------------------- \n');        
    fprintf('L1 and L2 norms of the error for file %s: \n', numeric_filename{ifile});    
    fprintf('----------------------------------------------------------- \n');        
    [L1(ifile,:),L2(ifile,:),n_cells(ifile,:)] = ...
            post_process_norm(numeric_filename{ifile},...
                              options,...
                              exact_filename, x_offset(2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(20)
plot(log(n_cells), log(L1) ,'-+'); hold all
plot(log(n_cells), log(L2) ,'-o');
legend(['L1';'L2']);
 
% if length(nquad_list) > 1
% figure(21)
% for iq=1:length(nquad_list)
%  plot(log(n_cells(:,iq)), log(L1(:,iq)) ,'-+'); hold all
% end
% figure(22)
% for iq=1:length(nquad_list)
%  plot(log(n_cells(:,iq)), log(L2(:,iq)) ,'-+'); hold all
% end
filename = 'error_norms.txt';
if exist(filename, 'file'), delete(filename); end

fileID = fopen(filename,'w'); 
for i_var=1:nb_var
    fprintf(fileID,'------------------------------- \n');
    fprintf(fileID,'------------------------------- \n');
    fprintf(fileID,'Variable number: %d\n',i_var);
    L1_rate = log(L1(1:end-1,i_var)./L1(2:end,i_var)) / log(2);
    L2_rate = log(L2(1:end-1,i_var)./L2(2:end,i_var)) / log(2);
    fprintf(fileID,'%12s \t %12s \t %12s \n', 'nb_cells', 'L1_norms','L2_norms');
    fprintf(fileID,'%d \t %12.8f \t %12.8f \r\n \n', [n_cells'; L1(:,i_var)'; L2(:,i_var)']);
    fprintf(fileID,'%12s \t %12s\n','L1_rate','L2_rate');
    fprintf(fileID,'%12.8f \t %12.8f\r\n \n', [L1_rate'; L2_rate']);
end

fprintf(fileID,'------------------------------- \n');
fprintf(fileID,'------------------------------- \n');
fprintf(fileID,'Polynomial order of the fitting curve. \n');
for i_var=1:nb_var
    p1(i_var,:) = polyfit(log(L1(2:end,i_var)), log(n_cells(2:end,1)), 1);
    p2(i_var,:) = polyfit(log(L2(2:end,i_var)), log(n_cells(2:end,1)), 1);
    fprintf(fileID,'%12s \t %12s \t %12s \n', 'Variable number', 'polynomial order L1 norm', 'polynomial order L2 norm');
    fprintf(fileID,'%d \t %12.8f \t %12.8f \r\n \n',i_var, abs(p1(i_var,1)), abs(p2(i_var,1)));
end
fprintf('Polynomial order of the fitting for L1 norm= %8.6e. \n', abs(p1(:,1)));
fprintf('Polynomial order of the fitting for L2 norm= %8.6e. \n', abs(p2(:,1)));

open(filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1
% L2
% figure;
% semilogy(nquad_list,L1);
% title('L1 norm versus nquad');
% figure;
% semilogy(nquad_list,L2);
% title('L2 norm versus nquad');
% 
% show_plot=true;
% [L1,L2] = post_process_norm(3,visit_filename,csv_file,...
%     update_vtk_data_with_csv,show_plot);
fprintf(fileID,'%6.2f %12.8f\r\n',A);
