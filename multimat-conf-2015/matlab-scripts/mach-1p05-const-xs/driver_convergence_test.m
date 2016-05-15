clear; close all; clc; format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb_var = 4; % number of variables to consider among [rho,epsilon,mach,T]
options.nquad = 20; % quadrature rule
options.interpolation_type = 'pchip'; % interpolation type: pchip, cubic, spline
indexes = [11;5;1;2;8];
options.index = indexes(1:nb_var+1,1); %(x,rho,eps,mach,mat_temp)
options.plot = true;
options.output = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute mach number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0; % (x,rho,radiation,mach,material temperature)
i=i+1; exact_filename{i}='x_1p05.txt';
i=i+1; exact_filename{i}='rho_1p05.txt';
i=i+1; exact_filename{i}='Er_1p05.txt';
i=i+1; exact_filename{i}='u_1p05.txt';
i=i+1; exact_filename{i}='T_1p05.txt';
eos_param = [5/3;0.14472799784454;1.372e-2]; % [\gamma, C_v, a]
compute_mach_exact( exact_filename, eos_param );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0;
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-10-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-20-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-30-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-40-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-50-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-60-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-70-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-80-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-90-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-100-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-110-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-120-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-130-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-140-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-150-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-160-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-170-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-180-points0.csv';
% i=i+1; numeric_filename{i}='cfl-0.02/mach-1p05-nel-190-points0.csv';
i=i+1; numeric_filename{i}='mach-1p05-nel-100-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-110-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-120-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-130-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-140-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-150-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-160-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-170-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-180-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-190-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-200-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-210-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-220-points0.csv';
%i=i+1; numeric_filename{i}='mach-1p05-nel-230-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-240-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-250-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-260-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-100-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-200-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-300-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-400-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-500-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-600-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-700-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-800-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-900-points0.csv';
% i=i+1; numeric_filename{i}='mach-1p05-nel-1000-points0.csv';

i=0; % (x,rho,radiation,mach,material temperature)
i=i+1; exact_filename{i}='x_1p05.txt';
i=i+1; exact_filename{i}='rho_1p05.txt';
if(i<=nb_var), i=i+1; exact_filename{i}='Er_1p05.txt'; end
if(i<=nb_var), i=i+1; exact_filename{i}='mach_1p05.txt'; end;
if(i<=nb_var), i=i+1; exact_filename{i}='T_1p05.txt'; end;

% options to minimize the mass
options_min.nquad = 4;
options_min.interpolation_type = 'spline'; % interpolation type
options_min.index = [11;5]; %(x,rho,eps,mach,mat_temp)
options_min.eos = [5/3;0.14472799784454;1.372e-2]; % [\gamma, C_v, a]
options_min.output = false;
i=0;
i=i+1; exact_filename_min{i}='x_1p05.txt';
i=i+1; exact_filename_min{i}='rho_1p05.txt';

% % % % % options to minimize the total energy
% % % % options_min.nquad = 20;
% % % % options_min.interpolation_type = 'pchip'; % interpolation type
% % % % options_min.index = [11;6;1]; % [11;5;7;8;1]; % [11;5;1;2;8]; %[11;6]; % [11;5;1;2;8]; %(x,rho,eps,mach,mat_temp)
% % % % options_min.eos = [5/3;0.14472799784454;1.372e-2]; % [\gamma, C_v, a]
% % % % options_min.output = false;
% % % % i=0;
% % % % i=i+1; exact_filename_min{i}='x_1p05.txt';
% % % % i=i+1; exact_filename_min{i}='rho_1p05.txt';
% % % % i=i+1; exact_filename_min{i}='Er_1p05.txt';
% % % % i=i+1; exact_filename_min{i}='u_1p05.txt'; % 'mach_1p05.txt';
% % % % i=i+1; exact_filename_min{i}='T_1p05.txt'; 

x_offset_init = -2.e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ifile=1:length(numeric_filename)
    fprintf('----------------------------------------------------------- \n');
    fprintf('----------------------------------------------------------- \n');
    fprintf('Offset value for the filename %s: \n', numeric_filename{ifile});
    fprintf('----------------------------------------------------------- \n');
    optimset_fminsearch = optimset('TolX',1e-6,'TolFun',1e-4,'Display','off');
    [x_offset(ifile,1), mass_diff, exitflag, output] = fminsearch(@(x_offset_iter)compute_shift_with_mass_conservation(numeric_filename{ifile},options_min,...
                            exact_filename_min, x_offset_iter), x_offset_init, optimset_fminsearch);

% % % %     optimset_fminsearch = optimset('TolX',1e-10,'TolFun',1e-10,'Display','iter');  
% % % %     [x_offset(ifile,1), mass_diff, exitflag, output] = fminsearch(@(x_offset_iter)compute_shift_with_total_energy_conservation(numeric_filename{ifile},options_min,...
% % % %                             exact_filename_min, x_offset_iter), x_offset_init, optimset_fminsearch);

    fprintf('----------------------------------------------------------- \n');                    
    fprintf('Offset value: x_offset = %8.6e. \n', x_offset(ifile,1));                    
    fprintf('Mass difference: mass_diff = %8.6e. \n', mass_diff);
    fprintf('----------------------------------------------------------- \n');                        
                        
%    [L1(ifile,:),L2(ifile,:),n_cells(ifile,:)] = ...
%           post_process_norm(numeric_filename{ifile},...
%                              options,...
%                              exact_filename, x_offset(ifile,1));
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% figure(20)
% x=[log(n_cells(1)),log(n_cells(end))]; 
% y1=[log(L1(1,end));-log(n_cells(end))+log(L1(1,end))+log(n_cells(1))];
% y2=[log(L2(1,end));-0.5*log(n_cells(end))+log(L2(1,end))+0.5*log(n_cells(1))];
% plot(log(n_cells), log(L1(:,end)) ,'-+',x,y1); hold all
% plot(log(n_cells), log(L2(:,end)) ,'-o',x,y2);
% legend('show')

variable_names = {'density'; 'radiation'; 'mach'; 'material temperature'};
for i_var=1:nb_var
    fig_nb=20+i_var;
    fig=figure(fig_nb);
    x=[log(n_cells(1)),log(n_cells(end))]; 
    y1=[-2*log(n_cells(1))+log(L1(end,i_var))+2*log(n_cells(end));log(L1(end,i_var))];
    y2=[-2*log(n_cells(1))+log(L2(end,i_var))+2*log(n_cells(end));log(L2(end,i_var))];
    plot(log(n_cells), log(L1(:,i_var)) ,'-+',x,y1); hold all
    plot(log(n_cells), log(L2(:,i_var)) ,'-o',x,y2);
    legend('L_1 norm', 'slope 2', 'L_2 norm', 'slope 2')
    title(variable_names{i_var})
    saveas(fig,strcat(variable_names{i_var},['_',options.interpolation_type]),'epsc')
end

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
% 
% fprintf(fileID,'------------------------------- \n');
% fprintf(fileID,'------------------------------- \n');
% fprintf(fileID,'Polynomial order of the fitting curve. \n');
% for i_var=1:nb_var
%     p1(i_var,:) = polyfit(log(L1(2:end,i_var)), log(n_cells(2:end,1)), 1);
%     p2(i_var,:) = polyfit(log(L2(2:end,i_var)), log(n_cells(2:end,1)), 1);
%     fprintf(fileID,'%12s \t %12s \t %12s \n', 'Variable number', 'polynomial order L1 norm', 'polynomial order L2 norm');
%     fprintf(fileID,'%d \t %12.8f \t %12.8f \r\n \n',i_var, abs(p1(i_var,1)), abs(p2(i_var,1)));
% end
% fprintf('Polynomial order of the fitting for L1 norm= %8.6e. \n', abs(p1(:,1)));
% fprintf('Polynomial order of the fitting for L2 norm= %8.6e. \n', abs(p2(:,1)));
% 
% open(filename);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % L1
% % L2
% % figure;
% % semilogy(nquad_list,L1);
% % title('L1 norm versus nquad');
% % figure;
% % semilogy(nquad_list,L2);
% % title('L2 norm versus nquad');
% % 
% % show_plot=true;
% % [L1,L2] = post_process_norm(3,visit_filename,csv_file,...
% %     update_vtk_data_with_csv,show_plot);
% fprintf(fileID,'%6.2f %12.8f\r\n',A);
