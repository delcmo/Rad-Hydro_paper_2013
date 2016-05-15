function compute_mach_exact( exact_filename, eos_param )
% This function computes the mach number for the exact solution from
% the velocity, temperature and EOS parameters.

% load data from exact and numeric filenames
nb_exact_file = size(exact_filename);
for ifile=1:nb_exact_file(2)
   file_id = fopen(exact_filename{ifile});
   exact_value(:,ifile) = textread(exact_filename{ifile}, '%f');
   fclose(file_id);
end

u=exact_value(:,4);
T=exact_value(:,5);
% eos parameters
gamma=eos_param(1);
Cv=eos_param(2);

% compute speed of sound: c^2 = gamma*(gamma-1)*Cv*T
c2=gamma*(gamma-1)*Cv*T;

% compute mach number
mach=abs(u)./sqrt(c2);

% save mach number in a file of the same format of u
filename = 'mach_1p05.txt';
if exist(filename, 'file'), delete(filename); end

fileID = fopen(filename,'w');
fprintf(fileID,'%12.8f\n',mach');
end

