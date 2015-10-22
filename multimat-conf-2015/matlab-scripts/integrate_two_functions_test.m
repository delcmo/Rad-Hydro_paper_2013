clear all; close all; clc;

% plot for debugging purposes
plot_dgb = true;

% functions are given by a set of coordinates and nodal values (x,y)
% function 1 is a surrogate for Jim's semi analytical GRH
% function 2 is a surrogate for Marco's FEM solution

% extremities of the domain
x_left = 0;
x_right= 2*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defining data 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of points
n_pts_1 = 1000;
% random points
r=rand(n_pts_1,1); 
% random mesh
x_1 = x_left + (sort(r,'ascend'))*(x_right-x_left);
% values
y_1 = cos(2*pi/(x_right-x_left).*x_1);
% check
if plot_dgb, figure(1); plot(x_1,y_1,'.-'); hold all; end

% making poly of that data 
poly_1 = interp1(x_1,y_1,'pchip','pp');
%         pp = interp1(t,y,'linear','pp');
f_1 = @(x) ppval(poly_1,x);
if plot_dgb, x=linspace(x_left,x_right,10); plot(x,f_1(x)); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defining data 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of points
n_pts_2 = 11; n_elem_2 = n_pts_2-1;
% mesh
x_2 = linspace(x_left,x_right,n_pts_2);
% mesh size 
dx = (x_right-x_left)/n_elem_2;
% values
y_2 = cos(2*pi/(x_right-x_left).*x_2);
% check
if plot_dgb, plot(x_2,y_2,'.-'); hold all; end

% making poly of that data 
poly_2 = interp1(x_2,y_2,'linear','pp');
f_2 = @(x) ppval(poly_2,x);
if plot_dgb, figure(2); x=linspace(x_left,x_right,10); plot(x,f_2(x)); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reset norm
norm = 0;
% select GL quadrature
n_quad = 2;
[xq,wq] = GLNodeWt(n_quad);

% loop over FEM mesh
for iel=1:n_elem_2
    % beginning/end of mesh interval
    x_beg = x_left + (iel-1)*dx;
    x_end = x_beg  + dx;
    % quad points on the physical mesh
    x_qp = (x_beg+x_end)/2 + dx*xq;
    % difference squared
    diff = (f_2(x_qp)-f_1(x_qp)).^2;
    % increment norm
    norm = norm + dot(wq*dx/2, diff);     
end
norm = sqrt(norm);
fprintf('norm %g using %d cells \n',norm,n_elem_2);