clear all; clc;

% eos gamma
g=5/3;   % not sure
% upstream Mach
M=1.05; 
M=3; 
% radiation constant
aR = 1.3720172e-2;
cv = 0.14472799784454;

% a-dimensionalized pre-shock values
T0=0.1; % not 1.0 but 0.1 keV = 100 eV  
rho0=1; 

% pressure
press0 = (g-1)*rho0*cv*T0
% radiatio energy
epsilon0 = aR * T0^4
% sound speed. the initial radiation energy is very small, so 
% the sound speed and the radiation-modified sound speed should be very
% close
cm0 = sqrt( (g*press0 + 4*epsilon0/9)/rho0 )
c0 = sqrt( g*press0/rho0 )
% velocity
% u0 = M * c0
u0 = M % Marco picks u0 as M. why?
u0 = M *c0
% 
% fprintf('this P0 should be about 8.8e-5 if this script is correct\n');
fprintf('this P0 should be about 8.5e-5 if this script is correct (This is what I have in my code and this is what Jim uses. \n');
P0 = aR*T0^4 / (rho0*c0^2)

% From my Moose code, the reference values are computed:
% Reference value for density: 0.9999780175
% Reference value for material temperature: 0.09999968746
% Reference value for material sound speed: 0.12681025


% initial guess for P0 paramter (goal is to recover 8.8e-5 from paper)
% P0=1e-4;
% initial guess for post shock T
T1 = ((1-g+2*g*M^2)*(2+(g-1)*M^2))/((g+1)^2*M^2)
% if P0>1
%     T1= (8/7*(M^2/(4/9*P0)-1))^0.25;
% end

% equations after eq.(12) in Lowrie/Rauenzahn
f1 = @(T) 3*(g+1)*(T-1)-P0*g*(g-1)*(7+T.^4);
f2 = @(T) 12*(g-1)^2*T.*(3+g*P0*(1+7*T.^4));
% eq.(12) in Lowrie/Rauenzahn
rho = @(T) (f1(T) + sqrt(f1(T).^2 + f2(T)))./(6*(g-1)*T);
% eq.(13) in Lowrie/Rauenzahn
residual = @(r,T) 3*r*(r*T-1) + g*P0*r*(T.^4-1) - 3*g*(r-1)*M^2;
dresdT = @(r,T) 3*r*r + 4*g*P0*r*T.^3 ;

% solve eq.(12) and eq.(13) for T1 and rho1:
for iter=1:100
    rho1 = rho(T1/T0);
    fprintf('iter %d, T1=%g, rho1=%g \n',iter,T1,rho1);
    res = residual(rho1,T1/T0);
    if abs(res)<1e-10
        break
    end
    slope = dresdT(rho1,T1/T0);
    T1= T1 - res/slope*T0;
end
        

% In Moose: T1=0.1494666288 if T0=0.1
% Marco, you must have a typo: not 0.1494... but 0.10494... (missing a 0)

% pressure
press1 = (g-1)*rho1*cv*T1
% radiatio energy
epsilon1 = aR * T1^4
% sound speed
% post shock, the soupdspeed and the radiation-modified soundspeed should be different
cm1 = sqrt( (g*press1 + 4*epsilon1/9)/rho1 )
c1 = sqrt( g*press1/rho1 )
% conservation of momentum gives the downstream velocity
u1 = rho0 * u0 / rho1
% finally, we can get a post-shock Mach
M1 = u1
M1 = u1/c1
M1_ = u1/cm1

radiation_energy_ratio=epsilon1/epsilon0
radiation_T_ratio=(epsilon1/epsilon0)^0.25
Mach_ratio=M1_/M

% attempt at getting sigma
L=1; c=3e8;
sigma = 1e6*cm0/c/L
sigma = 1e6*cm1/c/L
