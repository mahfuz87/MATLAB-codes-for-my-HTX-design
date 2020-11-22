function [height_mat, h_meth,U, area_ratio,v] = radiator_dim_rectangular_duct(n, n_t, h_air)
% radiator_dim_rectangular_duct finds the
% height_mat - heigth of the duct.
% h_meth - convective heat transfer coefficient of methanol.
% U- convective heat transfer coefficient.
% area_ratio - ratio of area required to dissipate the heat and
% area available in the eliptical duct.
% area_ratio will help to decide if the value of n and n_t is a practical option 
% because to satisfy the cooling demand which must be less than or equal to 1.
% v - velocity of the air which must not cross a certain value for a silent
% cooling experience.
% n - number of pipes.
% n_t - number of turn.
% h_air - convective heat transfer coefficient of air.


% ref1 :  Heat and mass transfer fundamentals and applications
% by Yunus A. Çengel, Afshin J. Ghajar - McGraw-Hill Education (2015)

format long;
h=zeros(31,1);
t=zeros(31,1);
a_in = zeros(31,1);
U=zeros(31,1);
a_req=zeros(31,1);
ratio=zeros(31,1);
v=zeros(31,1);
i=1; %iterative variable

%h_air=34;
q=150; %heat dissipated by processor
LMTD=20;
k_i=79.6; % thermal conductivity of cast iron
t_i=0.001; % Thickness of radiator wall
Q=2.023e-4; % volume flow rate of methanol
%all properties are at bulk/film temperatrure
rho_meth=.425; % density
miu_meth=.995e-5; % dynamic viscosity
k_meth=0.01721; % thermal conductivity of methanol vapor
Cp=1490; %Specific heat

l=n_t*.3; %length of the duct 30cm. n_t is number of turns
w=.03; %width of the duct 3cm.

for d= .002:.0001:.005 %d is the gap between the parallel plates
    A=w*d; % cross sectional area of pipe.
    d_h=2*w*d/(.03+d); %hydraulic diameter of the pipe
    v(i)=Q/(n*A); % velocity
    Re=(rho_meth*v(i)*d_h)/miu_meth; % Reynolds number
    Pr=(miu_meth*Cp)/k_meth; %Prandtl number
    if Re<3000
        Nu=7.54+.03*(d_h/l)*Re*Pr/(1+.016*((d_h/l)*Re*Pr)^(2/3));
        % entry region laminar for duct (ref1 eqN 8-64)
        h_var=Nu*k_meth/d_h; %Nu=h*d/k
    elseif Re>=3000 && Re<8000
        f=-6.38e-13*Re^3+1.17e-8*Re^2-6.69e-5*Re+.147;
        % smooth, parallel-plate channel (ref1 eqN 8-75)
        Nu=((Re-1000)*Pr*(f/8))/(1+12.7*(Pr^(2/3)-1)*(f/8)^.5);
        %Gnielinski eqN (ref1 eqN 8-70)
        h_var=Nu*k_meth/d_h;
    end
    a_in(i) =n*(d+.002+.03+.002)*l;
    U(i)=(h_air^-1+h_var^-1+t_i/k_i)^-1;
    a_req(i)=q/(U(i)*LMTD);
    ratio(i)=a_req(i)/a_in(i);
    h(i)=h_var;
    t(i)=d;
    i=i+1;
end
h_meth=h;
height_mat=t;
area_ratio=vpa(ratio);
end