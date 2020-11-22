function [h_air] = h_air_calc(v)

% h_air_calc calculates convective heat transfer coefficient of air given
% the velocity of air - v.
%
% All necessary air properties are embedded in the code in SATP and 60%
% humidity. They can be adjusted in the code if necessary.

% Reference book for necessary equations is 
% Heat and mass transfer fundamentals and applications
% by Yunus A. Çengel, Afshin J. Ghajar - McGraw-Hill Education (2015)


rho_a=1.18; % density of air
miu_a=1.86e-5; % dynamic viscosity
d=.005; % average spacing between two ducts through where air will flow. 5mm
w=.04; % width of the radiator. 4cm
l=.3; % length of the radiator. 30cm
d_h=2*l*d/(l+d); % hydraulic diameter of the pipe with rectangle cross section
k=.026; % thermal conductivity of air
Re=vpa(rho_a*v*d_h/miu_a);
Pr=0.72;
if Re<3000
    Nu=7.54+.03*(d_h/w)*Re*Pr/(1+.016*((d_h/w)*Re*Pr)^(2/3));
    % entry region laminar for duct (eqN 8-64)
    h_air=Nu*k/d_h; % Nu=h*d/k
elseif Re>=3000 && Re<8000
    f=-6.38e-13*Re^3+1.17e-8*Re^2-6.69e-5*Re+.147;
    % Smooth, parallel-plate channel (eqN 8-75)
    Nu=((Re-1000)*Pr*(f/8))/(1+12.7*(Pr^(2/3)-1)*(f/8)^.5);
    % Gnielinski eqation (eqN 8-70)
    h_air=Nu*k/d_h;
else
    f=(0.79*log(Re)-1.64)^-2; % smooth tubes (eqN 8-65)
    Nu=((Re-1000)*Pr*(f/8))/(1+12.7*(Pr^(2/3)-1)*(f/8)^.5);
    % Gnielinski eqation (eqN 8-70)
    h_air=Nu*k/d_h;
 end
end