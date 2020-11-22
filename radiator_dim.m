% This snip was used to call the functions for elliptical and rectangular
% duct to write the output variables into excel file

[h_air] = h_air_calc(10);
[height_mat, h_meth,U, area_ratio,v] = radiator_dim_elliptical_duct(5,2,h_air);
filename = 'elliptical_duct.xlsx';
area_ratio=double(area_ratio);
T = table(height_mat, h_meth, U, area_ratio,v);
writetable(T,filename);
[height_mat, h_meth,U, area_ratio,v] = radiator_dim_rectangular_duct(5,2,h_air);
filename = 'rectangular_duct.xlsx';
area_ratio=double(area_ratio);
T = table(height_mat, h_meth, U, area_ratio,v);
writetable(T,filename);