#This file is the "main function"
#type "source main.m" in the prompt to run

clear;

#d = 8 variant (functions with variable d)
#d = 4;
d = 8;
ds = [2,2,2]; #this means that the first subsystem d=2 and second d=4
rho0 = get_w();
rho1 = get_rho_1_d(d);

css = gilbert_d(rho0, rho1, ds)
#css = gilbert_d_bures(rho0, rho1, ds)
