#Bianka is acknowledged hereby by Omer for reminding us all about the sqrt().
#Thank you Bianka.

#This file is the "main function"
#type "source main.m" in the prompt to run

#d = 2 variant (functions with fixed d)
#rho_0 = get_rho_0();
#rho_1 = get_rho_1();
#css = gilbert(rho_0, rho_1);

#d = 8 variant (functions with variable d)
d = 8;
ds = [2,4]; #this means that the first subsystem d=2 and second d=4
rho0 = get_rho_0_d(d);
rho1 = get_rho_1_d(d);
css = gilbert_d(rho0, rho1, ds)
