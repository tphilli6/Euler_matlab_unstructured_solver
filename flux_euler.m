function flux = flux_euler( ul, ur, normal)


fl = euler_flux(ul, normal);
fr = euler_flux(ul, normal);

flux = (fl+fr)/2;