function a = sound_speed(p, rho)

%HARDCODE, gamma=1.4
gamma = 1.4;

a = sqrt(gamma*p./rho);
