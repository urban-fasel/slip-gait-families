function [A, B] = slip_lin(o, par)
%SLIP_LIN  The apex map linearised about the gait o, by complex step, at the energy of the gait, in the
%   states [forward speed; leg angle; leg rate].  A: how an error at one apex maps to the next; its
%   eigenvalues are the multipliers, and the largest in magnitude is rho.  B: how the swing input moves
%   the next apex.
h = 1e-30;  z = o.z;  m = par(1);  g = par(4);  H = 0.5*m*z(2)^2 + m*g*z(1);
Zc = repmat(z, 1, 4);  Zc = Zc + 1i*h*eye(4);
Zc(1,:) = (H - 0.5*m*Zc(2,:).^2)/(m*g);
A = imag(slip_hop(Zc, zeros(1,4), par))/h;  A = A(2:4,2:4);
B = imag(slip_hop(z, 1i*h, par))/h;         B = B(2:4);
end
