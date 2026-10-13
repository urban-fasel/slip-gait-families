function [K, P, rho] = slip_lqr(o, par, Q, R)
%SLIP_LQR  The multiplier rho = max|eig(A)| of the gait o, and the discrete LQR gain K and cost-to-go P
%   of the law u = -K dz for the cost sum dz'Q dz + R u^2.  Riccati iteration, no toolbox needed.
[A, B] = slip_lin(o, par);
rho = max(abs(eig(A)));
P = Q;
for it = 1:20000
    K  = (R + B.'*P*B) \ (B.'*P*A);
    Pn = Q + A.'*P*(A - B*K);  Pn = (Pn + Pn.')/2;
    if norm(Pn - P, 1) <= 1e-14*norm(Pn, 1), P = Pn; break; end
    P = Pn;
end
K = (R + B.'*P*B) \ (B.'*P*A);
end
