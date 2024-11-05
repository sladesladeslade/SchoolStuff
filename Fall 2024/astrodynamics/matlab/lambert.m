function [V1, V2] = lambert(R1, R2, t, string)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function solves Lambert's problem.
%
% Inputs:
%   R1, R2  - initial and final position vectors (km)
%   t       - time of flight from R1 to R2 (s)
%   string  - 'pro' for prograde or 'retro' for retrograde orbit
%
% Outputs:
%   V1, V2  - initial and final velocity vectors (km/s)
%
% Constants:
%   mu      - gravitational parameter (km^3/s^2)
%
% Variables:
%   r1, r2  - magnitudes of R1 and R2
%   c12     - cross product of R1 and R2
%   theta   - angle between R1 and R2
%   A       - constant (Equation 5.35)
%   z       - alpha * x^2, with alpha as the reciprocal of the semimajor axis
%   y(z)    - function of z (Equation 5.38)
%   F(z,t)  - function of z and t (Equation 5.40)
%   dFdz(z) - derivative of F(z,t) (Equation 5.43)
%   ratio   - F/dFdz
%   tol     - tolerance for convergence precision
%   nmax    - max iterations for Newton's method
%   f, g    - Lagrange coefficients
%   gdot    - derivative of g
%   C(z), S(z) - Stumpff functions

global mu
r1 = norm(R1);
r2 = norm(R2);
c12 = cross(R1, R2);
theta = acos(dot(R1, R2) / (r1 * r2));

% Determine if orbit is prograde or retrograde
if nargin < 4 || (~strcmp(string, 'retro') && ~strcmp(string, 'pro'))
    string = 'pro';
    fprintf('\n ** Prograde trajectory assumed.\n');
end
if strcmp(string, 'pro') && c12(3) <= 0
    theta = 2 * pi - theta;
elseif strcmp(string, 'retro') && c12(3) >= 0
    theta = 2 * pi - theta;
end

% Equation 5.35
A = sin(theta) * sqrt(r1 * r2 / (1 - cos(theta)));

% Determine starting value of z for Equation 5.45
z = -100;
while F(z, t) < 0
    z = z + 0.1;
end

% Set tolerance and max iterations
tol = 1.e-8;
nmax = 5000;

% Iterate on Equation 5.45 to find z
ratio = 1;
n = 0;
while (abs(ratio) > tol) && (n <= nmax)
    n = n + 1;
    ratio = F(z, t) / dFdz(z);
    z = z - ratio;
end

% Report if iteration limit exceeded
if n >= nmax
    fprintf('\n\n ** Number of iterations exceeds %g in ''lambert'' \n\n', nmax);
end

% Compute Lagrange coefficients
f = 1 - y(z) / r1;                       % Equation 5.46a
g = A * sqrt(y(z) / mu);                  % Equation 5.46b
gdot = 1 - y(z) / r2;                     % Equation 5.46d

% Compute initial and final velocities
V1 = (R2 - f * R1) / g;                   % Equation 5.28
V2 = (gdot * R2 - R1) / g;                % Equation 5.29
return

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Subfunctions used in the main body
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Function y(z) (Equation 5.38)
function dum = y(z)
    dum = r1 + r2 + A * (z * S(z) - 1) / sqrt(C(z));
end

% Function F(z,t) (Equation 5.40)
function dum = F(z, t)
    dum = (y(z) / C(z))^1.5 * S(z) + A * sqrt(y(z)) - sqrt(mu) * t;
end

% Derivative dFdz(z) (Equation 5.43)
function dum = dFdz(z)
    if z == 0
        dum = sqrt(2) / 40 * y(0)^1.5 + A / 8 * (sqrt(y(0)) + A * sqrt(1 / (2 * y(0))));
    else
        dum = (y(z) / C(z))^1.5 * (1 / (2 * z) * (C(z) - 3 * S(z) / (2 * C(z))) ...
            + 3 * S(z)^2 / (4 * C(z))) + A / 8 * (3 * S(z) / C(z) * sqrt(y(z)) ...
            + A * sqrt(C(z) / y(z)));
    end
end

% Stumpff function C(z)
function dum = C(z)
    dum = stumpC(z);
end

% Stumpff function S(z)
function dum = S(z)
    dum = stumpS(z);
end

end % lambert
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~