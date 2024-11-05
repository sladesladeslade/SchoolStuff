function [V2, ierr] = gibbs(R1, R2, R3)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function uses the Gibbs method of orbit determination to
% compute the velocity corresponding to the second of three
% supplied position vectors.
%
% Inputs:
%   R1, R2, R3 - three coplanar geocentric position vectors (km)
%
% Outputs:
%   V2         - the velocity corresponding to R2 (km/s)
%   ierr       - error flag
%                  = 0 if R1, R2, R3 are found to be coplanar
%                  = 1 otherwise
%
% Constants:
%   mu         - gravitational parameter (km^3/s^2)
%   tol        - tolerance for coplanarity check
%
% Variables:
%   r1, r2, r3 - magnitudes of R1, R2, R3 (km)
%   c12, c23, c31 - cross products of R1, R2, R3
%   N, D, S    - intermediate vectors in Gibbs' procedure

global mu
tol = 1e-4;
ierr = 0;

% Magnitudes of R1, R2, and R3
r1 = norm(R1);
r2 = norm(R2);
r3 = norm(R3);

% Cross products among R1, R2, and R3
c12 = cross(R1, R2);
c23 = cross(R2, R3);
c31 = cross(R3, R1);

% Check that R1, R2, and R3 are coplanar; set error flag if not
if abs(dot(R1, c23) / r1 / norm(c23)) > tol
    ierr = 1;
end

% Compute vectors using Gibbs' equations
N = r1 * c23 + r2 * c31 + r3 * c12;           % Equation 5.13
D = c12 + c23 + c31;                           % Equation 5.14
S = R1 * (r2 - r3) + R2 * (r3 - r1) + R3 * (r1 - r2); % Equation 5.21

% Compute velocity V2 corresponding to R2
V2 = sqrt(mu / norm(N) / norm(D)) * (cross(D, R2) / r2 + S);

end % gibbs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
