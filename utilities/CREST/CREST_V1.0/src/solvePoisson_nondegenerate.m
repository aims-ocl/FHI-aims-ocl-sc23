%% Analytical solution of Poisson's equation, applies only to the case that
% neither the bulk nor the surface Ef are within degenerate distance from
% either band edge. Yields only total space charge as a function of surface
% potential.
% Based on Sze 2006 chapter 4 pp. 200-203
function Qsc = solvePoisson_nondegenerate(psi_s, T, Ndop, Nv, Nc, epsilon, Eg)

% Constants
physConsts = initPhysicalConsts;
q = physConsts.q;     % C/e
e = physConsts.e;     % e
k = physConsts.k;     % J/K
epsilon_0 = physConsts.epsilon_0;     % e/(V-cm)

kT = k*T / q;     % J / C/e = CV / C/e = eV
beta = e / kT;     % e/eV

% Bulk electronic structure, all referenced to Ev (VBM)
Ev_bulk = 0;
Ec_bulk = Eg;

% Find Ef as determined by neutrality in the bulk
Ef = findEMBulkEf(T, Ndop, Nv, Nc, Ev_bulk, Ec_bulk);

% Bulk carrier concentrations
n0 = Nc * fermi_fo(1/2, beta * (Ef - Ec_bulk));     % 1/cm^3
p0 = Nv * fermi_fo(1/2, beta * (Ev_bulk - Ef));     % 1/cm^3

% Debye length for holes in material
LDp = sqrt((kT * epsilon * epsilon_0) / (p0 * e^2));     % sqrt((eV * e/V-cm) / 1/cm^3 * e^2) = sqrt(e^2/cm / e^2/cm^3) = sqrt(e^2/cm * cm^3/e^2) = cm
% F function dependent on surface potential
F = sqrt(...
    exp(-beta.*psi_s) + beta.*psi_s - 1 + ...
    n0/p0 .* (exp(beta.*psi_s) - beta.*psi_s - 1) ...
    );     % [beta*psi_s] = e/eV * V = 1
if min(F) < 0
    err = MException('MATLAB:solvePoisson_nondegenerate:NegativeF', ...
        'Something is wrong, F shouldn''t have been negative');
    throw(err);
end

% Total space charge per unit area
Qsc = -sign(psi_s) .* sqrt(2) .* epsilon .* epsilon_0 .* kT / (e .* LDp) .* F;     % e/V-cm * eV / (e * cm) = e^2/cm / e-cm = e/cm^2

end