%% Numerical solution of Poisson's equation
% Based on Sze 2006 chapter 4 pp. 200-203
function [Qsc, x, V, E, Qnet, trunc_est_from, errmsg] = solvePoisson_numeric(psi_s, ...
    T, Ndop, Nv, Nc, epsilon, Eg, degeneracyLimit_kT, Qsc_guess, varargin)
% For errors and estimations
trunc_est_from = [];
errmsg = '';

% Parse variable input
force_Qsc = false;
optionIndex = find(strcmp('force', varargin), 1);
if ~isempty(optionIndex)
    force_Qsc = true;
    varargin(optionIndex) = [];
end
output_grid_numpoints = 5001;
optionIndex = find(strcmp('outgridnpts', varargin), 1);
if ~isempty(optionIndex)
    output_grid_numpoints = varargin{optionIndex+1};
    varargin(optionIndex + [0 1]) = [];
end

% Constants
physConsts = initPhysicalConsts;
q = physConsts.q;     % C/e
e = physConsts.e;     % e
k = physConsts.k;     % J/K
epsilon_0 = physConsts.epsilon_0;     % e/(V-cm)

kT = k*T / q;     % J / C/e = CV / C/e = eV
beta = e / kT;     % e/eV

% Bulk electronic structure, all referenced to bulk Ev (VBM)
Ev_bulk = 0;
Ec_bulk = Eg;

% Find Ef as determined by neutrality in the bulk
Ef = findEMBulkEf(T, Ndop, Nv, Nc, Ev_bulk, Ec_bulk);

% Bulk carrier concentrations
n0 = Nc * fermi_fo(1/2, beta * (Ef - Ec_bulk));     % 1/cm^3
p0 = Nv * fermi_fo(1/2, beta * (Ev_bulk - Ef));     % 1/cm^3

% Estimate range for solver function from depletion region approximation
% and estimated Qsc given; if in depletion/inversion regime, use majority
% carrier concentration for estimation, if in accumulation, use minority
if Qsc_guess > 0
    % bands bent up ==> use electrons
    carrier_concentration = n0;
elseif Qsc_guess < 0
    % bands bent down ==> use holes
    carrier_concentration = p0;
else
    % Bands not bent, so it doesn't matter - the guess will be 0 anyway
    carrier_concentration = 1;
end
% Find formal depletion region depth
depletion_region_depth = abs(Qsc_guess / carrier_concentration);     % (e/cm^2) / (e/cm^3) = cm
% Set grid length. Min. 1 um (= 1e-4 cm) to ensure macroscopic dimensions
grid_total_depth = max(1e-4, 20*depletion_region_depth);     % cm
% Limit grid length to > 10 um (= 1e-3 cm) to ensure numerical convergence
grid_depth_max = 0.001;     % cm
if grid_total_depth > grid_depth_max; grid_total_depth = grid_depth_max; end

% Construct x squared grid (finer near 0)
xgrid = (linspace(0, 1, 1501)').^2 * grid_total_depth;     % cm

% Boundary values and initial guess
Va_bc = psi_s;
if force_Qsc
    % Find field
    dVdxa_bc = Qsc_guess / (epsilon_0 * epsilon);    % e/cm^2 / (e/(V-cm) = V/cm
else
    dVdxa_bc = [];
end
yinit = bvpinit(xgrid, @(x) y_initialGuess(x, 0, psi_s, depletion_region_depth));

% Error/warning checking
lastwarn('');
% Solve
try
    sol = bvp4c(@(x, y) diffVdVdx(x, y, n0, p0, beta, e, epsilon, epsilon_0, ...
        Nc, Nv, Ec_bulk, Ev_bulk, Ef, degeneracyLimit_kT / kT), ...
        @(ya, yb) PoissonEqBC(ya, yb, Va_bc, dVdxa_bc, ...
        T, Ndop, Nv, Nc, Eg, epsilon, epsilon_0), ...
        yinit);
catch err
    % Save error in errmsg
    errmsg = getReport(err);
end
% Warning checking
if ~strcmp(lastwarn, '')
    errmsg = ['Warning: ' lastwarn];
end

% Threshold against which to check results
Vthr = abs(psi_s / 10000);

if strcmp(errmsg, '')
    % Extract results (if no error occurred)
    xrange = linspace(0, 1, output_grid_numpoints)' * grid_total_depth;     % cm
    x = -flipud(xrange);
    y = deval(sol, xrange);
    % Potential
    V = flipud(real(y(1, :))');
    Vimaginary = flipud(imag(y(1, :))');
    % Field
    E = -flipud(real(y(2, :))');
    Eimaginary = -flipud(imag(y(2, :))');
    % Differentiate field for net volumetric charge
    Qnet = [0; diff(-E)/diff(x(1:2)) * (epsilon_0 * epsilon)];    % V/cm / cm * e/(V-cm) = V/cm^2 * e/(V-cm) = e/cm^3
    % Total charge in SCR, from last point in field
    Qsc = epsilon_0 * epsilon * (-E(end)); % e/(V-cm) * V/cm = e/cm^2
    
    % Check that results make physical sense:
    % V must be bounded by 0 and psi_s
    % V must be monotonically increasing/descreasing (within numerical error)
    % V must start at 0 if the solution range was not too much truncated
    % The imaginary parts of V and E must be in the range of numerical error
    Ethr = abs(max(abs(E)) / 10000);
    if max(V) > abs(psi_s) + Vthr || (psi_s > 0 && min(V) < -Vthr) || (psi_s < 0 && max(V) > Vthr) || ...
            (psi_s > 0 && any(diff(V) < -Vthr)) || (psi_s < 0 && any(diff(V) > Vthr)) || ...
            (grid_total_depth > 15*depletion_region_depth  && abs(V(1)) > Vthr) || ...
            max(abs(Vimaginary)) > Vthr || max(abs(Eimaginary)) > Ethr
        errmsg = 'Potential obtained is not physical';
    end
end

% If an error occurred, set all results to initial values so that later an
% approximation can be constructed.
if ~strcmp(errmsg, '')
    grid_total_depth = 0;     % cm
    x = 0;
    % Potential
    V = psi_s;
    % Field
    E = -Qsc_guess / (epsilon_0 * epsilon);    % e/cm^2 / (e/(V-cm) = V/cm
    % Net volumetric charge
    Qnet = 0;     % e/cm^3
    % Total charge in SCR
    Qsc = Qsc_guess;     % e/cm^2
end

% Check to see if the solution actually reached 0. If not, this may be due
% to truncation of the grid, or the numeric solution failed. Either way,
% complete the solution using a parabola with extremum at height 0 and
% continuity of differential at stitching point.
if abs(V(1)) > Vthr
    trunc_est_from = -grid_total_depth;
    y1_1 = V(1);
    y2_1 = -E(1);
    x0 = -2*y1_1/y2_1;
    par_a = y2_1^2 / (4*y1_1);
    xrange_completion = linspace(0, 1, output_grid_numpoints)' * x0;     % cm
    y1_completion = par_a .* (xrange_completion - x0).^2;
    y2_completion = 2*par_a .* (xrange_completion - x0);
    % Complete solutions
    x_completion = -flipud(grid_total_depth + xrange_completion);
    x = [x_completion; x(2:end)];
    % Potential
    V_completion = flipud(y1_completion);
    V = [V_completion; V(2:end)];
    % Field
    E_completion = -flipud(y2_completion);
    E = [E_completion; E(2:end)];
    % Net volumetric charge
    Qnet_completion = [0; diff(-E_completion)/diff(x_completion(1:2)) * (epsilon_0 * epsilon)];
    % V/cm / cm * e/(V-cm) = V/cm^2 * e/(V-cm) = e/cm^3
    Qnet = [Qnet_completion; Qnet(2:end)];
end

end


%% Second-order ODE bvp problem function - Poisson's equation
% Does NOT assume nondegenerate doping
function dydx = diffVdVdx(x, y, n0, p0, beta, e, epsilon, epsilon_0, ...
    Nc, Nv, Ec_bulk, Ev_bulk, Ef, degeneracyLimit)
% Unpack
V = y(1);     % V
dVdx = y(2);     % V/cm

% Check if the value of V is reasonable in view of the gap
if abs(e * V) > abs(Ec_bulk - Ev_bulk) * 10
    err = MException('MATLAB:solvePoisson_numeric:CrazyU', ...
        sprintf(['Reached electrotatic energy value of U = %0.10f, ' ...
        ', which is over ten times the SC gap Ec - Ev = %0.10f'], ...
        e * V, Ec_bulk - Ev_bulk));
    throw(err);
end

% x-dependent carrier concentrations. The band edges are shifted by V,
% down for V > 0:
Ec_x = Ec_bulk - e * V;
Ev_x = Ev_bulk - e * V;

% Use the same expression as for bulk situations. Note that if
% V = 0 this returns to the bulk situation.
if Ec_x - Ef < degeneracyLimit
    n = Nc * fermi_fo(1/2, beta * (Ef - Ec_x));     % 1/cm^3
else
    n = Nc * exp(beta * (Ef - Ec_x));     % 1/cm^3
end
if Ef - Ev_x < degeneracyLimit
    p = Nv * fermi_fo(1/2, beta * (Ev_x - Ef));     % 1/cm^3
else
    p = Nv * exp(beta * (Ev_x - Ef));     % 1/cm^3
end

% d2Vdx2 = -e / (epsilon * epsilon_0) .* (Nd - Na + p - n);
d2Vdx2 = -e / (epsilon * epsilon_0) .* (n0 - p0 + p - n);
% e / e/(V-cm) * 1/cm^3 = V-cm / cm^3 = V/cm^2

dydx = [dVdx; d2Vdx2];
end


%% Second-order ODE bvp boundary conditions function
% Nondegenerate-conditions relation between potential and
% derivative at b, potential at a, if derivative at a also provided then
% combined potential+derivative at a
function res = PoissonEqBC(ya, yb, Va_bc, dVdxa_bc, ...
    T, Ndop, Nv, Nc, Eg, epsilon, epsilon_0)
% Unpack
Va = ya(1);
dVdxa = ya(2);
Vb = yb(1);
dVdxb = yb(2);

% Find field assuming nondegenerate conditions at the edge of the range
psi_s_at_b = Vb;
Qsc_at_b = solvePoisson_nondegenerate(psi_s_at_b, T, Ndop, Nv, Nc, epsilon, Eg);
dVdxb_bc = Qsc_at_b / (epsilon_0 * epsilon);    % e/cm^2 / (e/(V-cm) = V/cm

% If dVdx bc given for a, find field at beginning of range and create a
% composite boundary condition
if ~isempty(dVdxa_bc)
    % Create composite bc for point a
    a_bc = norm([Va - Va_bc, dVdxa - dVdxa_bc]);
else
    % The given Qsc was just a guess
    a_bc = Va - Va_bc;
end

res = [a_bc; dVdxb - dVdxb_bc];
end


%% Second-order ODE bvp initial-guess construction function. Take simply
% a parabola beginning at (xa, Vinit) and reaching its extremum at (xa +
% xd, 0), and thereafter 0
function yinit = y_initialGuess(x, xa, Vinit, xd)

x_extr = xa + xd;

if x < x_extr
    % Find parabola that goes through (xa, Vinit) and has its extremum at
    % (x_extr, 0).
    coeffs = getParabolaFrom3P([x_extr - xd, Vinit], ...
        [x_extr, 0], [x_extr + xd, Vinit]);
    V = coeffs.a2 * x^2 + coeffs.a1 * x + coeffs.a0;
    dVdx = 2 * coeffs.a2 * x + coeffs.a1;
else
    V = 0;
    dVdx = 0;
end

yinit = [V; dVdx];

end