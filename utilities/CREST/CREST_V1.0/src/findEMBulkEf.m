%% Procedure that finds the Ef based on the Effective Mass model for bulk
% Based on Sze 2006.
function Ef_neutrality = findEMBulkEf(T, Ndop, Nv, Nc, Ev, Ec)

% Constants
physConsts = initPhysicalConsts;
q = physConsts.q;     % C/e
e = physConsts.e;     % e
k = physConsts.k;     % J/K

kT = k*T / q;     % J / C/e = CV / C/e = eV
beta = e / kT;     % e/eV

% Nd and Na - assume only one type of doping
if Ndop >=0
    Nd = abs(Ndop);
    Na = 0;
else
    Nd = 0;
    Na = abs(Ndop);
end

% Find Ef of the neutral bulk as the point where the
% concentration of electrons + donor contributions + holes = 0
initial_guess = (Ec - Ev) / 2;
Ef_neutrality = fzero(@(Ef) findNetCharge(Ef, Nc, Nv, Ec, Ev, Nd, Na, beta), initial_guess);

end


%% Helper function to provide net charge given bulk properties, band edges and doping density
function netCharge = findNetCharge(Ef, Nc, Nv, Ec, Ev, Nd, Na, beta)

% Concentration of electrons
n = Nc * fermi_fo(1/2, beta * (Ef - Ec));
% Concentration of holes
p = Nv * fermi_fo(1/2, beta * (Ev - Ef));

% Net concentration of charge under assumption of full ionization. Ionized
% donors are positive, ionized acceptors negative, holes positive and
% electrons negative.
netCharge = Nd - Na + p - n;

end