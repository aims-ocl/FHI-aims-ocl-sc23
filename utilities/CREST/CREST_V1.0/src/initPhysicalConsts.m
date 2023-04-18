%% Initialize physical constants in a structure
function physConsts = initPhysicalConsts

physConsts = struct;
physConsts.k = 1.3807e-023;     % J/K
physConsts.q = 1.60217733e-19;     % C/e
physConsts.e = 1;     % e (hole charge in e)
vacuum_permittivity = 8.854e-12;     % F/m = C/(V-m)
% vacuum permittivity in units of e, V, cm
% [C/(V-m)] / [C/e] = e/(V-m); e/(V-m) * m/cm = e/(V-cm)
physConsts.epsilon_0 = vacuum_permittivity / physConsts.q * 1.e-2;     % e/(V-cm)

end
