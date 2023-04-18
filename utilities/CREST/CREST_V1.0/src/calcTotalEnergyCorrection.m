%% Calculates total energy correction due to SCR "contraction"
function Utot_correction = calcTotalEnergyCorrection(epsilon, z_SCR, E_SCR, Qsheet, z_d, sheetZ)

% Constants
physConsts = initPhysicalConsts;
epsilon_0 = physConsts.epsilon_0;     % e/(V-cm)

% Total energy of SCR. Take the field, raise to power of 2 and
% integrate over all space.
Utot_SCR = 0.5 * epsilon * epsilon_0 * trapz(z_SCR, E_SCR.^2);
% e/(V-cm) * [(V/cm)^2]*cm = e/V * 1/cm * V^2/cm^2 * cm
% = e/V * V^2/cm^2 = eV/cm^2

% Total energy of parallel capacitor at field corresponding to Qsheet
% and with a width of z_d - sheetZ
capacitor_field = -Qsheet / epsilon_0;
% e/cm^2 / e/(V-cm) = e/cm^2 * (V-cm)/e = V/cm
capacitor_width = (z_d - sheetZ) * 1e-8;     % Ang * cm/Ang = cm
Utot_capacitor = 0.5 * epsilon_0 * ...
    capacitor_width * capacitor_field.^2;
% e/(V-cm) * cm * [(V/cm)^2] = eV/cm^2

% Correction: add SCR energy and remove capacitor energy
Utot_correction = Utot_SCR - Utot_capacitor;     % eV/cm^2

end

