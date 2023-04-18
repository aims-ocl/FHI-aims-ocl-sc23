%% Returns coefficients of line given two points on it.
% P1, P2 - points, each a 1x2 matrix.
function coeffs = getLineFromPAndSlope(P, slope)
% Unpack for more convenient reading
x0 = P(1);
y0 = P(2);

a1 = slope;
a0 = y0 - slope * x0;
% Pack for output
coeffs = struct;
coeffs.a0 = a0;
coeffs.a1 = a1;
% 
% % Plot result
% x = (0:0.1:11)';
% plot(x, a.*x + b);
% plot([P1(1) P2(1)]', [P1(2) P2(2)]', 'LineStyle', 'none', 'Marker', '^')
end