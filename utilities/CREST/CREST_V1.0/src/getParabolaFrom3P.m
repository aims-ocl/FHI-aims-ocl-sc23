%% Returns coefficients of parabola given three points on it. Based on
% Cramer's method.
% P1, P2, P3 - three points, each a 1x2 matrix.
function coeffs = getParabolaFrom3P(P1, P2, P3)

% Build matrix
A = [P1(1)^2 P1(1) 1; P2(1)^2 P2(1) 1; P3(1)^2 P3(1) 1];
% Result vector
B = [P1(2); P2(2); P3(2)];
% Determinants
dA = det(A);
dA2 = det([B A(:, 2:3)]);
dA1 = det([A(:, 1) B A(:, 3)]);
dA0 = det([A(:, 1:2) B]);
% Find coefficients
a2 = dA2 ./ dA;
a1 = dA1 ./ dA;
a0 = dA0 ./ dA;
% Pack for output
coeffs = struct;
coeffs.a0 = a0;
coeffs.a1 = a1;
coeffs.a2 = a2;
% 
% % Plot result
% x = (0:0.1:11)';
% figure; hold on
% plot(x, a.*x.^2 + b.*x + c);
% plot([P1(1) P2(1) P3(1)]', [P1(2) P2(2) P3(2)]', 'LineStyle', 'none', 'Marker', 'O')
end