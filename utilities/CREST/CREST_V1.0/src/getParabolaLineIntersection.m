%% Finds intersections of a line and a parabola
% Returns set of 2 points representing the intersection. If one of the
% points is empty, the lines intersect at one point. If both empty, no
% intersection. Input:
% parabola_coeffs (3 x 1 matrix)
% line_coeffs (2 x 1 matrix)
function [P1, P2] = getParabolaLineIntersection(parabola_coeffs, line_coeffs)
% Threshold taken as 0
thr = 1e-12;
% Unpack for ease of reading
a2 = parabola_coeffs.a2;
a1 = parabola_coeffs.a1;
a0 = parabola_coeffs.a0;
b1 = line_coeffs.a1;
b0 = line_coeffs.a0;

% Make sure the parabola does have some curvature
if abs(a2) < thr
    % Must be a straightline. Treat accordingly
    if abs(a1 - b1) < thr
        % Parallel lines
        P1 = [];
    else
        x1 = (b0 - a0) / (a1 - b1);
        y1 = b1 * x1 + b0;
        P1 = [x1 y1];
    end
    P2 = [];
else
    % Parabola with nonzero curvature
    d = (a1 - b1)^2 - 4 * a2 * (a0 - b0);
    % Check for 1 or no solutions
    if abs(d) < thr
        % One solution
        x1 = -(a1 - b1) / (2 * a2);
        x2 = [];
    elseif d < 0
        % No solutions
        x1 = [];
        x2 = [];
    else
        % 2 solutions
        x1 = (-(a1 - b1) - sqrt(d)) / (2 * a2);
        x2 = (-(a1 - b1) + sqrt(d)) / (2 * a2);
    end
    
    if ~isempty(x1)
        y1 = b1 * x1 + b0;
        P1 = [x1 y1];
        if ~isempty(x2);
            y2 = b1 * x2 + b0;
            P2 = [x2 y2];
        else
            P2 = [];
        end
    else
        P1 = [];
        P2 = [];
    end
end

end









