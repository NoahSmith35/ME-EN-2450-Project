function [t, y] = rk4(odefun, tspan, y0)
    n = length(tspan);
    y = zeros(5, n);
    y(1,1) = y0(1);
    y(2,1) = y0(2);
    y(3,1) = y0(3);
    y(4,1) = y0(4);
    y(5,1) = y0(5);

    for i = 2:n
     for j = 1:5
        h = tspan(i) - tspan(i-1);
        k1 = odefun(tspan(i-1), y(:,i-1),j,i-1);
        k2 = odefun(tspan(i-1) + 0.5 * h, y(:,i-1) + 0.5 * k1 * h,j,i-1);
        k3 = odefun(tspan(i-1) + 0.5 * h, y(:,i-1) + 0.5 * k2 * h,j,i-1);
        k4 = odefun(tspan(i-1), y(:,i-1) + k3 * h,j,i-1);

        phi = (k1 + 2*k2 + 2*k3 + k4) / 6;

        y(j,i) = y(j,i-1) + phi * h;
        
    end

    t = tspan;
    end
end
