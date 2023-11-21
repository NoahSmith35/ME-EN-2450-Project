function [t,y] = rk4(odefun, tspan, y0)
n = length(tspan);
y = zeros(1,n);
y(1) = y0;

for i=2:n
    h = tspan(i) - tspan(i-1);
    k1 = odefun(tspan(i-1),y(i-1));
    k2 = odefun(tspan(i-1)+0.5*h,y(i-1)+0.5*k1*h);
    k3 = odefun(tspan(i-1)+0.5*h,y(i-1)+0.5*k2*h);
    k4 = odefun(tspan(i-1)+h,y(i-1)+k3*h);
    
    phi = (k1+2*k2+2*k3+k4)/6;

    y(i) = y(i-1) + phi*h;
  
end
t = tspan;
end
