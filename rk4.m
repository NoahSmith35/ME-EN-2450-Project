function [t,y] = rk4(odefun, tspan, y0)
n = length(tspan);
y = zeros(1,n);
y(1)= y0(1)
for i=1:n
    h = tspan(i+1) - tspan(i)
    k1 = odefun(y(i),i)
    k2 = odefun(i+0.5*h,y(i)+0.5*k1*h);
    k3 = odefun(i+0.5*h,y(i)+0.5*k2*h);
    k4 = odefun(i+h,y(i)+k3*h);
    
    phi = (k1+2*k2+2*k3+k4)/6;

    y(i) = y(i-1) + phi*h;
  
end
t = tspan;
end