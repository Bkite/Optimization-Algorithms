function [x,y,obj] = SimulatedAnnealing

x0 = 0;
y0 = 1;

    function f_xy = Q3f(x,y)
        f_xy = 3*(1-x)^2*exp(-x^2 - (y+1)^2) - 10*(x/5 - x^3 - y^5)*exp(-x^2 - y^2) - exp(-(x+1)^2 - y^2)/3;
    end

    function p = Q3p(T_k,z,x)
        p = min(1,exp(-(z - x)/T_k));
    end


x = x0;
y = y0;
k = 0;

while k < 200
    
    T_k = 1000/(10^k);
    k = k+1;
    x_old = x;
    y_old = y;
    
    delta_x = -1+2*rand(1);
    delta_y = -1+2*rand(1);

    z_x = x + delta_x;
    z_y = y + delta_y;
    r = rand(1);
    
    Q3f(z_x,z_y);
    Q3f(x,y);

    if r <= Q3p(T_k,Q3f(z_x,z_y),Q3f(x,y))
        x = z_x;
        y = z_y;
    end

    if abs(Q3f(x,y)-Q3f(x_old,y_old)) <0.01 && x ~= x_old && y~=y_old
        break
    end

end

obj = Q3f(x,y);

f_xy = zeros(601,601);
[plotx,ploty] = meshgrid(-3:0.01:3);
for i = 1:601
    for j = 1:601
        f_xy(i,j) = Q3f(plotx(i,j),ploty(i,j));
    end
end
[C,h] = contour(plotx,ploty,f_xy,-7:1:8);
clabel(C,h)
zlabel('f(x, y)')
xlabel('x')
ylabel('y')
grid on
hold on

end