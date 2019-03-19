function [p_s,p_nl,p_l,Q2f_s,Q2f_nl,Q2f_l] = NMSimplex(x0)

    function f = Q2f(x)

    f = (x(2) - x(1))^4 + 12*x(1)*x(2) - x(1) + x(2) -3;

    end

plot_fx = zeros(231,231);
[x1,x2] = meshgrid(-0.9:0.01:1.4);
for i = 1:231
    for j = 1:231
        plot_fx(i,j) = (x2(i,j)-x1(i,j))^4 + 12*x1(i,j)*x2(i,j) - x1(i,j) + x2(i,j) - 3;
    end
end

[C,h] = contour(x1,x2,plot_fx,-5:1:5);
clabel(C,h)
xlabel('x1')
ylabel('x2')
hold on

lamda = 0.2;

p0 = x0;
p1 = p0 + lamda*[1;0];
p2 = p0 + lamda*[0;1];

p = [p0 p1 p2];
f_p = [Q2f(p0) Q2f(p1) Q2f(p2)];
[~,i] = sort(f_p);

p_l = p(:,i(3));
p_nl = p(:,i(2));
p_s = p(:,i(1));
fill3([p_s(1) p_nl(1) p_l(1)],[p_s(2) p_nl(2) p_l(2)],[Q2f(p_s) Q2f(p_nl) Q2f(p_l)],'white');
hold on
k = 0;

while k < 100 && polyarea([p_l(1) p_nl(1) p_s(1)],[p_l(2) p_nl(2) p_s(2)])> 0.000001
    
    k = k + 1;
    p_g = (p_s + p_nl)/2;

    p_r = p_g + 1*(p_g - p_l);

    if Q2f(p_r) < Q2f(p_nl) && Q2f(p_r) > Q2f(p_s)
        p_l = p_r;

    elseif Q2f(p_r) < Q2f(p_s)
        p_e = p_g + 2*(p_r - p_g);
    
        if Q2f(p_e) < Q2f(p_r)
            p_l  = p_e;
        else
            p_l = p_r;
        end
    
    elseif Q2f(p_r) < Q2f(p_l) && Q2f(p_r) > Q2f(p_nl)
        p_c = p_g + 0.5*(p_r - p_g);
    
        if Q2f(p_c) < Q2f(p_l)
            p_l = p_c;
        else
            p_l = p_s + 0.5*(p_l - p_s);
            p_nl = p_s + 0.5*(p_nl - p_s);
        end
    
    elseif Q2f(p_r) > Q2f(p_l)
        p_c = p_g + 0.5*(p_l - p_g);
    
        if Q2f(p_c) < Q2f(p_l)
            p_l = p_c;
        else
            p_l = p_s + 0.5*(p_l - p_s);
            p_nl = p_s + 0.5*(p_nl - p_s);
        end 

    end
    Q2f_s = Q2f(p_s);
    Q2f_nl = Q2f(p_nl);
    Q2f_l = Q2f(p_l);
    fill3([p_s(1) p_nl(1) p_l(1)],[p_s(2) p_nl(2) p_l(2)],[Q2f_s Q2f_nl Q2f_l],'yellow');
    grid on
    hold on
    
    p = [p_l p_nl p_s];
    f_p = [Q2f(p_l) Q2f(p_nl) Q2f(p_s)];
    [~,i] = sort(f_p);

    p_l = p(:,i(3));
    p_nl = p(:,i(2));
    p_s = p(:,i(1));
    
    disp('Simplex in this iteration ');
    disp(p_l);
    disp(p_nl);
    disp(p_s);
    
end

end