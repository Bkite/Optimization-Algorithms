function [k,obj] = PSO(w,c1,c2,d)

k = 0;
x0 = rand(1,3);
y0 = rand(1,3);
v0 = rand(3,2);
x = [x0.' y0.']; % 3*2
disp(x);
v = v0; % 3*2


    function f_x = Q4f(x)
        f_x = 3*(1-x(1))^2*exp(-x(1)^2 - (x(2)+1)^2) - 10*(x(1)/5 - x(1)^3 - x(2)^5)*exp(-x(1)^2 - x(2)^2) - exp(-(x(1)+1)^2 - x(2)^2)/3;
    end

f_x = zeros(1001,1001);
[x1,x2] = meshgrid(-5:0.01:5);
for i = 1:1001
    for j = 1:1001
        f_x(i,j) = Q4f([x1(i,j) x2(i,j)]);
    end
end

[C,h] = contour(x1,x2,f_x,-7:1:10);
clabel(C,h)
xlabel('x1')
ylabel('x2')
grid on
hold on

p = x; % 3*2
f_p = [Q4f(p(1,:)) Q4f(p(2,:)) Q4f(p(3,:))];
[~,i] = sort(f_p);

g = p(i(1),:); % 1*2

delta = ones(1,100);

plot3(x(:,1),x(:,2),[Q4f(x(1,:)) Q4f(x(2,:)) Q4f(x(3,:))],'.','MarkerSize',10);
hold on

while k < 100
    
    k = k+1;
    g_old = g;
    
    
    for i = 1:d
        r = rand(1,2);
        s = rand(1,2);
        v(i,:) = w*v(i,:) + c1*r.*(p(i,:) - x(i,:)) + c2*s.*(g - x(i,:));
        x(i,:) = x(i,:) + v(i,:);
        
        if Q4f(x(i,:)) < Q4f(p(i,:))
            p(i,:) = x(i,:);
            %g = x(i,:);
        end
        
        if Q4f(x(i,:)) < Q4f(g)
            g = x(i,:);
        end
    end 
    
    delta(k) = norm(abs(Q4f(g) - Q4f(g_old)));
    
    if k > 4
        delta5 = delta(k-4:k);
        
        if norm(delta5) <= 0.00000001
            break
        end
    end

    
    plot3(g(:,1),g(:,2),Q4f(g),'.','MarkerSize',10);
    hold on
    
    disp(x);
    disp(g);
    obj = Q4f(g);
    disp(obj);
end

disp(k);


end