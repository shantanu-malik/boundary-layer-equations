function blasiusEquationSolver()
% This function solves and plots the Blasius equation using shooting method
% and 4th order Runge-Kutta method.
% Blasius equation --> [f''' + 0.5*f*f'' = 0]

l1 = 0; l2 = 1;
h = 0.1;
eta = 0:h:10-h;
F = zeros(10/h,3);

while abs(F(end,2)-1)>0.00001
    shot = 0.5*(l1 + l2);
    F(1,:) = [0 0 shot];
    
    for i = 1:(10/h)-1
        f = F(i,:);
        k1 = h*([f(2) f(3) -0.5*f(1)*f(3)]);
        f = F(i,:) + k1/2;
        k2 = h*([f(2) f(3) -0.5*f(1)*f(3)]);
        f = F(i,:) + k2/2;
        k3 = h*([f(2) f(3) -0.5*f(1)*f(3)] + k2/2);
        f = F(i,:) + k3;
        k4 = h*([f(2) f(3) -0.5*f(1)*f(3)]);
        F(i+1,:) = F(i,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
    
    if(F(end,2)<1)
        l1 = shot;
    else
        l2 = shot;
    end
end

figure('name','Solution of Blasius Equation');
subplot(1,4,1);
plot(F(:,1),eta);
xlabel('f');
ylabel('\eta');
subplot(1,4,2);
plot(F(:,2),eta);
xlabel('f ''');
ylabel('\eta');
subplot(1,4,3);
plot(F(:,3),eta);
xlabel('f ''''');
ylabel('\eta');
subplot(1,4,4);
plot(-0.5*F(:,1).*F(:,3),eta);
xlabel('f ''''''');
ylabel('\eta');
end