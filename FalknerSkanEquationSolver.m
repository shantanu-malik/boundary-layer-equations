function FalknerSkanEquationSolver()
% This function solves and plots the Falkner-Skan equation using shooting method
% and 4th order Runge-Kutta method for m values ranging from 0 to 2
% Falkner-Skan equation --> [f''' + 0.5*(m+1)f*f'' + m*(1 - f'^2) = 0]

l1 = 0;
h = 0.1;
eta = 0:h:10-h;
f = zeros(10/h,3);
f_dPrime = zeros(1,21);
figure('name','Solution of Falkner-Skan Equation (for m = 0:0.1:2)','NumberTitle','off');
j=0;

for m = 0:0.1:2
    l2 = l1 + 0.5 -m/4;
    if m > 1.7
        l2 = l1 + 0.5 - m/4.5
    end
    %l2 = l1 + 0.5;
    while abs(f(end,2)-1)>0.00001
        shot = 0.5*(l1 + l2);
        f(1,:) = [0 0 shot];
        
        for i = 1:(10/h)-1
            F = f(i,:);
            k1 = h*([F(2) F(3) (-0.5*(m+1)*F(1)*F(3)-m*(1-(F(2))^2))]);
            F = f(i,:) + k1/2;
            k2 = h*([F(2) F(3) (-0.5*(m+1)*F(1)*F(3)-m*(1-(F(2))^2))]);
            F = f(i,:) + k2/2;
            k3 = h*([F(2) F(3) (-0.5*(m+1)*F(1)*F(3)-m*(1-(F(2))^2))]);
            F = f(i,:) + k3;
            k4 = h*([F(2) F(3) (-0.5*(m+1)*F(1)*F(3)-m*(1-(F(2))^2))]);
            f(i+1,:) = f(i,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
        end
        
        if(f(end,2)<1)
            l1 = shot;
        else
            l2 = shot;
        end
    end
    j = j + 1;
    f_dPrime(j) = f(1,3);
    
    subplot(1,6,1);
    plot(f(:,1),eta);hold on;
    xlabel('f');ylabel('\eta');
    subplot(1,6,2);
    plot(f(:,2),eta);hold on;
    xlabel('f ''');ylabel('\eta');
    subplot(1,6,3);
    plot(f(:,3),eta);hold on;
    xlabel('f ''''');ylabel('\eta');
    subplot(1,6,4);
    plot(-0.5*f(:,1).*f(:,3),eta);hold on;
    xlabel('f ''''''');ylabel('\eta');
    
    f = zeros(10/h,3);
end

subplot(1,6,6);
plot(0:0.1:2,f_dPrime,'*');
xlabel('m');ylabel('f ''''(0)');
end