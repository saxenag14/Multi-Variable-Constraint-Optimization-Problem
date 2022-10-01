% PROJECT PHASE 3

fprintf("Please uncomment the question that need to be executed and the corresponding constraints\n");

% User Inputs
ques = input("Please enter 1 for question 1 and 3, 2 for question 2\n");
n = input("Enter the number of variables\n");
fprintf('Enter the initial approximation\n');
x0 = [];                               % Taking initial approximation
for i = 1:n
   x0(i) = input("\n");
end
global R                               % Declaring R globally
global feval                           % Declaring feval globally to count number of function evaluation
R = 0.1;
feval = 0;
c = 10;
iter = 1;                              % Iter is used to count number of iterations
e1 = 0.0001;
xi = powell(x0,n);                     % Finding first point through Powell's Conjugate Direction Method
while 1
    P1 = Obj_fun(xi);                  % Calculating value of Penalty function at x1
    R = c*R;                         % Incrementing R
    xi = powell(x0,n);                 % Finding next point through Powell's Conjugate Direction Method
    P2 = Obj_fun(xi);                  % Finding Penalty Function Value at this point
    err = abs(P2-P1); 
    % For Plotting the results after each iteration
%     P3 = P2 - (R*constraint(xi));
%     fprintf('%f\n',xi);
%     fprintf('Iteration %d\n',iter);
%     fprintf('Optimum Function Value %f\n',P3);
%     fprintf('Function evaluation %d\n',feval);
    if err<e1                          % Termination Condition
        break;
    elseif iter>100                    % Maximum iteration value
        break;
    end
    iter = iter + 1;
end
if ques == 2                           % Changing the function value to its negative as 2nd question is maximization type
    P3 = -(P2 - (R*constraint(xi)));
else                                   % Getting optimum function value from Penalty Function Value
    P3 = P2 - (R*constraint(xi));
end

% Getting Output Values

fprintf("The minimum point is\n");
fprintf("%f\n",xi);
fprintf("The optimum function value is \n");
fprintf("%f\n",P3);
fprintf("Total number of iterations are \n");
fprintf("%d\n",iter);
fprintf("Total number of function evaluations are \n");
fprintf("%d\n",feval);

%% Objective Function
% Please uncomment the question that needs to be executed

function f = Obj_fun(x)
    global R
    global feval
    feval = feval + 1;
    %Problem 1
    f = (x(1) - 10)^3 + (x(2) - 20)^3;
    f = f + (R*constraint(x));

    %Problem 2
%     a1 = (((sin(2*pi*x(1)))^3)*(sin(2*pi*x(2))));
%     a2 = (x(1)^3)*(x(1) + x(2));
%     f = a1/a2;
%     f = -f + (R*constraint(x));

    %Problem 3
%      f = x(1) + x(2) + x(3);
%      f = f + (R*constraint(x));
end

%% Constraints
% Please uncomment the constraint corresponding to the question that needs to be executed

function constr = constraint(x)
    %Problem 1
    g(1) = min((x(1)-5)^2 + (x(2)-5)^2 - 100,0);
    g(2) = min(82.81 - ((x(1)-6)^2 + (x(2)-5)^2),0);
    g(3) = min((x(1)-13),0);
    g(4) = min((20-x(1)),0);
    g(5) = min(x(2),0);
    g(6) = min(4-x(2),0);
    constr = sum(g.^2);

    %Problem 2
%     g(1) = min(-(x(1)^2 - x(2) + 1),0);
%     g(2) = min(-(1 - x(1) + (x(2)-4)^2),0);
%     g(3) = min(10-x(1),0);
%     g(4) = min(10-x(2),0);
%     g(5) = min(x(1),0);
%     g(6) = min(x(2),0);
%     constr = sum(g.^2);

    %Problem 3
%     g(1) = min(1-0.0025*(x(4) + x(6)),0);
%     g(2) = min(1-0.025*(-x(4) + x(5) + x(7)),0);
%     g(3) = min(1-0.01*(-x(6) + x(8)),0);
%     g(4) = min(-(100*x(1) + (-x(6)*x(1)) + (833.33252*x(4)) + (-83333.333)),0);
%     g(5) = min(-(x(2)*x(4) + (-x(2)*x(7)) + (-1250*x(4)) + (1250*x(5))),0);
%     g(6) = min(-(x(3)*x(5) + (-x(3)*x(8)) + (-2500*x(5)) + 1250000),0);
%     g(7) = min(10000-x(1),0);
%     g(8) = min(10000-x(2),0);
%     g(9) = min(10000-x(3),0);
%     g(10) = min(x(1)-100,0);
%     g(11) = min(x(2)-1000,0);
%     g(12) = min(x(3)-1000,0);
%     g(13) = min(1000-x(4),0);
%     g(14) = min(1000-x(5),0);
%     g(15) = min(1000-x(6),0);
%     g(16) = min(1000-x(7),0);
%     g(17) = min(1000-x(8),0);
%     g(18) = min(x(4)-10,0);
%     g(19) = min(x(5)-10,0);
%     g(20) = min(x(6)-10,0);
%     g(21) = min(x(7)-10,0);
%     g(22) = min(x(8)-10,0);
%     constr = sum(g.^2);
end

%% Powell's Conjugate Direction Method

function min_point = powell(xn,n)
en = 0.001; p = 1;              % Defining accuracy and no of iterations
xn = xn';
while 1
    search1 = eye(n,n);         
    xn1 = zeros(n,n);
    for j=1:n
        alpha = Bounding_Search(xn,search1,j,n);       % Getting the alpha value
        for i=1:n
            xn1(:,j) = xn + (alpha*search1(j,:)');     % Getting new points
        end
        xn = xn1(:,j);
    end
    t = Bounding_Search(xn,search1,1,n);               % Calling Bounding Search Function
    x1 = xn + t*search1(:,1);
    A = x1;
    B = xn(:,1);
    d = A-B;
    s = d'*d;                   % Finding norm value
    m = sqrt(s);
    search2 = d/m;
    l = norm(search1);
    LI = acosd((search1*search2)/l);        %Checking linear independency
    if m<en                 % Terminating condition
        break;
    elseif LI==0            % Linear Independency condition
        break;
    elseif p>5000           % Maximum iteration value
        break;
    else
        for k = n:-1:2
            search1(:,k) = search1(:,k-1);
        end
        search1(:,1)=search2;
        r = Bounding_Search(xn,search1,1,n);
        min_point = xn + r*search1(:,1);            % Completion of first iteration
        p = p + 1;
    end
 M = p - 1;
%  fprintf('Iteration is %d\n',M);
%  fprintf('Minimum Value is %f\n',func_val(min_point,n));
end

% fprintf('Total Number of Iterations are %f\n',p);
% fprintf('Minimum point obtained by Powells Conjugate Direction Method is \n');
% fprintf('%f\n',min_point);
% fprintf('Minimum Value of the function is %f\n',func_val(min_point,n));

end

%% Bounding Search Method code for getting initial bracketing minima

function X = Bounding_Search(x,search,j,n)
    x0 = 0.25;
    delta = 0.04;
    feval = 0;
    x2 = x0; 
    x1 = x2 - abs(delta); 
    x3 = x2 + abs(delta);
    f1 = Single_Objective_Fun(x,search,x1,j,n);
    f2 = Single_Objective_Fun(x,search,x2,j,n);
    f3 = Single_Objective_Fun(x,search,x3,j,n);
    i = 1;
    feval= feval+3;
    if(f1<f2 && f2<f3)
        delta = -abs(delta);
    else
        delta = abs(delta);
        x1 = x3;
        f1 = f3;
    end
    while (f1<f2)
        f2 = f1;
        x3 = x2;
        x2 = x1;
        x1 = x2 + (2^i)*delta;
        f1 = Single_Objective_Fun(x,search,x1,j,n);
        %fprintf('The iteration number is:%d\n', i);
        %fprintf('The minimum value is: (%f,%f,%f)\n', f1,f2,f3);
        i = i+1;
        feval = feval + 1;
        if(f1>f2)
            break;
        end
    end
%     fprintf('The no of iterations using bounding phase method is: %d\n', i);
%     fprintf('The minimum point using bounding phase method lies between (%8.3f, %8.3f)', x1, x3);
%     fprintf('\nTotal number of function evaluation using bounding phase method is: %d\n', feval);
    X = GoldenSectionSearchMethod(x1,x3,x,search,j,n);
end

%% Golden Section Search Method

function min_alpha = GoldenSectionSearchMethod(x1,x3,x,search,j,n)
     a = x1;
     b = x3;
     enn = 0.001;
     aw = 0;
     bw = 1;
     Lw = 1;
     t = 1;
     feval = 1;
     while Lw > enn
         feval = feval + 1;
         w1 = aw + (0.618*Lw);
         w2 = bw - (0.618*Lw);
         a1 = (w1*(b - a)) + a;
         a2 = (w2*(b - a)) + a;
         y1 = Single_Objective_Fun(x,search,a1,j,n);
         y2 = Single_Objective_Fun(x,search,a2,j,n);
         if y1 > y2
             bw = w1;
             Lw = abs(bw - aw);
         elseif y2 > y1
             aw = w2;
             Lw = abs(bw - aw);
         else
             bw = w1;
             aw = w2;
             Lw = abs(bw - aw);
         end
         %fprintf('The iteration number is: %d\n',t);
         
         %fprintf('The minimum value is: %f\n',z)
         
     end  
     min_alpha = (a1+a2)/2;
 end
 
 function f = Single_Objective_Fun(x,search,a0,j,n)
    x = x + a0*search(j,:)';
     f = Obj_fun(x);
 end


