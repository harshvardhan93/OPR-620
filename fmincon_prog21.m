close all;
clear all;
% Construct a dialog box with three options
global budget;
load('60_workspace.mat');
choice = questdlg('Please select the risk level', ...
	'Risk Level', ...
	'Low Risk','Optimal Risk','High Risk','Optimal Risk');
switch choice
    case 'Low Risk'
        budget = 3000;
    case 'Optimal Risk'
        budget = 5000;
    case 'High Risk'
        budget = 7000;
end

ADJRatesAll=OPR620DatafinalS2;
global cov_matrix;
global average_row;
global no_stocks;

ADJRatesAll=ADJRatesAll*100; %Converting all the holding period rates to percentages
no_stocks=length(ADJRatesAll); 
average_row=zeros(1,no_stocks);
xx=ones(1,no_stocks);
x_initial=xx/(no_stocks); %setting the initial point for the algorithm
x_initial(no_stocks+1)=0;

cov_matrix=cov(ADJRatesAll); %calculating the covariance matrix which will be used to calculate the risk
var_matrix=diag(cov_matrix); 
average_row=mean(ADJRatesAll); %Average return of every stock
options = optimset('Algorithm','interior-point','Display', 'off'); %Setting up the algorithm
r = 1.48; % r denotes the risk free rate
Aeq=ones(1,no_stocks); %Setting up the equality constraint that all the weights should sum upto 1
Aeq(1,no_stocks+1)=0;
Beq=1;
lb=zeros(1,no_stocks+1); %Setting up lower bound equal to zero to avoid short selling
ub=ones(1,no_stocks+1); %Setting up the upper bound to the maximum value of sum of all the weights in the portfolio
ub(1,no_stocks+1)=budget; %Setting the upper bound for the transaction cost
average_row(1,no_stocks+1)=0;
cov_matrix=[cov_matrix zeros(no_stocks,1)];
cov_matrix=[cov_matrix; zeros(1,no_stocks+1)];
if no_stocks<=30
    param=0.01; %Setting the initial value for beta
    iterations=100;
    else 
        param=0.1;
        iterations=200;
    end
weights=zeros(no_stocks+1,iterations); %Setting up the weights matrix to hold the value of weights for various values of beta and calculating the corresponding values of risk and return
    
    fval=0; 
if iterations==100
    for i =1:100 
    [weights(:,i),fval(i),exitflag,output]=fmincon(@(x)((-x*average_row')+(param*(x*cov_matrix*x'))),x_initial,[],[],Aeq,Beq,lb,ub,@myfun21,options);
    i_weight=weights(:,i);
    return1(i)=i_weight'*average_row';  %Calculating the value of return & risk based on the weights obtained by running fmincon
    risk(i)=i_weight'*(cov_matrix*i_weight);
    psharpe(i)=((return1(i) - r)/sqrt(risk(i))); %Finding the sharpe's ratio
    param=param+0.005; 
    end
else   
    for i =1:200 
    [weights(:,i),fval(i),exitflag,output]=fmincon(@(x)((-x*average_row')+(param*(x*cov_matrix*x'))),x_initial,[],[],Aeq,Beq,lb,ub,@myfun21,options);
    i_weight=weights(:,i);
    return1(i)=i_weight'*average_row';  %Calculating the value of return & risk based on the weights obtained by running fmincon
    risk(i)=i_weight'*(cov_matrix*i_weight);
    psharpe(i)=((return1(i) - r)/sqrt(risk(i))); %Finding the sharpe's ratio
    param= param*1.2; 
    end
end
[C,I] = max(psharpe); %Calculating the maximum sharpe's ratio out of the values obtained in the matrix.
average_row(no_stocks+1)=[];
hold on;

plot(sqrt(risk),return1); %plotting the efficient frontier 
hold on;
title('Tangency Portfolio')
xlabel('Standard Deviation')
ylabel('Return');
plot(sqrt(risk),psharpe);
x = 0:1:7;
y = C*x+r;
line(x,y,'Color', 'r', 'LineWidth', 1); %Plotting the Capital Market Line
total_return=((10000-budget)*(100+r) + budget*(100+(-fval(I))))/100 
message = sprintf('Return at %s is $%f',choice,total_return);
msgbox(message) %Displaying the final value of return
hold off
figure
hold on
plot(risk,return1)
xlabel('Risk')
ylabel('Return')
scatter(var_matrix,average_row); %Plotting the Efficient frontier with the stock representation. 
hold off;