function [c,ceq]=myfun21(x)
global cov_matrix;
global average_row;
global budget;
global no_stocks;

c= [];
ceq=0;
for i=1:no_stocks
ceq=ceq+(((0.001)*x(i)*x(25))/(average_row(i)));
end
ceq=ceq+x(no_stocks+1)-budget; %Transaction cost constraint.
end