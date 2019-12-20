% Function f from exercise 1
function [f] = fun1(x)
    consts = [100, 1, 90, 1, 10, 0.1]; % coefficients
    f = consts(1)*(x(2)-x(1)^2)^2;
    f = f + consts(2)*(1-x(1))^2;
    f = f + consts(3)*(x(4)-x(3)^2)^2;
    f = f + consts(4)*(1-x(3))^2;
    f = f + consts(5)*(x(2)+x(4)-2)^2;
    f = f + consts(6)*(x(2)-x(4))^2;
    % add terms here
end

