% Function f (exercise 2)

function f = fun2(x)

    A = [ 6, 2, 1; 2, 5, 2; 1, 2, 4];
    b = [8; 3; 3];
    f = 0.5*transpose(x)*A*x - transpose(b)*x;
    % add terms here
end