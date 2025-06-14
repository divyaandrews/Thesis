function y=nogelik(z)
global x n
lam=zeros(1,n);
lam(1)=mean(x);
for i=2:n
    lam(i)=z(1)+z(2)*x(i-1)+z(3)*lam(i-1);
        % Calculate the indf function (indicator function)
        indicator = ind(x(i));
        % Compute the likelihood components
        f(i) = indicator * log(z(4)) + ...
               (1 - indicator) * (log(1 - z(4)) + ...
               (x(i) - 1) * log(1 - (1 / lam(i))) - log(lam(i)));
    end
    
    % Sum the likelihood components (ignoring NaN values)
    y = sum(-f, 'omitnan');
end
%y=-log(la(2:n))*x(2:n)'+sum(la(2:n))+sum(log(factorial(x(2:n))));