
function [invertedA detA]  = invertA(A)

a = A(1,1);
b = A(1,2);
sigmasq = a - b;
%sigmasq = 0.01;
n = size(A,1);

d = (1 - (a/sigmasq))/((n*b) + sigmasq);
%c = d + (1/sigmasq);

invertedA = (ones(n).*d) + (eye(n)*(1/sigmasq));
detA = det(A);
% Better way to do the above:
newsigmasq = sigmasq/b;

detA = (newsigmasq^n + (n*newsigmasq^(n-1)))*b^n;

end
