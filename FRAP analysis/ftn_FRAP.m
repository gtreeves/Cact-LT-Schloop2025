function y = ftn_FRAP(P,t,ycyt)
%
%
%
%
%

c0 = P(1);
k_in = P(2);
k_out = P(3);

T = t-t';
T(T < 0) = NaN;
dycyt = diff(ycyt);
dt = diff(t);
E = exp(-k_out*T);

A = dycyt./dt/k_out; % col vec
E1 = E.*(k_out*t'-1);
E1diff = E1(:,2:end) - E1(:,1:end-1);
E1diff(isnan(E1diff)) = 0;

B = ycyt(1:end-1) - dycyt./dt.*t(1:end-1); % col vec
Ediff = E(:,2:end) - E(:,1:end-1);
Ediff(isnan(Ediff)) = 0;
I = E1diff*A + Ediff*B;

y = c0*exp(-k_out*t) + k_in/k_out*I;












