function F = ftn_fit_FRAP_2cpt(P,t,ycyt,c0,Cnucbar,Ccytbar,KnucG)
%
%
% 
%
% This function is to fit the two-component Cact-LT FRAP recovery.

%
% Unpack parameters
%
kout = P(1);
koutG = P(2);
phi = P(3);

%
% Parameter relationships
%
phicyt = 1 - (1 - phi)*Cnucbar/Ccytbar/KnucG;
Knuc = phi*Cnucbar/phicyt/Ccytbar;
kin = Knuc*kout;
kinG = KnucG*koutG;

%
% Create param vectors for the two scenarios
%
pG = [(1-phi)*c0 kinG koutG];
p = [phi*c0 kin kout];



%
% Calculating the two recovery curves
%
GFP = ftn_FRAP(pG,t,ycyt*(1-phicyt));
CLG = ftn_FRAP(p,t,ycyt*phicyt);
F = CLG + GFP;