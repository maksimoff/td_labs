function  vle_2_3
%VLE_2_3 Summary of this function goes here
%   Detailed explanation goes here
    P_s_water = 149.9*1.3332E-3;
    P_s_propanol= 152.6*1.3332E-3;
    T = 273.15+60;
    R=8.31;
    X_exp=[0.983; 0.948; 0.910; 0.835; 0.765; 0.707; 0.590; 0.393; 0.249; ...
        0.122; 0.063; 0.057; 0.020];
    P_propanol_exp=[163.8;
172.6;
175.5;
172.4;
167.8;
165.4;
162.5;
159.1;
156.5;
139.2;
108.5;
93.3;
58.0;];



    a_propanol = @(g00, x) exp(x.^2.*g00./R./T).*x; 
    a_water= @(g00, x) exp((1-x).^2.*g00./R./T).*(1-x);
    P_propanol= @ (g00, x) P_s_propanol.*a_propanol( g00, x); 
    P_water= @(g00, x) P_s_water.*a_water(g00,x);
    y = @(g00, x) P_propanol(g00, x)./(P_propanol(g00,x)+P_water(g00,x)); 
    P= @(g00, x) P_propanol(g00,x)+P_water(g00,x);
    Fun = @(g00) P_propanol(g00, X_exp)-P_propanol_exp;
    P_exp = @(g00, x) P_propanol_exp+P_water(g00, x);
    g00_calc = lsqnonlin(Fun, 1);
    t=[0:0.01:1];
    plot(t.', P(g00_calc, t.'),'-b',  y(g00_calc, t.'), P(g00_calc, t.'),'-r',  X_exp, P_exp(g00_calc, X_exp), '+k'); 
end

