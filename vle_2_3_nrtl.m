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
    P_propanol_exp= P_propanol_exp.*1.3332E-3;
    G_ba= @(g00) exp(-g00(1)*g00(3));
    G_ab= @(g00) exp(-g00(1)*g00(2));
    opt = optimoptions('lsqnonlin', 'Display', 'iter-detailed', 'MaxFunctionEvaluations', 5000, 'MaxIterations', 5000);
    ln_Ga= @(x, g00)  x.^2.*(g00(3).*(G_ba(g00)./(1-x+x.*G_ba(g00))).^2+g00(2).*G_ab(g00)./(x+(1-x).*G_ab(g00)).^2);
    ln_Gb = @(x, g00) (1-x).^2.*(g00(2).*(G_ab(g00)./(1-x+x.*G_ab(g00))).^2+g00(3).*G_ba(g00)./(1-x+x.*G_ba(g00)).^2);
    a_propanol = @(g00, x) exp(ln_Gb(x, g00)).*x; 
    a_water= @(g00, x) exp(ln_Ga(x, g00)).*(1-x);
    P_propanol= @ (g00, x) P_s_propanol.*a_propanol( g00, x); 
    P_water= @(g00, x) P_s_water.*a_water(g00,x);
    y = @(g00, x) P_propanol(g00, x)./(P_propanol(g00,x)+P_water(g00,x)); 
    P= @(g00, x) P_propanol(g00,x)+P_water(g00,x);
    Fun = @(g00) P_propanol(g00, X_exp)-P_propanol_exp;
    P_exp = @(g00, x) P_propanol_exp+P_water(g00, x);
    [g00_calc, ~,resid,~,~,~,J] = lsqnonlin( Fun, [0.4;2;-0.8],[], [], opt);
    ci = nlparci(g00_calc,resid,'jacobian',J)
    g00_calc
    t=[0:0.01:1];
    plot(t.', P_propanol(g00_calc, t.'),'-b', X_exp, P_propanol_exp, '+k'); 
end

