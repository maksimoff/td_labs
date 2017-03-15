function vle_2_2
    R = 8.314;
    t =[0:0.01:1];
    g00 = 3350;
    X0=[0.5; 0.6];
    p_chexane_s= @(T) 10^(3.17125-780.637./(T-107.29));
    p_methethket_s= @(T) 10^(3.9894-1150.2079./(T-63.904));
    P= 2;
    G_g = @(T, y) R*T*log(P)+(1-y)*R*T*log(1-y)+y*R*T*log(y);
    G_l = @ (T, x) (1-x)*R*T*log(p_chexane_s(T))+x*R*T*log(p_methethket_s(T)) ...
        +R*T*((1-x)*log(1-x)+x*log(x))+x*(1-x)*g00;
    G = @(T, X, n) X(2).*G_g(T, (1-n-X(1))/X(2)+X(1))+(1-X(2)).*G_l(T, X(1));
    G_plot = @(n) fmincon(@(X) G(350, X, n), X0 );
    G_plot(0.2)
    
    
end
