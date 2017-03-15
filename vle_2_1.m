function vle_2_1
    f= @(x) (1-x(1))^2+100*(x(2)-x(1)^2)^2;
    x0=[10, 10];
    fminsearch( f, x0.')
end
