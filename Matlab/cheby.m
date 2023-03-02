syms x y
ezplot(chebyshevT(8,x))
axis([-1.5 1.5 -2 2])
grid on
ylabel('T_n(x)')
legend('T_0(x)','T_1(x)','T_2(x)','T_3(x)','T_4(x)','Location','Best')
title('Chebyshev polynomials of the first kind')
