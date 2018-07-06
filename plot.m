clc
p=linspace(0,1,10000);  %probability vector
[DoFMax,H,msg]=montecarlo(6000);

A1=2*p+(1-(1-p).^2+p.*(1-p).^3).*(1+(1-p).^2);
B = 3+(1+(1-p).^3).*(1-(1-p).^2+p.*(1-p).^3)+p.*(1+(1-p).^2);



DoF_Theorem_5 = 1/3*(1-p).*(1+(1-p).^3+B.*p); % analytical expression of Theorem 5

DoF_Theorem_4 = 4/5*(1-p)...
    + 2/5 * p.^2 .* (1-p)...
    + 1/5 * (1-p) .* (1-(1-p).*(1-p.*(1-p))).*(1-(1-p).^2.*(1-p.*(1-p))); % analytical expression of Theorem 4

figure;
plot(p, DoF_Theorem_4, 'green','LineWidth',2);hold all;
plot(p, DoF_Theorem_5, 'blue','LineWidth',2);
plot(p, DoFMax, 'red','LineWidth',2);
legend(...'1/6', '1/5', '3/10', '1/7', ...'numeric 1/6',...'analytic 2/6',  ... 'analytic 1/6 new', ...'analytic 4/5', ...'analytic 1/5', ...
     'K = 5, f(p) = 3/5', 'K \rightarrow \infty, f(p) = 0', 'maximum DoF'); xlabel('p'); ylabel('puDoF');
