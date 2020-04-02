pd = makedist('Gamma',0.001,0.001);
x = 0.1:20
y = pdf(pd, x);
figure
plot(x,y)
grid