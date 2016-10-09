clear;
clc
clf
A2 =load('solution2.dat');
A0 =load('solution0.dat');
A1 =load('solution1.dat');
A3 =load('solution3.dat');


dx0=A0(:,2);
dy0=A0(:,3);
dx1=A1(:,2);
dy1=A1(:,3);
dx2=A2(:,2);
dy2=A2(:,3);
dx3=A3(:,2);
dy3=A3(:,3);


figure(1)
scatter(dx0,dy0,'filled')
grid on
title('1000 particles, core 0')

figure(2)
scatter(dx1,dy1,'filled')
grid on
title('1000 particles, core 1')

figure(3)
scatter(dx2,dy2,'filled')
grid on
title('1000 particles, core 2')

figure(4)
scatter(dx3,dy3,'filled')
grid on
title('1000 particles core3')


