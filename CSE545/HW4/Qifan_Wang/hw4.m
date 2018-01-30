%% part c

alphaz=25;
betaz=6;
alphax=8;
y0=0;
x0=1;
z0=0;
N=10;
g=1;
c =[1.0000 0.6294 0.3962 0.2494 0.1569 0.0988 0.0622 0.0391 0.0246 0.0155];
sigma2=[41.6667 16.3934 6.5359 2.5840 1.0235 0.4054 0.1606 0.0636 0.0252 0.0252]/1000;


t=0:0.001:1;
y=zeros(size(t));
z=zeros(size(t));
ydotdot=zeros(size(t));
x=ones(size(t));




%% part c
w = [0 0 0 0 0 0 0 0 0 0];
gama_over_t=[];
for i=1:1000 % starting from t = 0.001
    % Nonlinear function
    ten=1:10;
    gama = exp(-1./(2*sigma2(ten)).*(x(i)-c(ten)));
    gama_over_t=[gama_over_t;gama];
    f=sum(gama.*w)/sum(gama)*x(i)*(g-y0);
    
    
    % Transformation System
    cur_t = t(i+1);
    zdot= alphaz*(betaz*(g-y(i))-z(i))+f;
    ydot=z(i);
    
    % Canonical System
    xdot=-alphax*x(i);
    
    %update
    x(i+1)=x(i)+xdot*0.001;
    y(i+1)=y(i)+ydot*0.001;
    z(i+1)=z(i)+zdot*0.001;
    ydotdot(i+1)=zdot;
end

fig1=figure;
for i=1:10
    plot(0:0.001:0.999,gama_over_t(:,i));
    hold on;
end
title('PSI graph over time')
print(fig1, 'c2','-dpng');

fig2=figure;
plot(t,y);hold on;
plot(t,z);hold on;
plot(t,ydotdot);hold on;
hold off
legend('y','ydot','yddot');
title('[y ydot yddot] graph over time')
print(fig2, 'c1','-dpng');

fig3=figure;
plot(t,x)
title('x graph over time')
print(fig3, 'c3','-dpng');

%% part d
w1=[0 0 0 0 1 1 1 1000 1000 100000];
w=w1;
gama_over_t=[];
for i=1:1000 % starting from t = 0.001
    % Nonlinear function
    ten=1:10;
    gama = exp(-1./(2*sigma2(ten)).*(x(i)-c(ten)));
    gama_over_t=[gama_over_t;gama];
    f=sum(gama.*w)/sum(gama)*x(i)*(g-y0);
    
    
    % Transformation System
    cur_t = t(i+1);
    zdot= alphaz*(betaz*(g-y(i))-z(i))+f;
    ydot=z(i);
    
    % Canonical System
    xdot=-alphax*x(i);
    
    %update
    x(i+1)=x(i)+xdot*0.001;
    y(i+1)=y(i)+ydot*0.001;
    z(i+1)=z(i)+zdot*0.001;
    ydotdot(i+1)=zdot;
end

fig1=figure;
for i=1:10
    plot(1:1000,gama_over_t(:,i));
    hold on;
end

fig2=figure;
plot(t,y);hold on;
plot(t,z);hold on;
plot(t,ydotdot);hold on;
hold off
legend('y','ydot','ydotdot');

%% part e
imitation = load('imitation.data');
fig3=figure;
plot(t,imitation(:,1));hold on;
plot(t,imitation(:,2));hold on;
plot(t,imitation(:,3));hold on;
hold off
legend('y','ydot','ydotdot');

tau = imitation(:,3)-alphaz*(betaz*(g-imitation(:,1))-imitation(:,2));
X=[];
for i=1:1001
    X = [X;exp(-1./(2*sigma2(ten)).*(x(i)-c(ten)))/sum(exp(-1./(2*sigma2(ten)).*(x(i)-c(ten))))*x(i)*(g-y0)];
end

imitated_w = [];
imitated_w = X\tau;

residual = sum(tau-X*imitated_w);
disp('residual');
disp(residual)

w=imitated_w';

gama_over_t=[];
for i=1:1000 % starting from t = 0.001
    % Nonlinear function
    ten=1:10;
    gama = exp(-1./(2*sigma2(ten)).*(x(i)-c(ten)));
    gama_over_t=[gama_over_t;gama];
    f=sum(gama.*w)/sum(gama)*x(i)*(g-y0);
    
    
    % Transformation System
    cur_t = t(i+1);
    zdot= alphaz*(betaz*(g-y(i))-z(i))+f;
    ydot=z(i);
    
    % Canonical System
    xdot=-alphax*x(i);
    
    %update
    x(i+1)=x(i)+xdot*0.001;
    y(i+1)=y(i)+ydot*0.001;
    z(i+1)=z(i)+zdot*0.001;
    ydotdot(i+1)=zdot;
end

fig1=figure;
for i=1:10
    plot(0:0.001:0.999,gama_over_t(:,i));
    hold on;
end
title('PSI graph over time')
print(fig1, 'e2','-dpng');

fig2=figure;
plot(t,y);hold on;
plot(t,z);hold on;
plot(t,ydotdot);hold on;
hold off
legend('y','ydot','yddot');
title('[y ydot yddot] graph over time')
print(fig2, 'e1','-dpng');

fig3=figure;
plot(t,x)
title('x graph over time')
print(fig3, 'e3','-dpng');
%%
figure;
plot(t,y);hold on;
plot(t,imitation(:,1));hold on;
hold off
legend('y','imitation y');
title('Comparing y graph over time')
print( 'compare_y','-dpng');

figure;
plot(t,z);hold on;
plot(t,imitation(:,2));hold on;
hold off
legend('ydot','imitation ydot');
title('Comparing ydot graph over time')
print( 'compare_ydot','-dpng');

figure;
plot(t,ydotdot);hold on;
plot(t,imitation(:,3));hold on;
hold off
legend('yddot','imitation yddot');
title('Comparing yddot graph over time')
print('compare_yddot','-dpng');