[b,a]=butter(2,0.1);

noisy=load('noisy.data');
y_filtered=filter(b,a,noisy(:,1));
hold on
plot(y_filtered);
plot(noisy(:,2));
hold off

legend('x_f','x_f true')
%%
A=0.5;
B=3.5;
C=1;
t=0:0.01:9.99;
Q=0.01;
R=1;
%%
yv=noisy(:,1);
u=noisy(:,3);
%%
P = B*Q*B';         % Initial error covariance
x = zeros(1,1);     % Initial condition on the state
ye = zeros(length(t),1);
ycov = zeros(length(t),1);
P_time=zeros(length(t),1);
K_time=zeros(length(t),1);

%%
for i = 1:length(t)
  % Measurement update
  Mn = P*C'/(C*P*C'+R);
  x = x + Mn*(yv(i)-C*x);   % x[n|n]
  P = (eye(1)-Mn*C)*P;      % P[n|n]
  ye(i) = x;
   P_time(i)=P;
   K_time(i)=P*C'*inv(C*P*C'+R);
  x = A*x + B*u(i);        % x[n+1|n]
  P = A*P*A' + Q;     % P[n+1|n]
end

%%
figure
plot(ye)
hold on
plot(noisy(:,2))
hold off
legend('Estimated X','True X')

%%
figure
plot(P_time)
hold on
plot(K_time)
hold off
legend('gain K','P')