M=10; h=2^M; N=10000;
K=0.6; b=1.4; T=1;

flag = zeros(N,1);
for i = 1:N
    delta = normrnd(0,sqrt(T/h),h,1);
    W = cumsum(delta);
    flag(i) = (W(h)>=K)&&(max(W)<b);
end
mean1 = sum(flag)/N
std1 = std(flag)/sqrt(N)
lower1 = mean1+norminv(0.025,0,1)*std1
upper1 = mean1+norminv(0.975,0,1)*std1

Nsteps = zeros(N,1);
flag = zeros(N,1);
for i = 1:N
    [Nsteps(i),flag(i)] = BMBridge(M,T,b,K);
end
mean2 = mean(flag)
std2 = std(flag) /sqrt(N)
lower2 = mean_BMB+norminv(0.025,0,1)*std2
upper2 = mean_BMB+norminv(0.975,0,1)*std2
mean_steps = mean(Nsteps)



