%% Part2
random_bits = round(rand(1,5))
rando_bit = ones(1,5);
rando_bit(random_bits == 0) = -1;
T = 0.1;
Fs = 1000;
t = 0:1/Fs:T-1/Fs;
tnew = 0:1/Fs:5*T-1/Fs;

A1 = zeros(1,100);
% A0 = zeros(1,100);
% sig = zeros(3,100);

A1(t<T/4) = t(t<T/4);
A1(t>=T/4 & t<T*3/4) = T/2 - t(t>=T/4 & t<T*3/4);
A1(t>=T*3/4 & t<=T) = t(t>=T*3/4 & t<=T) - T;

% A0(t<T/4) = t(t<T/4);
% A0(t>=T/4 & t<T*2/4) = T/2 - t(t>=T/4 & t<T*2/4);
% A0(t>=T*2/4 & t<T*3/4) = t(t>=T*2/4 & t<T*3/4) - T/2;
% A0(t>=T*3/4) = T -t(t>=T*3/4);
x1 = zeros(1,500);
x2 = zeros(1,500);

% s = zeros(5, 500);
for i = 1:5
    x1(((i-1)*100)+1:(i*100)) = rando_bit(i)*A1(1:100);
%     x2(((i-1)*100)+1:(i*100)) = rando_bit(i)*A0(1:100);
%     s(i,((i-1)*100)+1:(i*100)) = rando_bit(i)*A1(1:100);
end
figure();
plot(tnew, x1);
xlabel('t');
ylabel('x1(t)');

% figure();
% plot(tnew, x2);
% xlabel('t');
% ylabel('x2(t)');

E1 = sum(A1.^2);
% E2 = sum(x2.^2)/0.5;

Y1 = A1/(sqrt(E1));
figure();
plot(t, Y1);
xlabel('t'); ylabel('Y1');
title('Orthonormal Signal');

A11 = [sum(A1.*Y1) -sum(A1.*Y1)] ;

% d2 = x2 - x21.*Y1;

% Y2 = d2/(sqrt(sum(d2.^2)));
figure();
plot(A11,[0 0], '*b');
xlabel("Y1");
title('Orthonormal Signal Space');

%% Snr
figure()
plot(tnew,x1);
hold on;
rx1 = x1 + (10^(-4))*(randn(1,500));
plot(tnew, rx1);
legend('Original', 'Variance 10e-4');
xlabel('t');
ylabel('signals');
figure();
rx2 = x1 + (10^(-2))*(randn(1,500));
plot(tnew,x1);
hold on;
plot(tnew, rx2);
legend('Original', 'Variance 10e-2');
xlabel('t');
ylabel('signals');
rx3 = x1 + (10^(0))*(randn(1,500));
figure();
plot(tnew,x1);
hold on;
plot(tnew, rx3);
legend('Original', 'Variance 10^0');
xlabel('t');
ylabel('signals');
snr1 = 10*log10(sum(x1.^2)/sum((rx1-x1).^2));
disp('With variance 10^(-4) snr is ');
disp(snr1);
snr2 = 10*log10(sum(x1.^2)/sum((rx2-x1).^2));
disp('With variance 10^(-2) snr is ');
disp(snr2);
snr3 = 10*log10(sum(x1.^2)/sum((rx3-x1).^2));
disp('With variance 10^(0) snr is ');
disp(snr3);

%Pe

m  = zeros(1,5);
for i = 1:5
    m(i) = sum(rx1(((i-1)*100)+1:(i*100)).*Y1)>0;
end
    

Eb = A11(1)^2;
varianc = [10^-4 10^-2 1];

for i = 1:3
    N0  = 2*(varianc(i));
    rb = Eb/N0;
    pe = qfunc(sqrt(2*rb));
    disp('with variance ');
    disp(varianc(i));
    disp(' Pe is ');
    disp(pe);
end


%% optimal reciever
randomn_bits = round(rand(1,10^5));
newrando_bits = ones(1,10^5);
newrando_bits(randomn_bits == 0) = -1;
new_x = zeros(1,10^7);

for i = 1:10^5
    new_x(((i-1)*100)+1:(i*100)) = newrando_bits(i)*A1(1:100);
end


variance = [10^-4 5^-4 10^-3 5^-3 10^-2 5^-2 10^-1 5^-1 10^0 5 10];
pe_new = zeros(1,size(variance,2));
pe_actual = zeros(1,size(variance,2));
m2  = ones(1,10^5);
snrs = zeros(1,size(variance,2));
for i = 1:size(variance,2)
    N0  = 2*(variance(i));
    rx = new_x + (variance(i))*(randn(1,10^7));
%     snrs = snr(new_x,rx);
    snrs(i) = 10*log10(sum(new_x.^2)/sum((rx-new_x).^2));
    %actual probability of error of signal
    
    for j = 1:10^5
        m2(j) = sum(rx(((j-1)*100)+1:(j*100)).*Y1)>0;
    end
    a = sum(m2==randomn_bits);
    pe_actual(i) = 1- sum(m2==randomn_bits)/(10^5);
%     
    %theoratical probability of error
    rb = Eb/N0;
    pe_new(i) = qfunc(sqrt(2*rb));
%     disp('with variance ');
%     disp(variance(i));
%     disp(' Pe is ');
%     disp(pe);
end

figure();
semilogy(snrs, 20*log10(pe_new), '-o');
xlabel('SNR');
ylabel('Probability of Error');
title('Pe vs SNR');
figure();
semilogy(snrs, 20*log10(pe_actual), '-o');
xlabel('SNR');
ylabel('Ratio of Error in Actual Predictions');
title('Ratio of Error vs SNR');
%% part g
randomg_bits = rand(1,10^5);
randg = zeros(1,10^5);
randg(randomg_bits<=0.1) = 1;
randg(randomg_bits>0.1) = -1;
new_xg = zeros(1,10^7);


for i = 1:10^4
    new_xg(((i-1)*100)+1:(i*100)) = randg(i)*A1(1:100);
end

peg_new = zeros(1,size(variance,2));
peg_actual = zeros(1,size(variance,2));
m2g  = ones(1,10^5);
snrsg = zeros(1,size(variance,2));
for i = 1:size(variance,2)
    N0  = 2*(variance(i));
    rxg = new_xg + (variance(i))*(randn(1,10^7));
%     snrs = snr(new_x,rx);
    snrsg(i) = 20*log(sum(new_xg.^2)/sum(rxg.^2));
    %actual probability of error of signal
    
    for j = 1:10^5
        m2(j) = sum(rxg(((j-1)*100)+1:(j*100)).*Y1)>0;
    end
    a = sum(m2==randomn_bits);
    peg_actual(i) = 1- sum(m2==randomn_bits)/(10^5);
%     
    %theoratical probability of error
    rb = Eb/N0;
    peg_new(i) = qfunc(sqrt(2*rb));
%     disp('with variance ');
%     disp(variance(i));
%     disp(' Pe is ');
%     disp(pe);
end

figure();
plot(snrsg, 20*log10(peg_actual));
xlabel('SNR');
ylabel('Ratio of Error');
title('Ratio vs SNR');


%% Part 3
random_bits = round(rand(1,10));

Fs = 5000;
T = 0.1;
t = 0:1/Fs:T-1/Fs;
tx = 0:1/Fs:5*T-1/Fs;
s = zeros(4,500);

s(1,:) = cos(2*pi*250*t);
s(2,:) = -1*s(1,:);
s(3,:) = cos(2*pi*2*250*t);
s(4,:) = -1*s(3,:);

xa = zeros(1,2500);
r = zeros(1,5);

r(1) = random_bits(2)*1 + random_bits(1)*2 + 1;
r(2) = random_bits(4)*1 + random_bits(3)*2 + 1;
r(3) = random_bits(6)* + random_bits(5)*2 + 1;
r(4) = random_bits(8)*1 + random_bits(7)*2 + 1;
r(5) = random_bits(10)*1 + random_bits(9)*2 + 1;

for i = 1:5
    xa(((i-1)*500)+1:(i*500)) = s(r(i),:);
end

figure();
plot(tx, xa);
xlabel('t');
ylabel('x(t)');


Y1 = sqrt(2/T)*cos(2*pi*250*t);
Y2 = sqrt(2/T)*cos(2*pi*2*250*t);
figure()
plot(t,Y1)
xlabel('t')
ylabel('Y1')

figure()
plot(t,Y2)
xlabel('t')
ylabel('Y2')
%mapping Y1 and Y2 plane
%%
sy1 = [sqrt(T/2)  0];
sy2 = [-sqrt(T/2) 0];
sy3 = [0 sqrt(T/2)];
sy4 = [0 -sqrt(T/2)];
figure();
plot(sy1(1), sy1(2), '*b');
hold on;
plot(sy2(1), sy2(2), '+r');
plot(sy3(1), sy3(2), 'og');
plot(sy4(1), sy4(2), 'pc');
legend('s1','s2','s3','s4');
xlabel('Y1');
ylabel('Y2');
title('Orthogonal Signal Space');

variance = [10^-2 10^0 10^2];
rx =zeros(3,2500);
for i = 1:3
    rx(i,:) = xa + variance(i)*(randn(1,2500));
    snr(i) = 10*log10(sum(xa.^2)/sum((rx(i,:)-xa).^2))
end
disp(snr)
figure();
plot(tx, rx(1,:));
hold on;
plot(tx,xa);
legend('Original', 'Variance 10e-2');
xlabel('t');
ylabel('signals');
figure();
plot(tx, rx(2,:));
xlabel('t');
ylabel('signals');
hold on;
plot(tx,xa);
legend('Original', 'Variance 10e0');
figure();
plot(tx, xa);
hold on;
plot(tx, rx(3,:));
legend('Original', 'Variance 10^2');
xlabel('t');
ylabel('signals');
%optimal reciever
randomr_bits = round(rand(1,2*10^5));
s_r = zeros(1,10^5);
j = 1;
xf =zeros(1,500*10^5);
for i = 1:1*10^5
    s_r(i) = randomr_bits(j+1)*1 +randomr_bits(j)*2 + 1;
    j = j + 2;
    xf(((i-1)*500)+1:(i*500)) = s(s_r(i),:);
end



%%
%noisy signal
Eb = 1;

% new_x = zeros(1,10^7);

variance = [10^-2 5^-2 10^-1 5^-1 10^0 5 10 30 50 100];
% variance = [10^-4];
pe_new = zeros(1,size(variance,2));
pe_actual = zeros(1,size(variance,2));
m1  = ones(1,10^5);
m2  = ones(1,10^5);
snrs = zeros(1,size(variance,2));
for i = 1:size(variance,2)
    N0  = 2*(variance(i));
    rx = xf + (variance(i))*(randn(1,500*10^5)); %generating
    snrs(i) = 10*log10(sum(xf.^2)/sum((rx-xf).^2));
    %Uncoment for ratio of error
    %actual probability of error of signal
    for j = 1:10^5
        rr1(j) = sum(rx(((j-1)*500)+1:(j*500)).*Y1);
        rr2(j) = sum(rx(((j-1)*500)+1:(j*500)).*Y2);
    end
    lbl = rr1>=rr2 & rr1>=-rr2;
%     lbl(lbl == 1) = 1;
    lbl2 = 3*(rr1<rr2 & rr1>= -rr2);
    lbl(lbl2==1) = 3;
    lbl3 = 2*(rr1<rr2 & rr1<-rr2);
%     lbl(lbl2 == 1) = 2;
    lbl(lbl == 0) = 4;
    lbl = lbl+lbl2+lbl3;
    a == sum(lbl==s_r);
%     a = sum(m2==randomn_bits);
    pe_actual(i) = 1 - sum(a)/(10^5);
     
    %theoratical probability of error
    rb = Eb/N0;
    pe_new(i) = qfunc(sqrt(rb));
end

figure();
semilogy(snrs, 20*log10(pe_new), '-o');
xlabel('SNR');
ylabel('Probability of Error');
title('Pe vs SNR');
figure();
semilogy(snrs, 20*log10(pe_actual), '-o');
xlabel('SNR');
ylabel('Ratio of Error in Actual Predictions');
title('Ratio of Error vs SNR');

% m  = zeros(1,2*10^5);
% for i = 1:5
%     ds = sum(rx1(((i-1)*500)+1:(i*500)).*Y1)>0;
% end
