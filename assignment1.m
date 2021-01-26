%% Part 1
close all; 
clear all;
rand('seed', sum(100*clock));
% area13 = (1000000- (1000000/6))/2; % dividing samples according to probabilities found in part 1
% % y = linspace(-3,6,1000000);
% area1 = linspace(-3,0, area13); % -3<y<0
% area3 = linspace(3,6, area13); % 3<y<6
% area123 = 1000000 - (area13*2);
% area2 = linspace(0,3, area123);

% Generating random samples and histogram and avg src power
y = linspace(-3,5,1000000); 
fy = zeros(1,1000000);
fy(y<0.0000000 &y>=-3.00000) = (y(y<0.000000& y>=-3.000000)/15) +1/5;
fy(y>=0.0000000 &y<2.00000) = 1/5; 
fy(y<=5.0000000 &y>=2.00000) = -(y(y<=5.000000& y>=2.000000)/15) +1/3;
sampled_y = (randsample(1000000,1000000,true, fy)./1000000)*8 -3 ;
figure();
histogram(sampled_y)
ylabel('Number of samples');
xlabel('Sample position');
title('Histogram of y samples');

fys = zeros(1,1000000);%prob according to random samples
fys(sampled_y<0.0000000 &sampled_y>=-3.00000) = (sampled_y(sampled_y<0.000000& sampled_y>=-3.000000)/15) +1/5;
fys(sampled_y>=0.0000000 &sampled_y<2.00000) = 1/5; 
fys(sampled_y<=5.0000000 &sampled_y>=2.00000) = -(sampled_y(sampled_y<=5.000000& sampled_y>=2.000000)/15) +1/3;


avg_srcpower = mean(sampled_y.^2);
disp('Avg source power is ');
disp(avg_srcpower);
% sampled_y = zeros(1, 1000000);
% sampled_y(1,1:416666) = sampled_y1;
codes = {'0','100','110','1010','1011','1110','11110','11111'};
% sampled_y(1, 416666+166666:416665+416666+166666) = sampled_y3;
% sampled_y(1,416667:416666+166666) = sampled_y2;
% y_ones = ones(1,999996);


%% Quantization
n_level = zeros(1,7);
quantized_sampled = sampled_y;
start = 0;
for quant = -3:4
    if quant == -3
        start = quant;
    end
    stop = start + 1;
    quantized_sampled(quantized_sampled <= stop & quantized_sampled>start) = (start+stop)/2;
    n_level(1,quant+4) = sum(quantized_sampled <= stop & quantized_sampled>start,'all');
    start = stop;
end

figure();
histogram(quantized_sampled);
ylabel('Number of samples');
xlabel('Sample position');
title('Histogram of quantized y samples');

% error power
sqr_er_power =mean((quantized_sampled - sampled_y).^2);
% sum(((quantized_sampled - sampled_y).^2).*transpose(fys),'all')/1000000;
disp('Error power is ')
disp(sqr_er_power);
disp('SQNR is ');
disp(10*log10(sqr_er_power/5^2) +6*8 + 4.8);

%% PMF & Huffman coding/decoding

pmf = n_level/1000000;
figure();
pmf_x = linspace(-2.5,4.5,8);
bar(pmf_x,pmf);
xlabel('y');
ylabel('PMF');
title('Y vs PMF');

%coding
sorted_pmf = sort(pmf);
new = sorted_pmf;
n1 = '';%a
n2 = '';%b
n3 = '';%c
n4 = '';%d
n5 = '';%e
n6 = '';%f
n7 = '';%g
n8 = '';%h
% interconnected = zeros(1,8);
% interconnected2 = zeros(1,8);
% interconnected3 = zeros(1,8);
% interconnected4 = zeros(1,8);
% interconnected5 = zeros(1,8);
% interconnected6 = zeros(1,8);
% interconnected7 = zeros(1,8);
% interconnected8 = zeros(1,8);
pmf1 = sorted_pmf(1);
pmf2 = sorted_pmf(2);
pmf3 = sorted_pmf(3);
pmf4 = sorted_pmf(4);
pmf5 = sorted_pmf(5);
pmf6 = sorted_pmf(6);
pmf7 = sorted_pmf(7);
pmf8 = sorted_pmf(8);


% frst_min = find(new == min(new));
% new(frst_min) = 2;

% scnd_min = find(new == min(new));
% 
% if scnd_min == 1
%     n1 = append(n1,'0');
%     interconnected(i)
% elseif scnd_min == 2
%     n2 = append(n2,'0');
% elseif scnd_min == 3
%     n3 = append(n3,'0');
% elseif scnd_min == 4
%     n4 = append(n4,'0');
% elseif scnd_min == 5
%     n5 = append(n5,'0');
% elseif scnd_min == 6
%     n6 = append(n6,'0');
% elseif scnd_min == 7
%     n7 = append(n7,'0');
% elseif scnd_min == 8
%     n8 = append(n8,'0');
% end

for i = 1:8
    if sorted_pmf(i) == pmf1
        n1 = append(n1,codes{1});
    elseif sorted_pmf(i) == pmf2
        n2 = append(n2,codes{2});
    elseif sorted_pmf(i) == pmf3
        n3 = append(n3,codes{3});
    elseif sorted_pmf(i) == pmf4
        n4 = append(n4,codes{4});
    elseif sorted_pmf(i) == pmf5
        n5 = append(n5,codes{5});
    elseif sorted_pmf(i) == pmf6
        n6 = append(n6,codes{6});
    elseif sorted_pmf(i) == pmf7
        n7 = append(n7,codes{7});
    elseif sorted_pmf(i) == pmf8
        n8 = append(n8,codes{8});
    end
end
% new_p = new(frst_min) + new(scnd_min);
% new(scnd_min) = new_p;

%decoder


keysets = {'a','b','c','d','e','f', 'g', 'h'};
vals = {n1,n2,n3,n4,n5,n6,n7,n8};

decodee = containers.Map(keysets,vals);

strng = {'b','a','d','g'};
for i = 1:length(strng)
    decoded = decodee(strng{i});
    disp(decoded);
end
disp('are codes for b a d g');
%average code word length
avg_code_length = 0; 
withoutcode = 0;
for i = 1:length(vals)
    avg_code_length = avg_code_length + strlength(vals{i})*sorted_pmf(i);
%     disp(strlength(vals{i}));
    withoutcode = withoutcode + 3*sorted_pmf(i);
end
disp('Avg codeword length with Huffman is ');
disp(avg_code_length);
disp('Avg codeword length without Huffman is ');
disp(withoutcode);

%% Applying PCM on non-random y first for testing
%cdf and inv_cdf
F_cdf = zeros(1,1000000);
inv_cdf = zeros(1,1000000);
% syms t;
F_cdf(y<=0.0000000 &y>=-3.00000) = (y(y<=0.0000000 &y>=-3.00000) + 3).^2/30;
inv_cdf(y<=0.0000000 &y>=-3.00000) = real(30^(1/2).*y(y<=0.0000000 &y>=-3.00000).^(1/2) - 3);
initial_cdf = F_cdf(375000);
initial_inv = inv_cdf(375000);
F_cdf(y>0.0000000 &y<=2.00000) = initial_cdf + y(y>0.0000000 &y<=2.00000)./5;
inv_cdf(y>0.0000000 &y<=2.00000) = initial_inv + real(5.*y(y>0.0000000 &y<=2.00000));
initial_cdf2 = F_cdf(624994);
initial_inv2 = inv_cdf(624994);
F_cdf(y<=5.0000000 &y>2.00000) = initial_cdf2 - ((y(y<=5.0000000 &y>2.00000) - 2).*(y(y<=5.0000000 &y>2.00000) - 8))./30;
inv_cdf(y<=5.0000000 &y>2.00000) =initial_inv2 + real(15*(1/25 - (2*y(y<=5.0000000 &y>2.00000))/15).^(1/2) + 5);
qy = y.*F_cdf.*inv_cdf; 
figure();
plot(y, qy);
xlabel('y');
ylabel('Q(y)');
title('Samples are not random in this one');

%% Applying PCM on random samples
%PCM

%cdf and inverse calculated using int() function and inverse in the command
%window 
F_cdf = zeros(1,1000000);
inv_cdf = zeros(1,1000000);
% syms t;
F_cdf(sampled_y<=0.0000000 &sampled_y>=-3.00000) = (sampled_y(sampled_y<=0.0000000 &sampled_y>=-3.00000) + 3).^2/30;
initial_cdf = 0.3;
initial_inv = -3;
F_cdf(sampled_y>0.0000000 &sampled_y<=2.00000) = initial_cdf + sampled_y(sampled_y>0.0000000 &sampled_y<=2.00000)./5;
initial_cdf2 = 0.7;
initial_inv2 = 6.9997;
F_cdf(sampled_y<=5.0000000 &sampled_y>2.00000) = initial_cdf2 - ((sampled_y(sampled_y<=5.0000000 &sampled_y>2.00000) - 2).*(sampled_y(sampled_y<=5.0000000 &sampled_y>2.00000) - 8))./30;



compressor = F_cdf;
comprsd_sample = F_cdf .* transpose(sampled_y);
quantized_comsampled = comprsd_sample;
figure();
histogram(comprsd_sample);
xlabel('y');
ylabel('Number of samples');
title('Compressed Samples')

for quant = 0:8
    if quant == 0
        start = quant;
    end
    stop = start + 0.625;
    quantized_comsampled(quantized_comsampled <= stop & quantized_comsampled>start) = (start+stop)/2;
%     n_level(1,quant+4) = sum(quantized_comsampled <= stop & quantized_comsampled>start,'all');
    if quant == 0
        quantized_comsampled(quantized_comsampled <= 0) = (start+stop)/2;
    end
    start = stop;
    
end

figure();
histogram(quantized_comsampled);
xlabel('y');
ylabel('Number of samples');
title('Compressed Quantized Samples')

inv_cdf(quantized_comsampled<=0.0000000 &quantized_comsampled>=-3.00000) = real(30^(1/2).*quantized_comsampled(quantized_comsampled<=0.0000000 &quantized_comsampled>=-3.00000).^(1/2) - 3);
inv_cdf(quantized_comsampled>0.0000000 &quantized_comsampled<=2.00000) = initial_inv + real(5.*quantized_comsampled(quantized_comsampled>0.0000000 &quantized_comsampled<=2.00000));
inv_cdf(quantized_comsampled<=5.0000000 &quantized_comsampled>2.00000) =initial_inv2 + real(15*(1/25 - (2*quantized_comsampled(quantized_comsampled<=5.0000000 &quantized_comsampled>2.00000))/15).^(1/2) + 5);
expanded = quantized_comsampled.*inv_cdf;
error_qy = mean((transpose(expanded)-sampled_y).^2);
disp('Error from Q(y) is ');
disp(error_qy);
sqnr_qy = 10*log10(error_qy/56^0)+ 6*8 +4.8;
disp('SQNR from Q(y) is ');
disp(sqnr_qy);
figure();
expands = histogram(expanded);
xlabel('y');
ylabel('Number of samples');
title('Expanded Quantized Samples');
valuesss = expands.Values;
oncemore = zeros(1,8);
oncemore = valuesss(valuesss~=0);
total = sum(oncemore,'all');
pmfs_qy = oncemore/total;
q_x = [-0.5 1.5 7.5 26 33.5 41 48.5 56];

figure();
plot(sampled_y, expanded);
xlabel('y');
ylabel('Q(y)');


figure();
bar(q_x, pmfs_qy);
xlabel('output of expander');
ylabel('PMF');



%% PMF & Huffman coding/decoding for PCM


%coding
sorted_pmf = sort(pmfs_qy);
new = sorted_pmf;
n1 = '';%a
n2 = '';%b
n3 = '';%c
n4 = '';%d
n5 = '';%e
n6 = '';%f
n7 = '';%g
n8 = '';%h
% interconnected = zeros(1,8);
% interconnected2 = zeros(1,8);
% interconnected3 = zeros(1,8);
% interconnected4 = zeros(1,8);
% interconnected5 = zeros(1,8);
% interconnected6 = zeros(1,8);
% interconnected7 = zeros(1,8);
% interconnected8 = zeros(1,8);
pmf1 = sorted_pmf(1);
pmf2 = sorted_pmf(2);
pmf3 = sorted_pmf(3);
pmf4 = sorted_pmf(4);
pmf5 = sorted_pmf(5);
pmf6 = sorted_pmf(6);
pmf7 = sorted_pmf(7);
pmf8 = sorted_pmf(8);
% codes = {'0','100','110','1010','1011','1110','11110','11111'};

% frst_min = find(new == min(new));
% new(frst_min) = 2;

% scnd_min = find(new == min(new));
% 
% if scnd_min == 1
%     n1 = append(n1,'0');
codes = {'0','100','110','110','100','101','1110','11110','11111'};
%     interconnected(i)
% elseif scnd_min == 2
%     n2 = append(n2,'0');
% elseif scnd_min == 3
%     n3 = append(n3,'0');
% elseif scnd_min == 4
%     n4 = append(n4,'0');
% elseif scnd_min == 5
%     n5 = append(n5,'0');
% elseif scnd_min == 6
%     n6 = append(n6,'0');
% elseif scnd_min == 7
%     n7 = append(n7,'0');
% elseif scnd_min == 8
%     n8 = append(n8,'0');
% end

for i = 1:8
    if sorted_pmf(i) == pmf1
        n1 = append(n1,codes{1});
    elseif sorted_pmf(i) == pmf2
        n2 = append(n2,codes{2});
    elseif sorted_pmf(i) == pmf3
        n3 = append(n3,codes{3});
    elseif sorted_pmf(i) == pmf4
        n4 = append(n4,codes{4});
    elseif sorted_pmf(i) == pmf5
        n5 = append(n5,codes{5});
    elseif sorted_pmf(i) == pmf6
        n6 = append(n6,codes{6});
    elseif sorted_pmf(i) == pmf7
        n7 = append(n7,codes{7});
    elseif sorted_pmf(i) == pmf8
        n8 = append(n8,codes{8});
    end
end
% new_p = new(frst_min) + new(scnd_min);
% new(scnd_min) = new_p;

%decoder


keysets = {'a','b','c','d','e','f', 'g', 'h'};
vals = {n1,n2,n3,n4,n5,n6,n7,n8};

decodee = containers.Map(keysets,vals);

strng = {'b','a','d','g','f','c'};
for i = 1:length(strng)
    decoded = decodee(strng{i});
    disp(decoded);
end
disp('are codes for b a d g f c');
%average code word length
avg_code_length = 0; 
withoutcode = 0;
for i = 1:length(vals)
    avg_code_length = avg_code_length + strlength(vals{i})*sorted_pmf(i);
%     disp(strlength(vals{i}));
    withoutcode = withoutcode + 3*sorted_pmf(i);
end
disp('Avg codeword length with Huffman is ');
disp(avg_code_length);
disp('Avg codeword length without Huffman is ');
disp(withoutcode);


%% lloyd Max

% n_level = zeros(1,7);
lloyd_sampled = sampled_y;
start = -3.5;
% mid = zeros(1,10);
mse = ones(1,10);
i = 1;
% mid(1) =-3.5;
mid2 = -3.5;
% start = -3.5
j =1;
a= -3.5;
while j<10
for quant = a:a+8
%     if quant == -3.5 
%         start = quant;
%     end
    start = mid2;
    stop = start + 1;
%     mid(i) = (start+stop)/2
    if quant == a
        mid2 = (start+stop)/2;
%         a = mid2;
    end
    lloyd_sampled(lloyd_sampled <= start & lloyd_sampled>stop) = (start+stop)/2;
%     n_level(1,quant+4) = sum(lloyd_sampled <= stop & lloyd_sampled>start,'all');
    start = stop;
%     i = i + 1;
end
mse(j) = mean((lloyd_sampled - sampled_y).^2);
% figure();
j = j + 1;
 a = mid2;
% i = 1;
% histogram(lloyd_sampled);
end

%% Part2

%this str was generated using upper() and strrep()
% str = "THE MANIFESTO IS AIMED FOR THE PROLETARIANS IT AIMS TO CONVINCE THEM HOW THE COMMUNIST ARE THEIR TRUE REPRESENTATION AND HOW THEY CAN HELP THEMSELVES BY OVERTHROWING THE BOURGEOIS CLASS IT ALSO PRESENTS THE BOURGEOIS AS SOME EVIL CLASS WHOSE AIM IS TO PILE UP ITS OWN CAPITAL AND CONTINUE UNDERMINING THE PROLETARIANS  THE PROLETARIANS STRUGGLE CAN BE COMPARED WITH THE STRUGGLES OF THE CREATURE IN MARY SHELLEYS FRANKENSTEIN THOUGH THERE WAS NO POLITICAL ELEMENT IN THE NOVEL BUT THE SOCIAL STRUGGLES OF THE MONSTER ITS INTERACTION WITH THE OTHER HUMANS CAN BE COMPARED TO THE STRUGGLES FACED BY THE PROLETARIANS IF WE ARE REPRESENTING THE CREATURE AS A PROLETARIAT THEN WHO IS THE BOURGEOIS OF THE STORY WELL IT CAN BE SAID VICTOR IS A PART OF THE BOURGEOIS CLASS VICTORS COMPARISON TO THE BOURGEOIS CAN BE BOTH SOCIAL AND POLITICAL FIRSTLY HE IS FROM A WEALTHY FAMILY AND HIS FAMILY OWNS A HOUSE WHICH IS PRIVATE PROPERTY SECONDLY HE ENJOYS THE SERVICE OF SERVANTS WHICH MEANS HE HAS POWER OVER SOME WORKING CLASS PEOPLE OTHER THAN THESE FINANCIAL AND POLITICAL REASONS HIS TREATMENT OF THE CREATURE ALSO SHOWS HIM AS SOME EVIL PERSON WHO ISNT WILLING TO GIVE THE CREATURE A PROLETARIAT A CHANCE IN RELATION TO THE QUOTE THE QUOTE MENTIONS BOURGEOIS SOCIETY AS SOMEONE WHO IS NO LONGER ABLE TO CONTROL THE POWERS OF THE NETHER WORLD WHOM HE HAS CALLED UP BY HIS SPELLS SO THIS QUOTE CAN ALSO REPRESENT AS FRANKENSTEIN WHO IS A BOURGEOIS BEING UNABLE TO CONTROL SOMETHING WHICH HE HAS CALLED UP BY HIS OWN SPELLS HIS CREATION THE CREATURE AS THE CREATURE WAS CREATED HE WAS IMMEDIATELY ABANDONED BY HIS CREATOR FRANKENSTIEN HE LATER HAS EXPERIENCES WITH OTHER HUMANS WHICH SHOWS HIS STRUGGLES FOR EXAMPLE THE MONSTER SAYS THAT HE ENTERS A HUT FOR FOOD BUT THE RESIDENT AND OLD MAN RUNS AWAY WHEN HE SEES THE CREATURE AND ALSO PEOPLE FLEE AT THE SIGHT OF HIM WHEN HE ENTERS A VILLAGE THIS IS ONE OF HIS FIRST EXPERIENCES WHICH HE NARRATES TO VICTOR WHICH SHOWS HIS EARLY STRUGGLES SIMILARLY IT IS MENTIONED THE PROLETARIAT GOES THROUGH VARIOUS STAGES OF DEVELOPMENT WITH ITS BIRTH BEGINS ITS STRUGGLE WITH BOURGEOIS  THIS EXPERIENCE OF HIS CAN ALSO BE COMPARED TO HOW THE WORKERS THE PROLETARIATS ARE TREATED AS SLAVES AND THEY ARE GIVEN NO RESPECT THE CREATURE AT THIS MOMENT CAN BE IDENTIFIED AS PROLETARIAT BUT STILL NOT A COMMUNIST BECAUSE COMMUNISM PROMOTES THE USE OF ANTAGONISTIC BEHAVIOR FOR THE RIGHTS OF THE PROLETARIAT BUT THE CREATURE IS STILL GENTLE TO THE HUMANS WHO TREAT HIM BADLY I THINK HE COULD BE SEEN AS PERSON IDENTIFYING WITH THE SOCIALIST PRINCIPLES AS HE STILL WANTS TO BEFRIEND HUMANS WHICH HE TRIES TO DO HIS ATTEMPT GOES WRONG AND THIS PUSHES HIM OFF THE EDGE NOW HE SWEARS TO REVENGE HIMSELF AGAINST ALL HUMANS AND HE KILLS VICTORS BROTHER WILLIAM#";
str = "MY NAME IS ABDUL MOEED AND THIS IS THE SENTENCE I AM TESTING FOR THE REPORT. MY CODE IS WORKING #";
str2Test = "TEST TE# ";

valueSet = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
names = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ' '};
M = containers.Map(names,valueSet);
b = 0;
a = 1;
i = 1;
codecell = {};
 while 1
        count = 0;
       value = int64(length(valueSet));
        while value > 0
          count = count + 1;
           value = bitsra(value,1);
        end
%     
    while 1
        newStr = extractBetween(str,a,a+b);
        if  M.isKey(newStr) == 1
            b = b+1;
            newStr2 = extractBetween(str,a,a+b); 
            if M.isKey(newStr2) == 1
                b = b+1;
                continue               
            else
              valueSet(length(valueSet)+1) = length(valueSet)+1;
              names(length(names) + 1) = cellstr(newStr2(:));
              x = M.values(cellstr(newStr(:)));
              binVal = decimalToBinaryVector(x{1},count);
              
              binstr = num2str(binVal);
              binstr(isspace(num2str(binstr))) = '';
              codecell{i} = binstr;
              i =i + 1;
              break;
                
            end               
        else
        valueSet(length(valueSet)+1) = length(valueSet)+1;
        names(length(names) + 1) = cellstr(newStr(:));
        x = M.values(cellstr(newStr2(:)));
        binVal = decimalToBinaryVector(x{1},count);
        binstr = num2str(binVal);
        binstr(isspace(num2str(binstr))) = '';
        codecell{i} = binstr;
        i =i + 1;
        break;
        end
    end
      a = a+b;
      b = 0;
      M = containers.Map(names,valueSet);
     if  newStr == "#" | newStr2 == "#"
         break;
     end
 end

%% avg code word length
loopi = length(codecell);
codetot = 0;
for j = 1:loopi
    
 codetot = codetot + strlength(codecell{j});
end
avg_code = codetot/loopi;
disp('Avg length is ');
disp(avg_code);

%% decoding
a = {};
for j = 1:loopi
 posi = bin2dec(codecell{j});
%  a{1} = a{1}+ names(posi);
 disp(names(posi));
end