M = csvread('logicanalyzer_noheader.csv');
N = length(M);
dts = [];
reptime = 10^-8;
lastpulse = 0;
for i = 2:N
   if xor(M(i,2) > M(i-1, 2), M(i,3) > M(i-1, 3))
       if M(i) - lastpulse < reptime
           dts = [dts, M(i) - lastpulse];
       end
       lastpulse = M(i);
   else if M(i,2) > M(i-1, 2) && M(i,3) > M(i-1, 3)
           dts = [dts, 0];
           lastpulse = M(i);
       end
   end 
end
N
dts = dts * 10^9;
length(dts)
hist(dts)

