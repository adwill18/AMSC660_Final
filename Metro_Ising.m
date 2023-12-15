function [m,m_bar,mu,var_m]=Metro_Ising(beta)
N=30;
itermax=10^8;
% Generate inital magnetic configurations
% spin=rand(N);
% upi=(spin >= 0.5);
% downi=(spin<0.5);
% spin(upi)=1;
% spin(downi)=-1;
spin=ones(N);

msum=0;
vsum=0;
m = sum(sum(spin))/(N^2); % Find the magnetization
mu=m;
for iter=1:itermax
    %Randomly pick a site
   index=randi(numel(spin));
   [L ,k]=ind2sub(size(spin),index);

   %Find the four nearest neighbors
   above = mod(L -1 -1, size(spin,1)) +1;
   below = mod(L +1 -1, size(spin,1)) +1;
   left     = mod(k -1 -1, size(spin,2)) +1;
   right   = mod(k +1 -1, size(spin,2)) +1;
   neighbors = [spin(above,k); spin(L, left); spin(L,right); spin(below,k)];

   %Calculate the energy diffrence deltaH
   dH= 2* spin(L,k) * sum(neighbors);

   if dH <= 0
       accept = 1;
   else
       u = rand;
       prob = exp(-beta*dH);
       if  u <= prob
           accept = 1;
       else
           accept = 0;
       end
   end
   if accept ==1
       spin(L,k) = - spin(L,k); % Flip the spin that was proposed to flip
       m = sum(sum(spin))/(N^2); % Calculate the magnetization of the new state
   end
   mu = (iter * mu+m)/(iter+1); % Update the mean magnetization 

   % Calculate the the running mean
   msum = msum+m;
   m_bar = (1/iter )*msum;

  % Calculate the running variance 
  vsum=vsum+(m-m_bar)^2;
  var_m= (1/(iter-1))*vsum;
  
end
end
