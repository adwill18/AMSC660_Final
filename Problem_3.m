beta = [0.2:0.01:1];
mu_calc=zeros(length(beta),1);
mu=zeros(length(beta),1);
m=zeros(length(beta),1);
m_bar=zeros(length(beta),1);
var_m=zeros(length(beta),1);

for k= 1:length(beta)
    [m(k),m_bar(k),mu(k),var_m(k)]=Metro_Ising(beta(k));
    if beta(k)>0.4408
        mu_calc(k)=(1-sinh(2*beta(k))^(-4))^(1/8);
    else
        mu_calc(k)=0;
    end
end

figure('Name','Problem 3')
plot(beta,mu_calc,'b','LineWidth',2)
hold on
plot(beta,mu,'r*')
plot(beta,m_bar+sqrt(var_m),'m--')
plot(beta,m_bar-sqrt(var_m),'m--')
title("Problem 3","FontSize",20)
xlabel("Beta","FontSize",20)
ylabel("Mean Magnetization","FontSize",20)
legend('predicted mu','computed mu','magnetization + ST deviation','magnetization - ST deviation')
hold off
