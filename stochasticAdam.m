function [w, f, normgrad,k] = stochasticAdam(w,kmax,tol)
    N = 91; % the number of vertices
    x=w(1:N,1);
    y=w(N+1:2*N,1);
    A=load("Adjacency_matrix.csv");

    beta1=0.9;
    beta2=0.999;
    eps=10^-8;
    eta=0.001;
    normgrad = zeros(kmax,1);
    f = zeros(kmax + 1,1);
    k = 1;
    lowestNG=100;
    normgrad(k) = 1;
    w0=w;
    w=w0;
    v=0;
    m=0;
    for k=1:kmax
        g=-forces(x,y,A); % estimate the gradient based on the vector
        f(k)=energy(x,y,A); % estimate the value of the function
        normgrad(k)=norm(g); % norm of the estimated gradient
        m=beta1*m+(1-beta1)*g;
        v=beta2*v+(1-beta2)*(g.*g);
        m_hat=m/(1-beta1);
        v_hat=v/(1-beta2);
        w=w-eta*m_hat./(sqrt(v_hat)+eps);
        x=w(1:N,1);
        y=w(N+1:2*N,1);

        if mod(k,6000)==0
             fprintf(' k = %d,  f = %d, ||g|| = %d \n',k,f(k),normgrad(k));
        end
        
        if normgrad(k) < tol
            fprintf(' k = %d,  f = %d, ||g|| = %d\n',k,f(k),normgrad(k));
            break;
        end   
    end
   
    if normgrad(k)< lowestNG
        bestW=w;
        bestF=f;
    end
    
    w= bestW;
    f=bestF;
end

function f = forces(x,y,A)
% f = force = - grad U = column vector with 2*N components
% x, y are column vectors with N components
% A is an N-by-N adjacency matrix
N = length(x);
%% find pairwise distances between linked vertices
xaux = x*ones(size(x))';
yaux = y*ones(size(y))';
dx = A.*xaux - A.*(xaux'); 
dy = A.*yaux - A.*(yaux');
dxy = sqrt(dx.^2+dy.^2);

%% spring forces due to linked vertices
Aind = find(A == 1);
idiff = zeros(N);
idiff(Aind) = 1 - 1./dxy(Aind);
fx = -sum(idiff.*dx,2);
afx = min(abs(fx),1);
sfx = sign(fx);
fx = afx.*sfx;

fy = -sum(idiff.*dy,2);
afy = min(abs(fy),1);
sfy = sign(fy);
fy = afy.*sfy;

f_linked = [fx; fy];

%% repelling spring forces due to unlinked vertices
h = sqrt(3);
Aind = find(A==0);
A = ones(size(A))-A;
dx = A.*xaux - A.*(xaux'); 
dy = A.*yaux - A.*(yaux');
dxy = sqrt(dx.^2+dy.^2);
fac = zeros(N);
diff = dxy - h;
fac(Aind) = min(diff(Aind),0); 
fx = sum(fac.*dx,2);
fy = sum(fac.*dy,2);
f_unlinked = -[fx; fy];

f = f_linked + f_unlinked;
end
%%
function U=energy(x,y,A)
[iA, jA] = find(A == 1);
U_linked=sum((sqrt((x(iA)-x(jA)).^2+(y(iA)-y(jA)).^2) -1).^2);
[inE, jnE] = find(A==0);
hh = sqrt(3);
value=sqrt((x(inE)-x(jnE)).^2+(y(inE)-y(jnE)).^2) -hh;
Nind=length(iA);
U_unlinked=zeros(Nind,1);
for ii=1:Nind
    U_unlinked(ii)=min(value(ii),0);
end
U_unlinked=sum(U_unlinked.^2);
U= U_linked+U_unlinked;
end
