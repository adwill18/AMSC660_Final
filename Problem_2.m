%% Set up the initial configuration
N = 91; % the number of vertices
w=randn(2*N,1) *sqrt(N);
kmax = 1e6;
tol = 5*10^(-4);
%% Run Solver
[w, f, normgrad,k] = stochasticAdam(w,kmax,tol);

x=w(1:N,1);
y=w(N+1:2*N,1);
A=load("Adjacency_matrix.csv");

%figure(1)
plot_graph(x,y,A)

kvals = 1:k;

figure(2)
plot(kvals,normgrad(1:k),'LineWidth',2)
set(gca, 'YScale','log')
xlabel('Iteration #','FontSize',20)
ylabel('Norm of Force','FontSize',20)
title('Norm of Force')

%%
function plot_graph(x,y,A)
figure;
hold on
plot(x,y,'o','Markersize',15,'MarkerEdgeColor',[0.5,0,0],'MarkerFaceColor',[1,0,0]);
ind = find(A == 1);
[I,J] = ind2sub(size(A),ind);
for k = 1 : length(ind)
    plot([x(I(k)),x(J(k))],[y(I(k)),y(J(k))],'linewidth',4,'Color',[0,0,0.5]);
end
daspect([1,1,1])
axis off
hold off
end
