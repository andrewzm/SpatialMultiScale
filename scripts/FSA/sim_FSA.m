clc;
clear;

addpath('./nllik')
addpath('./pred')
addpath('./Cov/')

%%%Use an exponential model.
%%%2-dim example, [0 10]*[0  10]
n=1000;
ns=round(n*0.9); np=n-ns;

%%%Generate locations
X=10.*rand(n, 2);

d.alld=distance(X, X, 'euclidean');

%%%Generate responses
funname='exp';
theta.alpha=1; theta.beta=4; theta.epsilon=0.01;

MatC=CovFun(theta, d, funname)+theta.epsilon*eye(n);

z=mvnrnd(zeros(n,1),MatC)';


%%%Randomly partition the data into training and predition sets
[dum,I]    = sort(rand(n)); clear dum;
idx      = ones(n,1); idx(I(1:np)) = 0; 
idx=logical(idx);

%%%Xs: training locations, Xp: testing locations
Xs=X(idx,:);  Xp=X(~idx,:); 
zs=z(idx,:);  zp=z(~idx,:);

%%%Use K-means algorithm to find K clusters for training data
%%%K: number of blocks
K=10;

%%%%Knots
%%%%m: number of knots
m=100;
knots=10.*rand(m,2);

%%%Initial values
logtheta0=log([1 0.1 2]);
%%%specify tolerance, maximum number of iterations for matlab optimization function "fminunc"
options=optimset('TolX',1e-6, 'MaxIter', 2000, 'MaxFunEvals',2000);
 
[theta, MSPE_FSA, pred_zp, zp_sort]=Fun_FSAB(Xs, zs, Xp, zp, knots, K, funname, logtheta0);


       
   