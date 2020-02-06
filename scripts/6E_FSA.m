%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Authors: Bohai Zhang and Andrew Zammit-Mangion
%%% Date: 24 January 2020
%%% Details: Generates predictions using the FSA
%%%          This code comes with no warranty or guarantee of any kind.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Add paths
addpath('./FSA/nllik')
addpath('./FSA/pred')
addpath('./FSA/Cov')
addpath('./FSA')

%%% Field names are lon-lat-bias-sst-error-lat2-z
fileID = fopen('cache_folder.txt','r');
cache_folder = fscanf(fileID,'%s')
y_tot = csvread(strcat(cache_folder, "/y_tot.csv"),1);
y_pred = csvread(strcat(cache_folder, "/y_pred.csv"),1);

disp(strcat('We have to do',{' '}, num2str(length(y_pred)),' predictions'));

n = length(y_tot) + length(y_pred);
ns = length(y_tot); 
np= length(y_pred);

Xs = y_tot(:,1:2);
zs = y_tot(:,7);

Xp = y_pred(:,1:2);
zp = y_pred(:,7);

%%% Use an Matern1 model.
funname='matern1';

%%% Use K-means algorithm to find K clusters for training data
%%% K: number of blocks
K=round(length(zs)/100);

%%% m: number of knots
m=100;

%%% Make knots on sphere
phi = rand(m,1)*2*pi;
costheta = rand(m,1)*2-1;
theta = acos(costheta);

x = sin(theta) .* cos(phi);
y = sin(theta) .* sin(phi);
z = cos(theta);

lat = 90 - acos(z)*360/(2*pi);
lon = atan2(y,x)*360/(2*pi)+180;

knots(:,1) = lon;
knots(:,2) = lat;

%%% Initial values
logtheta0 = log([1 0.1 2]);

%%% Specify tolerance, maximum number of iterations for matlab optimization function "fminunc"
options = optimset('TolX', 1e-6, 'MaxIter', 2000, 'MaxFunEvals', 2000);
 
%%% Run FSA
[theta, MSPE_FSA, pred_zp, zp_sort, loc_p, pred_sd]=Fun_FSAB(Xs, zs, Xp, zp, knots, K, funname, logtheta0);

%%% Save results
save(strcat(cache_folder, '/FSA_results.mat'),'theta','MSPE_FSA','pred_zp','loc_p','pred_sd');
dlmwrite(strcat(cache_folder, '/FSA_results.csv'),[loc_p,pred_zp,pred_sd],'precision',7);   
