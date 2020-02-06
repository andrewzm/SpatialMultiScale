function [theta, MSPE, pred_zp, zp, loc_p, Pred_SD]=Fun_FSAB(Xs, zs, Xp, zp, knots, K, funname, logtheta0)


%%%Use K-means algorithm to find K clusters for training data
%%%K: number of blocks
[IDx, center]=kmeans(Xs, K,'MaxIter',800);
%Sort the block indices and obtain new indices of data
%data are grouped according to block indices
[Blockid, IX]=sort(IDx);
%Sort the training data from the smallest block id to the largest block id
loc_s=Xs(IX,:); zs=zs(IX,:);

%%%Find Block id of prediction locations
%%%Block ID based on distances to block centers of the training sites
Dis_p=distance(Xp, center, 'euclidean');
Blockid_p=zeros(size(Xp, 1), 1);
np=size(zp, 1);
 
 for i=1:np
     
    [~, temp_ind_p]=min(Dis_p(i,:));
     
     Blockid_p(i,:)=temp_ind_p;
     
 end
 
%Sort predictive locations in an increasing order of Block id. 
[Blockid_p, IX_p] = sort(Blockid_p);
loc_p=Xp(IX_p,:);
zp=zp(IX_p,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bsize=max(Blockid);
cB=hist(Blockid,1:Bsize);
Bindex=[0 cumsum(cB)]; 

%%%Compute Euclidean distances for each location pair in a block
for l=1:(size(Bindex,2)-1)

loc_s_list{l}.loc=loc_s((Bindex(l)+1):Bindex(l+1),:);
alld_list{l}.alld=distance(loc_s_list{l}.loc, loc_s_list{l}.loc, 'euclidean');
    

end


%%%%Knots
%%%%Distance matrix between training sites and knot set
alld_sk.alld=distance(loc_s, knots, 'euclidean');
%%%%Distance matrix of the knot set
alld_k.alld=distance(knots, knots, 'euclidean');

%%%specify tolerance, maximum number of iterations for matlab optimization function "fminunc"
options=optimset('TolX',1e-6, 'MaxIter', 2000, 'MaxFunEvals',2000);
   
%%%FSA-Block
    tic;
    
    [xval,fval,exitflag,output,grad,hessian] = fminunc(@(paras) nllik_Block_space(Blockid, paras, zs, alld_sk, alld_k,...
                                                              alld_list, funname), logtheta0, options);
    toc;
    
    paras=exp(xval);
    
    switch funname
        
        case 'exp'
    
    
    theta.beta=paras(1); theta.epsilon=paras(2); theta.alpha=paras(3);
    
    
        case 'matern'

%   For Matern model    
    theta.phi=paras(1); theta.nu=3*paras(2)/(1+paras(2)); theta.epsilon=paras(3); theta.sigma2=paras(4);
    
    end

    %%%%Prediction
    [~, ~, ~, U, ~, Cs, Q, prodVP1] = nllik_Block_space(Blockid, xval, zs, alld_sk, alld_k, alld_list, funname);
    
    %%%Cp is the cross covariance matrix between predictive locations and the
    %%%training locations, using the FSA-Block covariance function. 
    Cp=Cov_FSABlock(theta, loc_s, loc_p, knots, Blockid, Blockid_p, funname);

    [pred_zp, MSPE] = Fun_pred_FSABlock(Blockid, zs, zp, alld_sk, prodVP1, Cs, Q, U, Cp) ;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Predictive standard errors
    Cpp=Cov_FSABlock(theta, loc_p, loc_p, knots, Blockid_p, Blockid_p, funname)+theta.epsilon*eye(np);

    %%%%obtain the predictive variance
    Temp1_Cp=0; Temp2_Cp=0;

    for i=1:max(Blockid)

       C_resid=Cs{i}; 
       QC_re=chol(Cs{i});
       invQC_reCp=QC_re'\Cp((Bindex(i)+1):Bindex(i+1),:);
       Temp1_Cp=Temp1_Cp+invQC_reCp'*invQC_reCp;

       invQC_reU=QC_re'\U((Bindex(i)+1):Bindex(i+1),:);
       Temp2_Cp=Temp2_Cp+invQC_reCp'*invQC_reU;

    end

    Pred_SD=sqrt(diag(Cpp-(Temp1_Cp-Temp2_Cp*(Q\(Temp2_Cp')))));



    switch funname
        
        case 'exp'
    
    %%%Save the results
      sprintf('sigma2:%3.4f,phi:%3.4f,nugget:%2.4f,MSPE:%5.5f', theta.alpha, theta.beta, theta.epsilon, MSPE)
         
        case 'matern'
      
      sprintf('sigma2:%3.4f,phi:%3.4f,nu: %2.4f, nugget:%2.4f,MSPE:%5.5f', theta.sigma2, theta.phi, theta.nu, theta.epsilon, MSPE)

      
    end
    
    