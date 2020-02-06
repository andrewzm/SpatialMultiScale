%%%%Negative loglikelihood function for the FSA-Block approach
%%%%Outputs:
%%%%nll: negative loglikelihood value
%%%%zsz: quadratic term: t(z)*inv(C)*z
%%%%logdet: logrithm of determinant
%%%%U: cross-covariance matrix of the training set and the knot set
%%%%B: covariance of the knot set
%%%%Cs: Cs{i} is the residual covariance matrix of block i
%%%%Q:  an m by m matrix in Woodbury inversion formula: B + t(U)*inv(Cs)*U
%%%%prodVP1: U*inv(Cs)*z

function [nll, zsz, logdet, U, B, Cs, Q, prodVP1] = nllik_Block_space(blockid, para, z, alld_sk, alld_k, alld_list, funname)
 
paras=exp(para);

switch funname
    
    
    case 'matern'
        
         theta.phi=paras(1); theta.nu=3*paras(2)/(1+paras(2)); 
         theta.epsilon=paras(3); theta.sigma2=paras(4);

    case 'matern1'
        
         theta.phi=paras(1); theta.nu=1; 
         theta.epsilon=paras(2); theta.sigma2=paras(3);
        
    case 'exp'
        
        theta.beta=paras(1); theta.epsilon=paras(2); theta.alpha=paras(3);
          
        
end

 blockidd=blockid ;  

 
 Bsize=max(blockidd);
 ctB=hist(blockidd,1:Bsize);
 temp=[0 cumsum(ctB)]; 
 
 Q_term2=0;
 

                 
    U=CovFun(theta, alld_sk, funname);
          
                  
    B=CovFun(theta, alld_k, funname);
  

 
    B=B+exp(-10).*eye(size(B,1));

  prodVP1=0;Cs={};
  hAy=[];logdetA=0;
hB = chol(B)';


 for i=1:Bsize 
    
 yB=z(temp(i)+1:temp(i+1),:);
 UB=U(temp(i)+1:temp(i+1),:); 
        


 d_XB=alld_list{i}; 
 
 cXB=CovFun(theta, d_XB, funname);

 %%%%Covariance matrix of ith data block
C_B=cXB+theta.epsilon*eye(size(cXB,1));  %Add nugget

hBU = hB\UB'; 
Appr_str=hBU'*hBU; 

%%%%AB: residual covariance matrix of ith data block
AB=C_B-Appr_str; 

 hAB = chol(AB)'; 
 hAUB=hAB\UB;
 hAy(temp(i)+1:temp(i+1),:)=hAB\yB;
 prodVP1=prodVP1+hAUB'*hAy(temp(i)+1:temp(i+1),:);
 logdetA=logdetA+(2 * sum(log(diag(hAB))));
 
 Q_term2=Q_term2+hAUB'*hAUB; 
 Cs{i}=AB;
 
 end
  
 Q=B+Q_term2; hQ=chol(Q)';hQWy=hQ\prodVP1;
  
 zsz=sum(-sum(hQWy.^2)+sum(hAy.^2));
 
 logdetQ = 2 * sum(log(diag(hQ)));
 logdetB= -2 * sum(log(diag(hB)));
 logdet=logdetA+logdetQ+logdetB;
 
 nll =(zsz  + size(z,2).*logdet)/2;
 

 
 
 
