%%%The cross-covariance matrix by the FSA-Block approach

function Cp=Cov_FSABlock(theta, loc_1, loc_2, knots, Blockid, Blockid_p, funname)

loc_s=loc_1; loc_p=loc_2;

np=size(loc_p,1);

ID_B=zeros(np, size(loc_s,1));
for i=1:np

    ID_B(i,Blockid==Blockid_p(i))=1;

end


   alld_pk=distance(loc_p, knots,'euclidean');
  
   alld_sp=distance(loc_s,loc_p,'euclidean');
  
   alld_k=distance(knots,knots,'euclidean');
  
   alld_sk=distance(loc_s,knots,'euclidean');
  

   d_pm.alld=alld_pk; 
   d_mm.alld=alld_k; 
   d_nm.alld=alld_sk; 
   d_np.alld=alld_sp; 
     

    
%%For other covariance functions
%%Cpm: cross-covariance between predictive sites and the knot sites
%%Cmm: covariance of the knot set
%%Cnm: cross-covariance beween training and knot sites
%%Cnp: cross-covariance between training and predictive sites

Cpm=CovFun(theta, d_pm, funname);
Cmm=CovFun(theta, d_mm, funname);
Cnm=CovFun(theta, d_nm, funname);
Cnp=CovFun(theta, d_np, funname);


Cp_s=Cpm*(Cmm\Cnm');
Cp=(Cp_s+(Cnp'-Cp_s).*ID_B)';


