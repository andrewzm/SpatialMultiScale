function[Cmat]=CovFun(theta,d,funname)
 %%%% 
switch funname
                
    case 'exp' 
        
        d=d.alld;
        Cmat=theta.alpha.*exp(-d./theta.beta);        
            
          
            
    
    case 'matern'
            
            d=d.alld/theta.phi;
            
            Cmat=(theta.sigma2/gamma(theta.nu)*2^(1-theta.nu)).*d.^(theta.nu)...
             .*besselk(theta.nu, d);
         
            Cmat(d==0)=theta.sigma2;
        
    case 'matern1'
            
            d=d.alld/theta.phi;
            theta.nu = 1;
            
            Cmat=(theta.sigma2/gamma(theta.nu)*2^(1-theta.nu)).*d.^(theta.nu)...
             .*besselk(theta.nu, d);
         
            Cmat(d==0)=theta.sigma2;
            
    
   
end
