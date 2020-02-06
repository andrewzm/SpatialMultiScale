function D = distance(x1, x2, type)

n1=size(x1,1);
n2=size(x2,1);
D=zeros(n1,n2);

switch type
    
    case 'euclidean'
        
       for j=1:n2
           
           D(:,j)= sqrt(sum((x1-repmat(x2(j,:), n1,1)).^2,2));
           
       end
        
        
        
end