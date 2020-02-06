%%%Prediction function of the FSA-Block approach
%%%UCsy: U*inv(Cs)*zs
%%%Cs: Cs{i} is the residual covariance matrix of ith block
%%%Q: an m by m matrix in Woodbury inversion formula: B + t(U)*inv(Cs)*U
%%%U: cross-covariance matrix of the training set and the knot set
%%%Cp: the FSA-Block cross covariance matrix between predictive locations and the
%%%training locations.



function [pred_zp, SPE] = Fun_pred_FSABlock(Blockid, zs, zp, alld_sk, UCsy, Cs, Q, U, Cp) 


          Bsize=max(Blockid);
          cB=hist(Blockid,1:Bsize);
          Bindex=[0 cumsum(cB)]; 
          
          invCsy=zeros(size(zs,1),1);
          invCsU=zeros(size(alld_sk.alld,1), size(alld_sk.alld,2));
          
          

          for i=1:(size(Bindex,2)-1)


                  invCsy((Bindex(i)+1):Bindex(i+1),:)=Cs{i}\zs((Bindex(i)+1):Bindex(i+1),:);

                  invCsU((Bindex(i)+1):Bindex(i+1),:)=Cs{i}\U((Bindex(i)+1):Bindex(i+1),:);



          end
          
          
          
          pred_zp=Cp'*(invCsy-invCsU*(Q\UCsy));
          
          SPE=(zp-pred_zp).^2;
          
          
          