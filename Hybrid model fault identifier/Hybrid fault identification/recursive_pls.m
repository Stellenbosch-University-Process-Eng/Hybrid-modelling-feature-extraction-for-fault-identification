function [b1,Q,P,W,T,U,Bcat,Qcat,Yt,Wcat,Pcat] = recursive_pls (X,Y,ncomp,lambda,bsize)

M = floor(size(X,1)/bsize); %Number of blocks which satisfy blocksize

for i = 1:1:M
    
    %[~,~,b,~,P,Q,W] = npls(X(1+bsize*(i-1):bsize*i,:),Y(1+bsize*(i-1):bsize*i,:), ncomp);   
    [~,~,b1,~,P,Q,W] = npls(X(1+bsize*(i-1):bsize*i,:),Y(1+bsize*(i-1):bsize*i,:), ncomp);
    
    B = diag(b1);
    
    
    if i>1
        
    %X0 = [lambda*Pold';P'];
    %Y0 = [lambda*Bold*Qold';B*Q'];
     X0 = [P';lambda*Pold'];
     Y0 = [B*Q';lambda*Bold*Qold'];
    
    [T,U,b1,~,P,Q,W] = npls(X0,Y0,ncomp);
    
    B = diag(b1);
    %B = b1;
        
    end
    
    Bold = B;
    Pold = P;
    Qold = Q;
    
    Bcat(i,:) = b1;
    Qcat(i,:) = Q;
    %Wcat(bsize*(i-1)+1:bsize*i,:) = W;
    Wcat(:,:,i) = W;
    Pcat(:,:,i) = P;
    
    % predictions for the block used in the regression
    
    tscale = xscale(X(1+bsize*(i-1):bsize*i,:),W,P,ncomp);
    
    for j=1:1:size(tscale,1)
    	Ytemp(j,:) = tscale(j,:).*b1.*Q;
    end

    Yt(1+bsize*(i-1):bsize*i,:)= sum(Ytemp,2);
    Yt(1+bsize*(i-1):bsize*i,:)= Yt(1+bsize*(i-1):bsize*i,:)+mean(Y);
    
end
end
function[Tscale] =  xscale (X,W,P,ncomp)

X0 = X - mean(X);

for i = 1:1:ncomp
    
    Tscale(:,i) = X0*W(:,i);
    X0 = X0-Tscale(i)*P(i);
end

end
