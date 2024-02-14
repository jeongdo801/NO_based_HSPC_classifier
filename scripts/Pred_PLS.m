function [ypred, that1] = Pred_PLS(xtrain1, ytrain, xtest1, lv)
[m1,ssq1,p1,q1,w1,t1,u1,b1] = pls(xtrain1,ytrain,lv);
B = zeros(lv,lv);
    for l=1:lv
        B(l,l)=b1(l);that1(:,l) = xtest1*w1(:,l);
    end
        
    for j = 1:lv
        ypred=that1(:,1:j)*B(1:j,1:j)*q1(:,1:j)';
    end
end