function [RMS,E] = bundle_adjustment(beta,p)
%beta would have tformto1.T(i)
for n=2:4
    for m=n+1:4
        e=0;
        for i=1:4
            e=e+(norm(transformPointsForward(beta(n-1),p{n,m})-...
                transformPointsForward(beta(m),p{m,n})))^2;
            E(n,m)=sqrt(0.25*e);
        end
    end
end

% Getting RMS:
%initialise RMS = 0 (counter)
RMS=0;
for n=1:4
    for m=n+1:4
        RMS=sqrt((1/6)*sum(sum(E.^2)));
    end
end 