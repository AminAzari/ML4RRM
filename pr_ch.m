function [y,Pout,Poutd]=pr_ch(rrb,rrb_d,Rf, SINRth_ur,Pt_ur,Out_ur_th,NP,BW,PLur)
                       
Pout=1;
% PLZ=10^(-(133+38.3*log10(Rf/1000) )/10);
PLZ=10^(-(PLur )/10);
% NP=10^(-174/10);
% BW=180000;
Gth1=(NP*BW*SINRth_ur)/Pt_ur;
for i=1:rrb_d 
    Pout=Pout*(1-exp(-1*Gth1/PLZ)); 
end
Poutd=Pout;
 
for j=1:length(rrb) 
    if(isnan(rrb(j,2)))
        Pout=Pout*(1-exp(-1*Gth1/PLZ));
    else
        IP=rrb(j,2);
        Gth2=((NP*BW+IP)*SINRth_ur)/Pt_ur;
        Pout=Pout*(1-exp(-1*Gth2/PLZ));
    end 
end



y=1;
if(Pout>Out_ur_th)
    y=-1;
end


% + 10 ? log 10 (X) + Y

% X = exprnd (1, W, 1) ; % make the fast-fading part
% Y = normrnd (0, s, 1, 1); % make lognormal shadowing part
%
% s = 11;
