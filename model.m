%% This function is used as an observer to give the next state and the next reward using the current state and action
function [next_state,r,ep,rrb,Rfc,SR,OP,OPs,SRd,OPd] =...
    model(xi,Rfc,state,SINRth_ur,Out_ur_th,rrb_d,rrb,Ur_loc,NP,BW,Pt_em,Pt_ur,state_f,Nem,PLur,k)

ik=state; 
g_ine=1;
% NP=10^(-174/10);
% BW=180000;

% Pt=10^(21/10);
Pg=10^(-(133+38.3*log10(Rfc(rrb(ik,1),1)/1000)+10*log10(exprnd(1,1,1)))/10);
SINR_em=0.9*Pt_em*Pg/(NP*BW);
if(rrb(ik,1)~=Nem)
    rrb(ik,2)=Pt_em*Pg;
    RP=Pt_em*Pg;
    V=g_ine*SINR_em/(1+g_ine*SINR_em);
    rrb(ik,3)=max(0,BW*(log2(1+g_ine*SINR_em)-sqrt(V/1)*qfuncinv(0.0001)));
end
[risk,Pout,Poutd]=pr_ch(rrb,rrb_d,Ur_loc,SINRth_ur,Pt_ur,Out_ur_th,NP,BW,PLur);
 
r=(1-xi*Pout/Out_ur_th)*BW*log2(1+g_ine*SINR_em);
Rfc(rrb(ik,1),3)=max(0,Rfc(rrb(ik,1),3)-BW*log2(1+g_ine*SINR_em));
SR=nan;
OPs=nan; 
OP=Pout;
if(state==1 || sum(Rfc(1:Nem-1,3))==0 || and(risk==-1,xi~=0))
    next_state=state_f; 
    ep=1;
    OP=Pout;
 OPd=Poutd;
    if(risk~=-1)
        SR=nansum(rrb(:,3)); 
        OPs=Pout;
        SRd=nansum(rrb(:,3))*(length(rrb(:,3))/sum(rrb(:,3)>0));
    end
    %     sum(Rfc(1:Nem-1,3))
else
    next_state = state-1;
    ep=0;
end
if(rrb(ik,1)==Nem)
    r=0; 
end

