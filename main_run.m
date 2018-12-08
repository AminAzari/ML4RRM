%IHN
% clc
% clear all
% close all
function main_run(PLur,itn)
% PLur=110;
for it=itn
    Uc=0;
    XT=['data',num2str(PLur),num2str(it),'T.mat'];
    X=['data',num2str(PLur),num2str(it),'.mat'];
    if (exist(XT) == 2 || exist(X) == 2) %#ok<*EXIST>
        break
    else
        save(XT,'Uc')
    end
    for U=1:0.1:7
        Uc=Uc+1;
        % learning parameters
        gamma = 0;    % discount factor  % TODO : we need learning rate schedule
        alpha = 0.1;    % learning rate    % TODO : we need exploration rate schedule
        epsilon = 0.85;  % exploration probability (1-epsilon = exploit / epsilon = explore)
        % states
        rrb_s=4;
        rrb_d=1;
        state = [1:rrb_s];
        state_f=state(end);
        % actions
        % AC=[408   577   707   816   912];
        %         AC=[0.7071    1.0000    1.2247    1.4142    1.5811    1.7321    1.8708]*1000;
        AC=[0.9045    1.2792    1.5667    1.8091    2.0226    2.2156    2.3932    2.5584    2.7136    2.8604]*1000;
        action_all = [1:length(AC)+2];
        % initial Q matrix
        Q = zeros(length(state),length(action_all));
        K = 10000;     % maximum number of the iterations
        NP=10^(-174/10);
        BW=180000;
        SINRth_ur=1;
        Out_ur_th=10^(-U);
        Pt_ur=(10^((23)/10))/5;
        Pt_em=10^(21/10);
        state_idx = state(end);  % the initial state to begin from
        
        demU=200000;
        ep=1;
        Ur_loc=3001;
        Cel_ro=3000;
        Cel_ri=50;
        Nem=100;
        AA=0;
        SR=nan;kc=1;OP=nan;OPs=nan;
        %% the main loop of the algorithm
        for k = 1:K
            if(k>2500 && k<5000)
                epsilon=0.95;
            elseif(k>5000)
                epsilon=1;
            end
            
            xi=1;% variable for search in feasible rooths
            
            if (ep==1)
                DR(kc,1)=SR;
                OU(kc,1)=OP;
                OUs(kc,1)=OPs;
                %         DR(kc,2)=k;
                kc=kc+1;
                Rfc=zeros(Nem,3);
                rrb=nan(rrb_s,3);
                %         nodes=datasample(action_all,3);
                Rfc(1:Nem,1)=Cel_ri+(Cel_ro-Cel_ri)*sqrt(rand(1,Nem));
                Rfc(Nem,1)=Ur_loc;
                for ij=1:Nem
                    Rfc(ij,2)=funC2(Rfc(ij,1),[AC,Cel_ro]);
                end
                %                 -(133+38.3*log10(Ur_loc/1000))
                Pg=10^(-(133+38.3*log10(Cel_ro/1000))/10);
                SINR_em=0.9*Pt_em*Pg/(NP*BW);
                g_ine=0.5;
                Rfc(:,3)=demU*BW*log2(1+g_ine*SINR_em)*ones(Nem,1);
                Rfc(Nem,3)=1e100;
                
            end
            
            %     disp(['iteration: ' num2str(k)]);
            r=rand; % get 1 uniform random number
            x=sum(r>=cumsum([0, 1-epsilon, epsilon])); % check it to be in which probability area
            % choose either explore or exploit
            if x == 2   % exploit
                Qt=Q(state_idx,:);
                wc=1;
                while(wc)
                    [~,umax]=max(Qt);
                    for id=1:Nem
                        if(Rfc(id,3)>0)
                            if(Rfc(id,2)==umax)
                                current_action = umax;
                                current_node = id;
                                wc=0;
                            end
                        end
                    end
                    if(wc==1)
                        Qt(umax)=-10000000000;
                    end
                end
            else        % explore
                current_node = datasample([1:Nem],1);
                current_action = Rfc(current_node,2);
            end
            %             x
            %             Qt
            %             current_action
            %             current_node
            rrb(state_idx,1)=current_node;
            %             if(k>2000)
            %                 if(state_idx==2)
            %                     2
            %                 end
            %             end
            %%action_idx = find(action_all==current_action); % id of the chosen action
            % observe the next state and next reward ** there is no reward matrix
            action_idx=current_action;
            
            [next_state,next_reward,ep,rrb,Rfc,SR,OP,OPs] = model(xi,Rfc,state_idx,SINRth_ur,Out_ur_th,rrb_d,rrb,...
                Ur_loc,NP,BW,Pt_em,Pt_ur,state_f,Nem,PLur,k);
            
            next_state_idx=next_state;
            
            %%next_state_idx = find(state==next_state);  % id of the next state
            % print the results in each iteration
            %     disp(['current state : ' num2str(state(state_idx)) ' next state : ' num2str(state(next_state_idx)) ' taken action : ' num2str(action_idx)]);
            %     disp([' dem : ' num2str(sum(Rfc(:,3))),' epis : ' num2str(ep),' next reward : ' num2str(next_reward)]);
            % update the Q matrix using the Q-learning rule
            Q(state_idx,action_idx) = Q(state_idx,action_idx) + alpha * (next_reward + gamma* max(Q(next_state_idx,:)) - Q(state_idx,action_idx));
            
            state_idx = next_state_idx;
            
            if(mod(k,1000)==0)
                [it,U,k/K]
            end
            [~,I]=max(Q,[],2);                              % finding the max values
            AA1=[];
            for i1=1:length(state)
                AA1=[AA1;action_all(I(i1,1))];
            end
            AA(k,1:length(state))=AA1';%[action_all(I(1,1));action_all(I(2,1));action_all(I(3,1));action_all(I(4,1));action_all(I(5,1))]';
            
            
        end
        UDR(it,Uc)=nanmean(DR(floor(length(DR)/2):end));
        UOU(it,Uc)= nanmean(OU(floor(length(DR)/2):end));
        UOUs(it,Uc)=nanmean(OUs(floor(length(DR)/2):end));
        %     mean(UOU)
    end
    if(mod(it,1)==0)
        X=['data',num2str(PLur),num2str(it),'.mat'];
        XT=['data',num2str(PLur),num2str(it),'T.mat'];
        save(X,'UDR','UOU')
        delete (XT)
    end
end
% % % % display the final Q matrix
% % % disp('Final Q matrix : ');
% % % disp(Q)
% % % [C,I]=max(Q,[],2);                                % finding the max values
% % % disp('Q(optimal):');
% % % disp(C);
% % % disp('Optimal Policy');
% % % disp('*');
% % % % AO=[action_all(I(1,1));action_all(I(2,1));action_all(I(3,1));action_all(I(4,1));action_all(I(5,1));action_all(I(6,1));action_all(I(7,1))];
% % % % AO=[action_all(I(1,1));action_all(I(2,1));action_all(I(3,1));action_all(I(4,1));action_all(I(5,1))];
% % % AO=[];
% % % for i1=1:length(state)
% % %     AO=[AO;action_all(I(i1,1))];
% % % end
% % % disp(AO);
% % % disp('*');

% % X=['data',num2str(PLur),'.mat'];
% % save(X,'UDR','UOU')