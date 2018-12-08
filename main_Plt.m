%IHN
clc
clear all
close all
pc=0;
for PL=70:5:120% 75:5:125
    pc=pc+1;
    for i=1:40
        X=['data',num2str(PL),num2str(i),'.mat'];
        load(X)
        UDRt(i,1:length(UDR))=UDR(i,:);
        UOUt(i,1:length(UDR))=UOU(i,:);
    end
    
    figure(1)
%     subplot(2,1,1);
    plot([1:0.1:7],nanmean(UDRt))
    A(pc,:)=nanmean(UDRt);
    hold on
%     subplot(2,1,2);
%     plot([1:0.1:7],nanmean(UOUt))
%     B(pc,:)=nanmean(UOUt);
%     hold on
end
figure(2)
% subplot(2,1,1);
plot([1:0.1:7],mean(A))

% subplot(2,1,2);
% plot([1:0.1:7],mean(B))
