%||AUM||
function [false_pos_rate,true_pos_rate]=roc_pval(pval,grnd_truth,stps)
pval=pval(:);grnd_truth=grnd_truth(:);
cnt=1;
ST=linspace(0,1,stps);false_pos_rate=[];true_pos_rate=[];
for thr=ST
    false_pos_rate(cnt)=sum((pval<=thr)&(~grnd_truth))/sum(grnd_truth);
    true_pos_rate(cnt)=sum((pval<=thr)&(grnd_truth))/sum(grnd_truth);
    cnt=cnt+1;
end

