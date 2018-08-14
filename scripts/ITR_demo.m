% Examples of ITR calculation
t=2;%time per selection
p=0:0.1:1;%accuracy
for k=1:length(p)
itr32(k)=itr(32,p(k),t);
itr16(k)=itr(16,p(k),t);
itr4(k)=itr(4,p(k),t);
end

figure;
plot(p*100,itr32,'b');
hold on
plot(p*100,itr16,'k');
plot(p*100,itr4,'r');


xlabel('Accuracy (%)')
ylabel('ITR (bpm)')
legend('n=32 (5 bits)','n=16 (4 bits)','n=4 (2 bits)');

h=gca;
set(h,'YTick',10:10:150)
title('Speed: 2 seconds per selection (30 selections per minute)')
box off