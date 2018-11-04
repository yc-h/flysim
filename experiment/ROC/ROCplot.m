clear all;
close all;


TabLen=380;
A=zeros(TabLen,152);
fileID = fopen('../../../Network_parameters/n380/NP380x152.txt');
ConFiles = textscan(fileID, '%f','delimiter',',','EmptyValue',0); %for .cvs file
fclose(fileID);
T=ConFiles{1};
for i=1:TabLen
    for j=1:152
        A(i,j)=T((i-1)*152+j);
    end
end
clear T;

x=zeros(length(find(A)),1);
y=x;
k=0;
for i=1:152
    for j=1:TabLen
        if A(j,i)~=0
            k=k+1;
            x(k)=i;
            y(k)=j;
        end
    end
end
%png=plot(x,y,'.');
%title('NP152xn380');
%saveas(gcf,'NP152xn380.png','png');

Exp380x380=zeros(TabLen);
for i=1:TabLen
    for j=1:152
        if A(i,j)==2
            for m=1:TabLen
                if A(m,j)==1
                    Exp380x380(i,m)=1;
                end
            end
        end
    end
end

b_FPR=cell(50,1);
b_TPR=cell(50,1);
b_wth=cell(50,1);
Max_TPR=zeros(50,1);


CCrt=300;

%seg_mode=1 for No/INum,...
%seg_mode=2 for Ni/INum,...
%seg_mode=3 for No/ONum,...
%seg_mode=4 for Ni/ONum,...
%seg_mode=5 for No/(INum+ONum),...
%seg_mode=6 for Ni/(INum+ONum)
seq_mode=5;

%for i=1:50
parpool(4);
parfor i=1:50
    i
    [a_TPR,a_FPR,a_wth]=ROCCmpCPUtoExpc6('Wgt',seq_mode,i,CCrt,Exp380x380);

    b_FPR{i}=a_FPR;
    b_TPR{i}=a_TPR;
    b_wth{i}=a_wth;
end
delete(gcp);


cla
png=plot(b_FPR{20},b_FPR{20},'Color',[0.6 0.6 0.6],'lineWidth',4);
hold;

plot(b_FPR{1},b_TPR{1},':k',b_FPR{13},b_TPR{13},'-k',b_FPR{20},b_TPR{20},'--k','lineWidth',4);
set(gca,'Fontsize',25);
set(gca,'FontWeight','bold');
set(gca,'lineWidth',3.5);
set(gca,'TickLength',[0.02 0.035]);
ylabel('True postive rate');
xlabel('False postive rate');
legend('Random selection','Distance criterion=1um','Distance criterion=13um','Distance criterion=20um','Location','east');
fname='ROC_space_seq5';
%fname_arb=sprintf('%s.%d.png',fname,i);
fname_arb=[fname '.png'];
%title(fname_arb);

legend boxoff;
box off;

figure_size = get(gcf, 'position');
set(gcf,'PaperPosition',figure_size/20);
print(gcf,'-dpng','-r300', fname_arb);

%{
set(gca,'xlim',[0.04 0.1]);
set(gca,'ylim',[0.65 0.75]);
figure_size = get(gcf, 'position');
set(gcf,'PaperPosition',figure_size/60);
print(gcf,'-dpng','-r300','D1-13-20_local');
%}


