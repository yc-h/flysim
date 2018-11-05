function [TPR,FPR,wth]=ROCCmpCPUtoExpc6(con_mode,seq_mode,k,CCrt,Exp380x380)

%con_mode for non-weighted connection('one') or weighted connection('wgt')
%seg_mode=1 for No/INum,...
%seg_mode=2 for Ni/INum,...
%seg_mode=3 for No/ONum,...
%seg_mode=4 for Ni/ONum,...
%seg_mode=5 for No/(INum+ONum),...
%seg_mode=6 for Ni/(INum+ONum)



fname=sprintf('./ConWgt_n380v1-2-1/ConWgt.txt.d%d.th0',k);

fileID = fopen(fname);
ConFiles = textscan(fileID, '%f %f %f');
fclose(fileID);

i=ConFiles{1};
i=i+1;
i=i';
j=ConFiles{2};
j=j+1;
j=j';
w_=ConFiles{3};
w=ones(1,length(w_));

clear ConFiles;
if strcmp(con_mode,'one')
    C=sparse(i,j,w,380,380);
else
    C=sparse(i,j,w_,380,380);
end

CPU380x380=full(C);


%png=imagesc(CPU380x380);
%map = colormap;
%map(1,:) = [1,1,1];
%colormap(map);
%colorbar;

%clear C;


%fname='CPU380x380';
%fname_arb=sprintf('%s.%d.png',fname,k);
%title(fname);
%saveas(gcf,fname_arb,'png');

[m,n]=size(CPU380x380);


fname=sprintf('./ConWgt_n380v1-2-1/ConWgt.deg.d%d.th0',k);

fileID = fopen(fname);
ConFiles = textscan(fileID, '%f %f %f %f %f %f %f %f');
fclose(fileID);

TDeg=ConFiles{1};
TNum=ConFiles{2};
IDeg=ConFiles{3};
INum=ConFiles{4};
ODeg=ConFiles{5};
ONum=ConFiles{6};
BDeg=ConFiles{7};
BNum=ConFiles{8};


if seq_mode==1
    Num=INum;
elseif seq_mode==2
    Num=INum;
elseif seq_mode==3
    Num=ONum;
elseif seq_mode==4
    Num=ONum;
else
    Num=(INum+ONum);
end

dwth=Num/CCrt;
wth=zeros(length(Num),CCrt);

for i=1:m
    if Num(i)~=0
        k=0:dwth(i):Num(i)-dwth(i);
        wth(i,:)=k;
    end
end


CPUth=zeros(m,n);

TPR=zeros(length(CCrt),1);
FPR=zeros(length(CCrt),1);

for k=1:CCrt
    TP=0;
    FP=0;
    FN=0;
    TN=0;


    if seq_mode==1 || seq_mode==3 || seq_mode==5
        for i=1:m  % No compare to Num
            for j=1:n
                if CPU380x380(i,j)>wth(i,k)
                    CPUth(i,j)=1;
                end
            end
        end
    else
        for j=1:n  % Ni compare to Num
            for i=1:m
                if CPU380x380(i,j)>wth(j,k)
                    CPUth(i,j)=1;
                end
            end
        end
    end



    for i=1:length(CPUth(:))
    
        if CPUth(i)==1 && Exp380x380(i)==1
            TP=TP+1;
        elseif CPUth(i)==1 && Exp380x380(i)==0
            FP=FP+1;
        elseif CPUth(i)==0 && Exp380x380(i)==1
            FN=FN+1;
        else
            TN=TN+1;
        end
        CPUth(i)=0;
    end
    
    TPR(k)=TP/(TP+FN);
    FPR(k)=FP/(FP+TN);
end


end
