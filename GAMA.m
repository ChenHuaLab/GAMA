function M=GAMA(obsdata_ancestor,obsdata_descendant,position,rr,binSize,varargin)
%����ʶ�������Ƭ�Ρ�
%obsdata_ancestor:���ȹ۲����С������ȸ�����cell��ɣ�ÿ��cell����һ��
%���ȵ����ݣ���һ��������m�е�cell���д�����壬�д���Ⱦɫ�塣
%obsdata_descendant:����۲����У���һ��n��m�е�cell���д�����壬�д���Ⱦɫ�塣
%disdata:m����������ɵ�cell��Ӧÿ����������ӦȾɫ���snp���롣
%r: ÿ��ÿ��snp����ϵ����
%
%����һ��strcture,����ÿ��Ⱦɫ��������Լ���Ӧ���ƶ�Ƭ�Ρ����ز�����ֵ�����һ�ͼ��
%20161223 �޸ģ�����waitbar
% r�������Ŵ������ֵ�ļ���ע�⻻���ÿ���ľ��룬ע���ֵ��ע�������Ŵ��������1 �������
% Cval, ÿ��λ�������ʶԱȶ�֮��
% 20191104�޸ģ���hidstate����е�һ�м�������λ����Ϣ�������۵�ǰ�Ľ����

global Fmat
global Bmat
global Tmat
global obslik
global snp_distance
global r
global Eijl
global Omat
global TTmat
global PP
global difnum
global binLen
global numdescendant
global state_distrib
global admix_generation

binLen=binSize;

difnum=zeros(2^binSize);
for i=1:2^binSize
    for j=1:2^binSize
        b1=double(dec2bin(i-1,binSize));
        b2=double(dec2bin(j-1,binSize));
        difnum(i,j)=sum(b1~=b2);
    end
end



if nargin>6||nargin<5
    error('The number of input parameters between 5-6');
end

%
if nargin==5
    maxitersteps=20; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tolvalue=1*10^-2;
else
    options=varargin{1};
    maxitersteps=options.maxitersteps;
    tolvalue=options.tolvalue;
end

% h=waitbar(0,'competed 0%');

if size(rr,1)>size(rr,2)
    rr=rr';
end

r=rr;
%*************************Ⱦɫ�����Լ��۲����еĻ�����Ϣ********************

if size(obsdata_descendant,1)>size(obsdata_descendant,2)  % 20180330
    obsdata_descendant=obsdata_descendant';
end

if size(position,2)>size(position,1)  % positionת��Ϊ������
    position=position';
end

disdata=position(2:end)-position(1:end-1);

if ~isscalar(r)
    disr=r(2:end)-r(1:end-1);
end

if size(disdata,1)>size(disdata,2)
    disdata=disdata';
end


% ���ݸ�ʽת��
A1=obsdata_ancestor{1}(1,:);

for i=1:length(obsdata_ancestor)
    temp=obsdata_ancestor{i};
    tempind=temp==repmat(A1,size(temp,1),1);
    obsdata_ancestor{i}(tempind)=1;
    obsdata_ancestor{i}(~tempind)=0;
end

tempind=obsdata_descendant==repmat(A1,size(obsdata_descendant,1),1);
obsdata_descendant(tempind)=1;
obsdata_descendant(~tempind)=0;

% �ض϶��������
numMod=mod(size(obsdata_descendant,2),binSize);
numsnp0=floor(size(obsdata_descendant,2)/binSize);
obsdata_descendant=obsdata_descendant(:,1:end-numMod);
%���۲������۵�
temp=zeros(size(obsdata_descendant,1),numsnp0);
for i=1:size(temp,2)
    temp(:,i)=obsdata_descendant(:,binSize*(i-1)+1:binSize*i)*(2.^(binSize-1:-1:0))';
end
obsdata_descendant=temp;
% �����������۵�
disdata=disdata(1:end-numMod);
if ~isscalar(r)
    disr=disr(1:end-numMod);
end

temp=zeros(1,numsnp0-1);
for i=1:size(temp,2)
    temp(i)=sum(disdata(binSize*(i-1)+1:binSize*i-1))/2+sum(disdata(binSize*i+1:binSize*(i+1)-1))/2+disdata(binSize*i);
end
disdata=temp;

if ~isscalar(r)
    temp=zeros(1,numsnp0-1);
    for i=1:size(temp,2)
        temp(i)=sum(disr(binSize*(i-1)+1:binSize*i-1))/2+sum(disr(binSize*i+1:binSize*(i+1)-1))/2+disr(binSize*i);
    end
    disr=temp;
end



for i=1:length(obsdata_ancestor)
    obsdata_ancestor{i}=obsdata_ancestor{i}(:,1:end-numMod);
end
numancestor=length(obsdata_ancestor);
% �����������۵�
for k=1:numancestor
    temp=zeros(size(obsdata_ancestor{k},1),numsnp0);
    tempancestor=obsdata_ancestor{k};
    for i=1:size(temp,2)
        temp(:,i)=tempancestor(:,binSize*(i-1)+1:binSize*i)*(2.^(binSize-1:-1:0))';
    end
    obsdata_ancestor{k}=temp;
end



[numdescendant,numsnp0]=size(obsdata_descendant);
snp_distance=repmat([disdata,Inf],1,numdescendant);  % 20170213
snp_distance=snp_distance(1:end-1);
obsdata_descendant=obsdata_descendant';
obslik=obsdata_descendant(:)'; % 20170213

if ~isscalar(r)
    disr=repmat([disr,1],1,numdescendant);  % 20200103
    disr=disr(1:end-1);
    r=disr;
end
% alphaset=0:0.02:0.2;

numsample=0;
for i=1:length(obsdata_ancestor)
    numsample=numsample+size(obsdata_ancestor{i},1);
end

countSample=cell(1,numancestor);
for i=1:numancestor
    temp=zeros(2^binSize,numsnp0);
    for j=1:2^binSize   %�������Щ�п��ܵ�haplotype����
        temp(j,:)=sum(obsdata_ancestor{i}==j-1);  
    end
    countSample{i}=temp;
end

allSample=zeros(numsample,numsnp0);
count=1;
for j=1:length(obsdata_ancestor)
    allSample(count:count+size(obsdata_ancestor{j},1)-1,:)=obsdata_ancestor{j};
    count=count+size(obsdata_ancestor{j},1);
end

for i=1:numsnp0
    tempval=unique(allSample(:,i));
    for j=1:length(obsdata_ancestor)
%         countSample{j}(tempval+1,i)=countSample{j}(tempval+1,i)+0.5;   % 20180719 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        countSample{j}(tempval+1,i)=countSample{j}(tempval+1,i)+size(obsdata_ancestor{j},1)/numsample;  % 20180720
    end
end
        
   

%����ÿһ������Ⱥ���P�������һ��cell����
PP=cell(1,numancestor);  
for i=1:numancestor
    temp=zeros(2^binSize,numsnp0);
    
    tempnum=sum(countSample{i});
    for j=1:2^binSize   %�������Щ�п��ܵ�haplotype����
        
        temp(j,:)=0.95*countSample{i}(j,:)./tempnum+0.05/2^binSize;
%         temp(j,:)=0.5*countSample{i}(j,:)./tempnum+0.5/2^binSize;
%         temp(j,:)=countSample{i}(j,:)./tempnum;
        
    end
    PP{i}=temp;
end


% �����P�����ǳ�ʼֵ����������Ʋ⵽�������������и���
P=cell(1,numancestor); 

%���ǵ��������⣬��P����ת����
for i=1:numancestor
    temp=0.05.^difnum*PP{i};
    P{i}=temp./repmat(sum(temp),size(temp,1),1)+eps;
end

% for i=1:numancestor
%     P{i}(P{i}<=0.0001)=0.0001;
%     P{i}(P{i}>=0.9999)=0.9999;
% end





for i=1:numancestor
    P{i}=repmat(P{i},1,numdescendant); %20170213 , % 20180619,P��һ��cell��ÿ��cell�洢һ��Ⱥ���Ƶ�ʾ���
end




numsnp=numsnp0*numdescendant;




%***************************����������*************************************

%����������ʼֵpar
a=ones(1,numancestor);  %%%%%%%%%%%%%%%%%%%%%%%%
a=a/sum(a);
state_distrib=a;
admix_generation=20;

% iter_EM=20;
parstore=zeros(maxitersteps+1,numancestor+2);
parstore(1,:)=[state_distrib,admix_generation,0.05];


Fmat=fwdmat_power(state_distrib,admix_generation,snp_distance,r,obslik,P);
Bmat=backmat_power(state_distrib,admix_generation,snp_distance,r,obslik,P);
Tmat=transmat(state_distrib,snp_distance,r,admix_generation);

%��Ϊ���飬������15�Σ�����������׼��Ҫ�������ơ�

for j=1:maxitersteps

    %Ϊ�Ż���������20160914
    F=Fmat{1};
    B=Bmat{1};
    Omat=obsmat(obslik,P);   



Eijl=zeros(numancestor,(numsnp-1)*numancestor);
for l=1:numsnp-1
    temp=(F(:,l)*(Omat(:,l+1).*B(:,l+1))').*Tmat{l};
    temp=temp-diag(diag(temp));
    Eijl(:,(l-1)*numancestor+1:l*numancestor)=temp;
end

TTmat=zeros((numsnp-1)*numancestor,numancestor);
for ii=1:numsnp-1
    TTmat((ii-1)*numancestor+1:ii*numancestor,:)=Tmat{ii};
end

   
    
    par0=parstore(j,end);

    if isscalar(r) % r��һ���������һ���Ŵ���������

        parnow=fmincon(@E_fun_ga_power,par0,-1,0,[],[],0,0.5);

    else
        parnow=fmincon(@E_fun_ga_power_gd2,par0,-1,0,[],[],0,0.5);

    end

    parstore(j+1,end)=parnow;  % ���һ��������ʹ����ţ�ٷ����
    
    %******************************   veterbi    *****************************

    Vmat=zeros(numancestor,numsnp);
    Vmat(:,1)=state_distrib;
    indmat=zeros(numancestor,numsnp);
    % postp2=zeros(numancestor,numsnp);

    Omat=obsmat(obslik,P);

    for l=2:numsnp
        for jj=1:numancestor
            Vmat(jj,l)=max(Vmat(:,l-1).*Tmat{l-1}(:,jj)*Omat(jj,l));
        end

        Vmat(:,l)=Vmat(:,l)/min(Vmat(:,l)); %20160829��ֵ���㾫�����⣡

        for k=1:numancestor
            temp=Vmat(:,l-1).*Tmat{l-1}(:,k);
            [~,indmat(k,l)]=max(temp);
        end

    end

    hid_state=zeros(1,numsnp);
    [~,hid_state(numsnp)]=max(Vmat(:,numsnp));

    for l=numsnp-1:-1:1
        hid_state(l)=indmat(hid_state(l+1),l+1);
    end
    
    for aiter=1:numancestor
        state_distrib(aiter)=mean(hid_state==aiter);
    end
    
    state_distrib=(state_distrib+eps)/sum(state_distrib+eps); % 20180918 revised, ��ֹ��E_fun�г���log0

    % chrLen=(position(end)-position(1))/10^8;
%     chrLen=sum(1-exp(-r*(position(2:end)-position(1:end-1))))*numdescendant;
%     temphid_state=reshape(hid_state,numsnp0,numdescendant)';
%     tempdif=temphid_state(:,1:end-1)~=temphid_state(:,2:end);
%     tempnum=sum(tempdif(:));
%     tempt=sum(state_distrib.*(1-state_distrib));
%     admix_generation=tempnum/(tempt*chrLen);

%     admix_generation=100;
%     admix_generation=1/(min(state_distrib)*chrLen/(0.5*tempnum))/(1-min(state_distrib)); 
    

%     hidstate=reshape(repmat(hid_state,binSize,1),binSize*numsnp0,numdescendant);  % ��binSizeչ������position��Ӧ
%     if size(hidstate,1)>10000
%         interval=floor(size(hidstate,1)/5000);
%         hidstate=hidstate(1:interval:end,:);
%         pos=position(1:interval:end);
%     else
%         pos=position;
%     end  % Ϊ�˼ӿ������ٶȣ�������ȡ��
%     prime=primes(10000);
%     prime=[1,prime(1:numancestor-1)];
%     for i=numancestor:-1:1   %����������Ϊ����������ظ�
%         hidstate(hidstate==i)=prime(i);
%     end
%     
%     dis=abs(repmat(pos,1,length(pos))-repmat(pos',length(pos),1));
%     dis=round(dis/10^6*5);% ��1cMΪ��λ
%     curvedata=cell(numancestor,numancestor);
%     disset=unique(dis(:));
%     disnum=zeros(size(disset));
%     for i=1:length(disnum)
%         disnum(i)=sum(dis(:)==disset(i));
%     end
%     
%     
%     for i=1:numancestor
%         for j2=1:numancestor
%             curvedata{i,j2}=zeros(size(disset));
%         end
%     end
% 
%     for i=1:numdescendant
%         temphidstate=hidstate(:,i);
%         A=temphidstate*temphidstate';
%         
%         for j2=1:numancestor
%             for k=1:numancestor
%                 temp=zeros(size(disset));
%                 multistate=prime(j2)*prime(k);
%                 ind1=A==multistate;
%                 for h=1:length(disset)
%                     temp(h)=sum(sum(ind1&(dis==disset(h))));
%                 end
%                 curvedata{j2,k}=curvedata{j2,k}+temp./disnum;
%             end
%         end
%          
%     end
% %     figure
% %     plot(disset,curvedata{1,1}/numdescendant)
%     curvefitdata=[];
%     curvefitdata(:,1)=disset(1:20)/500;
%     curvefitdata(:,2)=curvedata{1,1}(1:20)/numdescendant;
% 
%     lambda=fmincon(@(d)sum((d(2).*(1-d(2)).*exp(-d(1).*curvefitdata(:,1))+d(2).^2-curvefitdata(:,2)).^2),[2,0.1]);
%     admix_generation=lambda(1);
    
%     curvefitdata=[];
%     curvefitdata(:,1)=disset(1:end)/1000;
%     curvefitdata(:,2)=curvedata{2,2}(1:end)/numdescendant;
% 
%     lambda=fmincon(@(d)sum((d(2).*(1-d(2)).*exp(-d(1).*curvefitdata(:,1))+(1-d(2)).^2-curvefitd8ata(:,2)).^2),[2,0.1]);
%     admix_generation=lambda(1);
%     
%     
%     
%     curvefitdata=[];
%     curvefitdata(:,1)=disset(1:end)/1000;
%     curvefitdata(:,2)=curvedata{1,2}(1:end)/numdescendant;
% 
%     lambda=fmincon(@(d)sum((-d(2).*(1-d(2)).*exp(-d(1).*curvefitdata(:,1))+d(2)*(1-d(2))-curvefitdata(:,2)).^2),[2,0.1]);
%     admix_generation=lambda(1);
%     
    
   
    
    
    hidstate=reshape(repmat(hid_state,binSize,1),binSize*numsnp0,numdescendant);  % ��binSizeչ������position��Ӧ
    
    %��¼ÿ��chunk�Ĵ�С��λ��
    chunk=cell(numancestor,numdescendant); %��¼�����ӽ�����������в�ͬ��Դchunk�����ݣ�����chunk��ʼ����ֹ��λ�á�����
    chunkUnit=10^4;
    for iterd=1:numdescendant
        temp=hidstate(:,iterd);
        for itera=1:numancestor
            tempind=temp==itera;
            indS=find(tempind(1:end-1)==0&tempind(2:end)==1)+1;
            indE=find(tempind(1:end-1)==1&tempind(2:end)==0);
            if tempind(1)==1
                indS0=[1;indS];
                indS=indS0; %ʵ�ڲ��뿴��С���
            end
            
            if tempind(end)==1
                indE0=[indE;size(hidstate,1)];
                indE=indE0; %ʵ�ڲ��뿴��С���
            end
            
            % ע�⣡  ���������chunk�����п�����0
            tempchunk=[round(position(indS)/chunkUnit),round(position(indE)/chunkUnit),round((position(indE)-position(indS))/chunkUnit),round((position(indE)+position(indS))/2/chunkUnit)];
            chunk{itera,iterd}=tempchunk;
        end
    end
    
    chunksize=[];
    for i=1:size(chunk,2)
        chunksize=[chunksize;chunk{1,i}(:,3)];
    end
    lenchr=round((position(end)-position(1))/chunkUnit);
    MaxChunk=max(chunksize);
%     curve=zeros(1,lenchr);
    curve=zeros(1,lenchr+1); % 20191224
    for i=1:length(chunksize)
        tempchunksize=chunksize(i);
        for jj=1:tempchunksize-1
            curve(jj)=curve(jj)+tempchunksize-jj;
        end
    end
    
    for i=1:size(chunk,2)
        chunks=chunk{1,i};
        for jj=1:size(chunks,1)-1
            for k=jj+1:size(chunks,1)
                chunk1=chunks(jj,:);
                chunk2=chunks(k,:);
                dchunk=chunk2(1)-chunk1(2);
                L1=chunk1(3);
                L2=chunk2(3);
                curve(dchunk+1:dchunk+min(L1,L2))=curve(dchunk+1:dchunk+min(L1,L2))+(1:min(L1,L2));
                curve(dchunk+min(L1,L2)+1:dchunk+max(L1,L2))=curve(dchunk+min(L1,L2)+1:dchunk+max(L1,L2))+min(L1,L2);
                curve(dchunk+max(L1,L2)+1:L1+L2+dchunk-1)=curve(dchunk+max(L1,L2)+1:L1+L2+dchunk-1)+(min(L1,L2)-1:-1:1);
%                 curve(dchunk+max(L1,L2)+1:min(L1+L2+dchunk-1,length(curve)))=curve(dchunk+max(L1,L2)+1:L1+L2+dchunk-1)+(min(L1,L2)-1:-1:1); % 20191224
            end
        end
    end

    curve=curve(1:MaxChunk-1);
    curve=curve./(lenchr-1:-1:lenchr-MaxChunk+1)/numdescendant;
    curvefitdata=[];
    curvefitdata(:,1)=(1:MaxChunk-1)/(10^8/chunkUnit);
    curvefitdata(:,2)=curve;
    for k=2:length(curve)
        if (curve(k)-min(curve(1:k-1)))/min(curve(1:k-1))>0.01
            cutpoint=k;
            break
        end
        cutpoint=length(curve);
    end
%     
    %��curveƽ��������alpha^2�Ժ󣩣�Ч�����ܸ���
     curvefitdata0=curvefitdata(1:cutpoint,:);   
%     ind=find(curve(1:cutpoint-1)==min(curve(1:cutpoint-1)));

%     curvefitdata0=curvefitdata(1:end,:);
%     if admix_generation<20
%         curvefitdata0=curvefitdata(1:round(size(curvefitdata,1)/1),:);
%   
%     elseif admix_generation<50
%         curvefitdata0=curvefitdata(1:round(size(curvefitdata,1)/10),:);
%     else
%        curvefitdata0=curvefitdata(1:round(0.05/(chunkUnit/10^8)),:);      
%     end
    
%     ind=find(curvefitdata(:,2)<=state_distrib(1)^2,1);
%     curvefitdata0=curvefitdata(1:ind,:);
    
%     curvefitdata0=curvefitdata(1:1000,:);
%     lambda=fmincon(@(d)sum((d(2).*(1-d(2)).*exp(-d(1).*curvefitdata0(:,1))+d(2)^2-curvefitdata0(:,2)).^2),[2,0.1]);
    lambda=fmincon(@(d)sum((d(2).*(1-d(2)).*exp(-d(1).*curvefitdata0(:,1))+d(2)^2-curvefitdata0(:,2)).^2./(1:length(curvefitdata0))'),[20,0.5]);
    
%     lambda=fmincon(@(d)-sum(log((d(2).*(1-d(2)).*exp(-d(1).*curvefitdata0(:,1))+d(2)^2).^curvefitdata0(:,2))),[2,0.1]);
    admix_generation=lambda(1);
%     admix_generation=20;
%     admix_generation=0.5;
%     curve=zeros(10000000,1);
%     p=1;
%     for iterd=1:numdescendant
%         tempchunk=chunk{1,iterd};
%         for jj=1:size(tempchunk,1)
%             curve(p)=tempchunk(jj,3);
%             p=p+1;
% %             for kk=jj+1:size(tempchunk,1)
% %                 curve(p,1)=tempchunk(kk,4)-tempchunk(jj,4);
% %                 curve(p,2)=tempchunk(kk,3)*tempchunk(jj,3);
% %                 p=p+1;
% %             end
%         end
%     end
%     curve=curve(1:p-1,:);
%     
%     % ����
%     curvefitdata=zeros(1000,1);
%     curvefitdata(:,1)=(0.1:0.1:100)';
%     alpha=state_distrib(1);  %%%%%%%%%%%%%% 
%     for jj=1:1000
%         tempcurve=curve>curvefitdata(jj,1);
%         curvefitdata(jj,2)=sum(curve(tempcurve)-curvefitdata(jj,1))/((position(end)-position(1))*alpha/10^6)/numdescendant;
%     end
%     
% %     curve(:,1)=round(curve(:,1)/0.1)*0.1; % 0.1�����񣬶�Ӧrecent grid
% %     
% %     alldis=sort(unique(curve(:,1)));
% %     alldis=alldis(alldis<100);
% %     curvefitdata=zeros(length(alldis),2);
% %     curvefitdata(:,1)=alldis;
% %     for jj=1:length(alldis)
% %         curvefitdata(jj,2)=sum(curve(curve(:,1)==alldis(jj),2));
% %     end
%     
% %     lambda=fmincon(@(d)sum((d(3).*exp(-d(1).*curvefitdata(:,1)/100)+d(2)-curvefitdata(:,2)).^2),[1,10,1]);
% %     admix_generation=lambda(1);
% 
% %     alpha=state_distrib(1);
% %     lambda=fmincon(@(d)sum((alpha.*(1-alpha).*exp(-d(1).*curvefitdata(:,1)/100)+alpha.^2-curvefitdata(:,2)).^2),1);
% %     admix_generation=lambda(1);
%     
%     lambda=fmincon(@(d)sum((d(2).*(1-d(2)).*exp(-d(1).*curvefitdata(:,1)/100)+d(2).^2-curvefitdata(:,2)).^2),[2,0.1]);
%     admix_generation=lambda(1);
% 
%     plot(0.1:0.1:100,lambda(2).*(1-lambda(2)).*exp(-lambda(1).*curvefitdata(:,1)/100)+lambda(2).^2)
%     hold on
%     plot(0.1:0.1:100,curvefitdata(:,2),'r')

    parstore(j+1,1:numancestor)=state_distrib;
    parstore(j+1,numancestor+1)=admix_generation;
    
    if tempnum==0
        warning('No admixture')
        break 
    end
    
    if all(abs(parstore(j+1,1:end-1)-parstore(j,1:end-1))./parstore(j,1:end-1)<=tolvalue) %       ????????????????????????????? 20190124���һ����������Ϊ��ֹ����
        break
    end 

    
    %���²�����������һ��EM����
    alpha=parnow;
%     alpha=0.05;
    P=cell(1,numancestor);
    for i=1:numancestor
        temp=alpha.^difnum*PP{i};
        P{i}=temp./repmat(sum(temp),size(temp,1),1);
    end

%     for i=1:numancestor
%         P{i}(P{i}<=0.0001)=0.0001;
%         P{i}(P{i}>=0.9999)=0.9999;
%     end


    for i=1:numancestor
        P{i}=repmat(P{i},1,numdescendant); %20170213 , % 20180619,P��һ��cell��ÿ��cell�洢һ��Ⱥ���Ƶ�ʾ���
    end

    Fmat=fwdmat_power(state_distrib,admix_generation,snp_distance,r,obslik,P);
    Bmat=backmat_power(state_distrib,admix_generation,snp_distance,r,obslik,P);
    Tmat=transmat(state_distrib,snp_distance,r,admix_generation);

%     waitbar(j/maxitersteps,h,['completed ',num2str(j*100/maxitersteps),'%']);




% waitbar(1,h,'completed 100%');
% pause(2)
% close(h)



end

parstore=parstore(~all(parstore==0,2),:);%remove the all 0 row


%%%%%%%%%%%%%%%%%%%%%%%% 2022-08-14 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postp=Vmat./repmat(sum(Vmat),size(Vmat,1),1); % 20220814 postp�޸�
postp_norm=zeros(size(Vmat));
postp_raw=zeros(size(Vmat));

for l=2:numsnp-1
    tempProb=Tmat{l-1}(hid_state(l-1),:)'.*Omat(:,l).*Tmat{l}(:,hid_state(l+1)); 
    postp_raw(:,l)=tempProb;
    postp_norm(:,l)=tempProb/sum(tempProb);
end
% ��һ��
tempProb=state_distrib'.*Omat(:,1).*Tmat{1}(:,hid_state(2)); 
postp_raw(:,1)=tempProb;
postp_norm(:,1)=tempProb/sum(tempProb);
% ���һ��
tempProb=Tmat{end}(hid_state(end-1),:)'.*Omat(:,end); 
postp_raw(:,end)=tempProb;
postp_norm(:,end)=tempProb/sum(tempProb);
%%%%%%%%%%%%%%%%%%%%%%%% 2022-08-14 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    



% 20180629. ʹ�ýṹ��洢�����
hidstate=reshape(hid_state,numsnp0,numdescendant)'; %20191104
hidstate=reshape(repmat(hidstate,binSize,1),size(hidstate,1),binSize*size(hidstate,2)); %�����۵�
hidstate=[position(1:size(hidstate,2))';hidstate]; %������λ�÷���hidstate�ĵ�һ�С�20191104

% 20220814
temp=zeros(numdescendant,numsnp0*binSize,numancestor);
for i=1:numancestor
    pp=Omat(i,:);
    pp=reshape(pp,numsnp0,numdescendant)';
    pp=reshape(repmat(pp,binSize,1),size(pp,1),binSize*size(pp,2));
    temp(:,:,i)=pp;
end
observeProb=temp;

temp=zeros(numdescendant,numsnp0*binSize,numancestor);
for i=1:numancestor
    pp=postp_raw(i,:);
    pp=reshape(pp,numsnp0,numdescendant)';
    pp=reshape(repmat(pp,binSize,1),size(pp,1),binSize*size(pp,2));
    temp(:,:,i)=pp;
end
postp_raw=temp;


temp=zeros(numdescendant,numsnp0*binSize,numancestor);
for i=1:numancestor
    pp=postp_norm(i,:);
    pp=reshape(pp,numsnp0,numdescendant)';
    pp=reshape(repmat(pp,binSize,1),size(pp,1),binSize*size(pp,2));
    temp(:,:,i)=pp;
end
postp_norm=temp;
   
    
    
tempM.hidstate=hidstate;
tempM.par=parstore(end,:);
tempM.iterpar=parstore;
tempM.postp_raw=postp_raw; %���ǵ�ǰ��recombination
tempM.postp_norm=postp_norm;
tempM.observeProb=observeProb; % ֻ���ǹ۲����
M=tempM;



end
            