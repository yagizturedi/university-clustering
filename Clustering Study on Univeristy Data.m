
%% Yaðýz Türedi 2095719

clear
clc
warning('off')
%% Reading the data

[DATANUM,~,DATACELL] = xlsread("Universities",1,"A4:W1305");

Private=DATANUM(:,1)-1;

DATANUM(:,1)=[];

States=unique(DATACELL(:,2));

%% Normalizing the data

for i=1:20
    maximum=max(DATANUM(:,i));
    minimum=min(DATANUM(:,i));
    for j=1:length(DATANUM)
        DATANUM(j,i)=(DATANUM(j,i)-minimum)/(maximum-minimum);
    end
end



%% Cleaning the data
n=1;
m=1;
for i=1:length(DATANUM)
    t=0;
    for j=1:20
        if ismissing(DATANUM(i,j))
            
            t=1;
        end
    end
    
    if t==1
        Missing(n,:)=DATANUM(i,:);
        if i==476
            TuftsIndex=n;
        end
        n=n+1;
    else
        Complete(m,:)=DATANUM(i,:);
        CompleteCell(m,:)=DATACELL(i,:);
        CompletePrivate(m)=Private(i);
        
        m=m+1;
    end
    
end

%% Applying k-means to the clean data

for k=1:70
    minError=Inf;
    for ttt=1:4
        [IDX,C,SUMD]=kmeans(Complete,k,'Start','uniform','Replicates',10);

        %silhouette calculation for best k
        
        S=silhouette(Complete,IDX);
        error=sum(S)/length(Complete);


%     sse calculation
%         error=0;
%         for t=1:length(DATANUM)
%             if sum(ismissing(DATANUM(t,:)))
%                 continue;
%             end
%             
%             error=error+(EuclideanDist(DATANUM(t,:),C(IDX(t),:)))^2;
%         end
%         if error<minError
%             minError=error;
%         end
%     end

    kErrors(1,k)=k;
    kErrors(2,k)=error;
    
    if max(kErrors(2,:))==error
    
        IDXbest=IDX;
        Cbest=C;
    end
    end
end
maxerror=max(kErrors(2,:));
bestK=find(kErrors(2,:)==maxerror);
plot(kErrors(1,:),kErrors(2,:));

ttt=ones(length(Complete),1);
Complete=[[ttt],Complete];

for i=1:length(Complete)
    Complete(i,1)=i;
end

n=ones(max(IDXbest),1);
for i=1:max(IDXbest)
    for datapoint=1:length(IDXbest)
        if IDXbest(datapoint)==i
            cluster(i,n(i))=datapoint;
            n(i)=n(i)+1;
        end

    end
end
%% Getting statistics for a cluster

clusternumber=[3];
for i=1:length(clusternumber)
    clear ClusterPoints
    points=cluster(clusternumber(i),:);
    for t=1:length(points)
        if points(t)~=0
        ClusterPoints(t,:)=Complete(points(t),:);
        end
    end
    
    for t=1:length(ClusterPoints(1,:))
        InCluster(1,t)=max(ClusterPoints(:,t));
        InCluster(2,t)=mean(ClusterPoints(:,t));
        InCluster(3,t)=min(ClusterPoints(:,t));
    end

end

%% Comparison with categorical data

    % Percentage of private schools in cluster i
    
for i=1:max(IDXbest)

    
        points=cluster(i,:);
        clear ClusterPoints
    for t=1:length(points)
        if points(t)~=0
        ClusterPoints(t,:)=CompleteCell(points(t),:);
        end
    end
    for e=1:length(ClusterPoints(:,1))
        
        ClusterPrivate(e)=ClusterPoints{e,3}-1;
        
    end
    PrivatePerc(i)=mean(ClusterPrivate);
end

scatter(1:max(IDXbest),PrivatePerc)


    % Clusters of schools in each state
    
for i=1:length(Complete)
    clear idx
    idx = all(ismember(States,CompleteCell{i,2}),2);
    StateNumb(i)=find(idx);
    
end

scatter(StateNumb,IDXbest);

Tufts=Missing(408,:);

for i=1:length(Cbest(:,1))
    n=1;
    for j=1:length(Tufts)
        if ismissing(Tufts(j))
            continue;
        else
            diffvector(i,n)=Tufts(j)-Cbest(i,j);
            n=n+1;
        end
    end
    
    distance(i)=(diffvector(i,:)*diffvector(i,:)')^0.5;
end
