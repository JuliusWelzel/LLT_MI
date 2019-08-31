function SpacFilt = L1_SVD_CSP(EEGSignals,Label,FilterNum)
%
% <The L1_SVD_CSP 
%this function learn a L1_SVD_CSP (L1 Norm based Singular Value Decomposition Common Spatial Patterns) filters, using L1_Norm based SVD,
%to discriminate two mental states in EEG signals.
%
%Input:
%EEGSignals: the training EEG signals, composed of 2 classes. These signals
%are a structure such that:
%   EEGSignals: the EEG signals as a [Ns * Nc * Nt] Matrix where
%       Ns: number of EEG samples per trial
%       Nc: number of channels (EEG electrodes)
%       Nt: number of trials
%Label: a [1 * Nt] vector containing the class labels for each trial
%FilterNum: The number of CSP filter numbers
%Output:
%SpacFilt: the learnt CSP filters (a [(2*Nc) * Nc] matrix with the spatial filters as rows)
%
%created: 06/07/2013
%last revised: 12/23/2013
%


%check and initializations
%
ClassNum = length(unique(Label));
ChanNum=size(EEGSignals,2);
LengTrailWin=size(EEGSignals,1);
TrialNum=size(EEGSignals,3);
lab = Label;
StNum = 0;

classData={};
for c=1:ClassNum
    switch c
        case 1
            flag=-1;
        case 2
            flag=1;
        otherwise
            flag=0;
    end
    class{c}=find(lab==flag);
    classSize(c)=length(class{c});
    classData{c}(:,:,1:classSize(c))= EEGSignals(1:LengTrailWin,:,class{c});
end

% Ramoser equation (1)
for c=1:ClassNum
    Cmean{c}=zeros(ChanNum);
    for trailN=1:classSize(c)
        timeSeq= squeeze(classData{c}(:,:,trailN));
        timeSeq= timeSeq-repmat(mean(timeSeq),LengTrailWin,1);% centered
        tempC=timeSeq'*timeSeq;%covariance matrices
        C=tempC./trace(tempC);
        Cmean{c}=Cmean{c}+C;
    end
    Cmean{c}=Cmean{c}./classSize(c);
end
VsPair = 3;
fprintf('Classes VS pairs("0" means the others): \n');
for i=0:ClassNum
    for j=i+1:ClassNum
        vs{c}=[i j];
        fprintf('%2d: %d vs %d.\n',c,[i j]);
        c=c+1;
    end
end
fprintf('Choose one VS pair.\n');

num_vs=size(vs,1);


% Ramoser equation (2)
Ccomposite=zeros(ChanNum);
switch vs{VsPair}(1)
    case 0
        for c=1:ClassNum
            Ccomposite=Ccomposite+Cmean{c};
        end
    otherwise
        Ccomposite=Cmean{vs{VsPair}(1)}+Cmean{vs{VsPair}(2)};
end

% Sort eigenvalues in descending order
[Ucomposite,Lambdacomposite] = eig(Ccomposite);
[Lambdacomposite,ind] = sort(diag(Lambdacomposite),'descend');
Ucomposite = Ucomposite(:,ind);

% Ramoser equation (3) - Whitening transform
P=sqrt(inv(diag(Lambdacomposite)))*Ucomposite';

% Ramoser equation (4)
if vs{VsPair}(1)==0
    S{1}=P*Cmean{vs{VsPair}(2)}*P';
    S{2}=P*(Ccomposite-Cmean{vs{VsPair}(2)})*P';
else
    S{1}=P*Cmean{vs{VsPair}(1)}*P';
    S{2}=P*Cmean{vs{VsPair}(2)}*P';
end

Big_Matrix1=inv( S{2} + 0.001* eye(size(S{2})) )*S{1};
Big_Matrix2=inv( S{1} + 0.001* eye(size(S{1})) )*S{2};

[B,D] = eig(Big_Matrix1);
[D,ind] = sort(diag(D),'descend');
B = B(:,ind);
B1_ini = B;
B1=L1_SVD(Big_Matrix1,B1_ini);
B_max=max(B1);
for i=1:1:length(B_max)
    B1(:,i)=B1(:,i)/B_max(i);
end;  
W1=(B1'*P);
for i=1:length(ind)
    W1(i,:)=W1(i,:)./norm(W1(i,:));
end


B = []; D = [];
[B,D] = eig( Big_Matrix2 );
[D,ind] = sort(diag(D),'descend');
B = B(:,ind);
B2_ini = B;
B2=L1_SVD(Big_Matrix2,B2_ini);
B_max=max(B2);
for i=1:1:length(B_max)
    B2(:,i)=B2(:,i)/B_max(i);
end;
W2=(B2'*P);
for i=1:length(ind)
   W2(i,:)=W2(i,:)./norm(W2(i,:));
end


SpacFilt=[];

for i=1:FilterNum
    SpacFilt=[SpacFilt; W1(i,:);W2(i,:)];
end;

SpacFilt=SpacFilt';

end


function output=L1_SVD(Feature_Matrix,w_ini)
% <The L1_SVD
%this function learn a L1 Based SVD feature eigenvectors
%Input:
%Feature_Matrix: the training EEG signals, composed of 2 classes. These signals
%are a structure such that:
%   EEGSignals: the EEG signals as a [Nc * Nc] Matrix where
%       Nc: number of channels (EEG electrodes)
%       Nc: number of channels (EEG electrodes)
%w_ini: The L2 Norm Based SVD eigenvectors
%
%Output:
%output: L1 Norm Based SVD eigenvectors
%
%created: 06/07/2013
%last revised: 12/23/2013
ini_x=Feature_Matrix;         
C=ini_x;               
[m,n]=size(C);
sum1=0;
for i=1:1:n
    sum1=norm(C(:,i))^2+sum1;
end;
s=sum1/n;             
error=[];
%%%%%%%%%%%%%%%%%%取L2模最大的一列向量为初始值%%%%%%%%%%%%%%%%%%%%%
w=w_ini(:,1);
w_store=[];
judge_2=[];
judge_2(1)=0;
feature_number=size(w_ini,2);
s_1=zeros(1,feature_number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:1:feature_number
    w=w(:,j);
    w(:,j+1)=L1_Norm_Maximized(w,C);
    w(:,j+1)=w(:,j+1)/norm(w(:,j+1));
    w_store(:,j)=w(:,j+1);

    sum2=0;
    for i=1:1:n
        sum2=sum2+(w_store(:,j)'*ini_x(:,i))^2;
    end;
    feature_x(j)=sum2/n;
    for num=1:1:j
       s_1(j)=s_1(j)+feature_x(num);
    end;
    judge_2(j+1)=s_1(j)/s; 
    feature(j,:)=w(:,j+1)'*ini_x;

    for i=1:1:n
        C(:,i)=C(:,i)-w(:,j+1)*(w(:,j+1)'*C(:,i));
    end;
end; 
output=w_store(:,1:end);
end

function output=L1_Norm_Maximized(w,C)

[m,n]=size(C);
middle_w(:,1)=w;
for i=1:1:20
    flag1=0;
    x_sum=zeros(m,1);
    for k=1:1:n
        if middle_w(:,i)'*C(:,k)<0               %Polarity check
            P(k)=-1;
        else
            P(k)=1;
        end;
    end;
    for j=1:1:length(P)
        x_sum=x_sum+P(j)*C(:,j);
    end;
    middle_w(:,i+1)=x_sum/norm(x_sum);
    for j=1:1:n
        if abs(middle_w(:,i+1)'*C(:,j))<=1.0e-8
            flag1=1;
            break;
        end;
    end;
    if middle_w(:,i+1)~=middle_w(:,i)
        continue;
    elseif flag1==1
        vvv=ones(m,1)*0.000001;
        middle_w(:,i+1)=(middle_w(:,i+1)+vvv)/(norm(middle_w(:,i+1)+vvv));
        continue;
    else 
        break;
    end;
end;
output=middle_w(:,i+1);
end
