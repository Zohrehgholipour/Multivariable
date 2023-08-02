clc
clear
close all
%% Data
Data=readtable('data.xlsx');
p=50;
n=100;
V=[ Data.z1	    Data.z2	    Data.z3 	Data.z4	    Data.z5	    Data.z6	 Data.z7	    Data.z8	    Data.z9     Data.z10	Data.z11	Data.z12  Data.z13	Data.z14	Data.z15	Data.z16	Data.z17	Data.z18  Data.z19	Data.z20	Data.z21	Data.z22	Data.z23	Data.z24  Data.z25	Data.z26	Data.z27	Data.z28	Data.z29	Data.z30  Data.z31	Data.z32	Data.z33	Data.z34	Data.z35	Data.z36  Data.z37	Data.z38	Data.z39	Data.z40	Data.z41	Data.z42  Data.z43	Data.z44	Data.z45    Data.z46	Data.z47	Data.z48  Data.z49	Data.z50];
%% (1),(2)
Vbar=mean(V);                               % Population mean vector
Covp=(n-1)/n*cov(V);                        % Population Covariance Matrix 
%% (3)
                                            % Population Correlation matrix  
for i=1:p
    for j=1:p
        Rp(i,j)=Covp(i,j)/sqrt(Covp(i,i)*Covp(j,j));
    end
end
%% (4),(5)
Covs=cov(V);                                % Sample Covariance Matrix 
Rs=corr(V);                                 % Sample correlation Matrix

[EVe,EVa]=eig(Rs); 
                                 
landa=diag(EVa);                            % Eigenvalue    
[Landa,sortorder]=sort(landa,'descend');    % SORT Eigenvalue
e=EVe(:,sortorder);                         % Update Eigenvector
%% (6)
plot(sortorder,Landa,'r')                   % Scree plot

xlabel('(i)')
ylabel('Landa(i)')
title('scree plot')
%% (7)
S=sum(landa);                               % Sum of Eigenvalue
CS=cumsum(landa);                           % Cumulative sum of Eigenvalue
J=CS./S;                                    % Discription of dispersion
ind=find(J>0.95);                           % find of 0.95 dispersion

%% (8)
z=zscore(V);                                % Normal Standard
Y=z*e;                                      % Principle scores
%% (9),(10)
e(:,16:50)=[];                              % Eliminate redundant principle components
Y(:,16:50)=[];
Z=Y*e';                                     % Debugged observation matrix
%% (11)
Y(:,2:15)=[];
qqplot(Y);                                  % Q-Q plot

%% (12)
[EVe,EVa]=eig(Rs);                   
landa=diag(EVa);                             
[Landa,sortorder]=sort(landa,'descend');    
e=EVe(:,sortorder);                         
Y=z*e;
Y(:,16:50)=[];
Landa(16:50)=[];


tsquared=(((Y(:,1)).^2)/Landa(1))+(((Y(:,2)).^2)/Landa(2))+(((Y(:,3)).^2)/Landa(3))+(((Y(:,4)).^2)/Landa(4))+(((Y(:,5)).^2)/Landa(5))+(((Y(:,6)).^2)/Landa(6))+(((Y(:,7)).^2)/Landa(7))+(((Y(:,8)).^2)/Landa(8))+(((Y(:,9)).^2)/Landa(9))+(((Y(:,10)).^2)/Landa(10))+(((Y(:,11)).^2)/Landa(11))+(((Y(:,12)).^2)/Landa(12))+(((Y(:,13)).^2)/Landa(13))+(((Y(:,14)).^2)/Landa(14))+(((Y(:,15)).^2)/Landa(15));



L=chi2inv(0.95,15);                           % upper limit
UCL=L*ones(1,100);                   

plot(UCL,'r.');                               % UCL plot

hold on

plot(tsquared)                                % T^2 Hotelling plot

xlabel('i')
ylabel('T^2')
title('T^2 Hotelling plot')





