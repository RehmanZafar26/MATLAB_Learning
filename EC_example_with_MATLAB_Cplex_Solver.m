%This example is coded by REHMAN ZAFAR for under standing about MATLAB solver and CEPLEX solvers
% The example 3C  of Econimic Dispatch as LP from book "Power Generation, Operation and Control" by Allenj. Wood et al. is implemented.

clear all

n = input('For MatlabinputZ solver=1/For CPLEX =2: ');



cost= [561 7.92 0.001562; 310 7.85 0.00194; 78 7.97 0.00482];
Glim = [150 600;100 400;50 200];
numSegment =input('input the number of segment = ');
Pload = 850 % input('input the load demand = ');
num_of_gen=3;

Aeq =ones(1,numSegment*size(cost,1));
beq = Pload - sum(Glim(:,1));

Pedge = zeros(size(cost,1),numSegment+1);
slop = zeros(size(cost,1),numSegment);
delPgen=zeros(size(cost,1),numSegment);
for i=1:size(cost,1)
    Pedge(i,:)= linspace(Glim(i,1),Glim(i,2),numSegment+1);
    for j=1:numSegment
        limCost= cost(i,:)*[1 1; Pedge(i,j+1) Pedge(i,j); Pedge(i,j+1)^2 Pedge(i,j)^2];
        delCost = limCost(1)-limCost(2);
        delPedge = Pedge(i,j+1)-Pedge(i,j);
        delPgen(i,j) = delPedge;
        slop(i,j)= delCost/delPedge;
    end
end
fmin=[];
b=[];
for i=1:size(slop,1)
    fmin =[fmin slop(i,:)] ;
    b= [b delPgen(i,:)];
end
A = zeros(size(cost,1),numSegment*size(cost,1));
j=1;
for i=1:numSegment*size(cost,1)
    A(i,j)=1;
    j=j+1;
end
A = [A;-A];

b = [b zeros(1,length(b))];

switch n
    case 1
        tic
        [P,fval,exitflag,output] = linprog(fmin,A,b',Aeq,beq);
        toc
        disp ('MATLAB solver solution')
    case 2
      tic
        [P,fval,exitflag,output] = cplexlp(fmin,A,b',Aeq,beq);
        toc
    disp ('CPLEX solution')
end

%


j=1;
for i=1:size(cost,1)
    Pgen(i)= sum(P(j:j+numSegment-1,1))+Glim(i,1);
    j=j+numSegment;
    
end
Total_Power=Pgen;
nHours=1;
num_of_gen=3;
  sub_Hourly_Cost = zeros(nHours, num_of_gen);
for h = 1:1:nHours
    
    for num = 1:1:num_of_gen
        
        sub_Hourly_Cost(h, num) = ((cost(num, 3) * Total_Power(h, num) * Total_Power(h, num)) + (cost(num, 2) * Total_Power(h, num)) + cost(num, 1));
        
    end
    
end
sub_Hourly_Cost = sub_Hourly_Cost;
%Find Total Cost of Operation
Hourly_Cost = zeros(nHours, 1);
for h = 1:1:nHours
    for num = 1:1:num_of_gen
        Hourly_Cost(h, 1) = Hourly_Cost(h, 1) + sub_Hourly_Cost(h, num);
    end
end
totalcost=sum(Hourly_Cost)
P1=Pgen(:,1)
P2=Pgen(:,2)
P3=Pgen(:,3)
