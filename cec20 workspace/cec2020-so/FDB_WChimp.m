% FDB-WChimp
%function [bestSolution, bestFitness, fes] = MRFO(fhd, dimension, maxFE, fNumber)
function [bestSolution, bestFitness,fes]=FDB_WChimp(fhd, dimension, maxFE, fNumber)

setting;
MaxIt=maxFE; 
nPop=30;
%F_index=fNumber;
Low=lbArray;
Up=ubArray;
Dim=dimension;
fes=1;
It=fes;
%ss=maxFE+1000;

%%%%SearchAgents_no=30; % Number of search agents%%% it is dimension
% initialize Attacker, Barrier, Chaser, and Driver
Attacker_pos=zeros(1,Dim); %//best solution
Attacker_score=inf; %change this to -inf for maximization problems // best fitness

Barrier_pos=zeros(1,Dim);
Barrier_score=inf; %change this to -inf for maximization problems

Chaser_pos=zeros(1,Dim);
Chaser_score=inf; %change this to -inf for maximization problems

Driver_pos=zeros(1,Dim);
Driver_score=inf; %change this to -inf for maximization problems


Convergence_curve=zeros(1,maxFE);

l=0;% Loop counter
%%
        for i=1:nPop   
        
            %Initialize the positions of search agents (population of chimp)
            PopPos(i,:)=rand(1,Dim).*(Up-Low)+Low;
                
            PopFit(i)= testFunction(PopPos(i,:)', fhd, fNumber);
        end
            BestF=inf;
            BestX=[];
        
        for i=1:nPop
            if PopFit(i)<=BestF
                BestF=PopFit(i);
                BestX=PopPos(i,:);
            end
        end
 
% Main loop
while (fes<maxFE)

     M = rand;
     c=rand;
     f=2-l*((2)/maxFE); % a decreases linearly from 2 to 0

     %Group 1
    C1G1=1.95-((2*l^(1/3))/(maxFE^(1/3)));
    C2G1=(2*l^(1/3))/(maxFE^(1/3))+0.5;
        
    %Group 2
    C1G2= 1.95-((2*l^(1/3))/(maxFE^(1/3)));
    C2G2=(2*(l^3)/(maxFE^3))+0.5;
    
    %Group 3
    C1G3=(-2*(l^3)/(maxFE^3))+2.5;
    C2G3=(2*l^(1/3))/(maxFE^(1/3))+0.5;
    
    %Group 4
    C1G4=(-2*(l^3)/(maxFE^3))+2.5;
    C2G4=(2*(l^3)/(maxFE^3))+0.5;
    r11=C1G1*rand(); % r1 is a random number in [0,1]
    r12=C2G1*rand(); % r2 is a random number in [0,1]
            
    r21=C1G2*rand(); % r1 is a random number in [0,1]
    r22=C2G2*rand(); % r2 is a random number in [0,1]
            
    r31=C1G3*rand(); % r1 is a random number in [0,1]
    r32=C2G3*rand(); % r2 is a random number in [0,1]
            
    r41=C1G4*rand(); % r1 is a random number in [0,1]
    r42=C2G4*rand(); % r2 is a random number in [0,1]
    A1=2*f*r11-f; % Equation (4)
    C1=2*r12; % Equation (5)
    m=chaos(3,1,1); % Equation (6) with Gauss
    
 % Update the Position of search agents including omegas
    for i=1:size(PopPos,1)
        for j=1:size(PopPos,2)     

            D_Attacker=abs(C1*Attacker_pos(j)-m*PopPos(i,j)); % Equation (8)
            X1=Attacker_pos(j)-A1*D_Attacker; % Equation (9)
                       
            A2=2*f*r21-f; % Equation (4)
            C2=2*r22; % Equation (5)
            
                   
            D_Barrier=abs(C2*Barrier_pos(j)-m*PopPos(i,j)); % Equation (8)
            X2=Barrier_pos(j)-A2*D_Barrier; % Equation (9)     
            
        
            
            A3=2*f*r31-f; % Equation (4)
            C3=2*r32; % Equation (5)
            
            D_Driver=abs(C3*Chaser_pos(j)-m*PopPos(i,j)); % Equation (8)
            X3=Chaser_pos(j)-A3*D_Driver; % Equation (9)      
            
            A4=2*f*r41-f; % Equation (4)
            C4=2*r42; % Equation (5)
            
            D_Driver=abs(C4*Driver_pos(j)-m*PopPos(i,j)); % Equation (8)
            X4=Chaser_pos(j)-A4*D_Driver; % Equation (9)
            
            W1=abs(X1)/(abs(X1)+abs(X2)+abs(X3)+abs(X4)); % Equation (10)
            W2=abs(X2)/(abs(X1)+abs(X2)+abs(X3)+abs(X4)); % Equation (11)
            W3=abs(X3)/(abs(X1)+abs(X2)+abs(X3)+abs(X4)); % Equation (12)
            W4=abs(X4)/(abs(X1)+abs(X2)+abs(X3)+abs(X4)); % Equation (13)
           

            PopPos(i,j)=1./(W1+W2+W3+W4).*(W1.*X1+W2.*X2+W3.*X3+W4.*X4)/4; % Equation (14)
            
        end
    end
    %FDB_index = fitnessDistanceBalance(PopPos, PopFit);

       if M<0.5                       
            newPopPos(i,:) = PopPos(i,:) + c.*(PopPos(i,:)-PopPos(1,:));
       else
           if 1/3<M && M<2/3
               newPopPos(i,:) = PopPos(i,:) + c.*(PopPos(i,:)- PopPos(randi([2 i-1]),:));
           else
               M1 = rand;
               newPopPos(i,:) = PopPos(i,:) + c.*(M1.*(PopPos(i,:)- PopPos(1,:))+(1-M1).*(PopPos(i,:)- PopPos(randi([2 i-1]),:)));
           end
       end
         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:nPop        
               newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
               %newPopFit(i)=BenFunctions(newPopPos(i,:),F_index,Dim);  

               newPopFit(i) = testFunction(newPopPos(i,:)',fhd,fNumber);
               fes = fes + 1;

              if newPopFit(i)<PopFit(i)
                 PopFit(i)=newPopFit(i);
                 PopPos(i,:)=newPopPos(i,:);
              end

              % Update Attacker, Barrier, Chaser, and Driver
            if  PopFit(i)<Attacker_score 
                Attacker_score= PopFit(i); % Update Attacker
                Attacker_pos=PopPos(i,:);
            end
            
            if  PopFit(i)> Attacker_score &&  PopFit(i)<Barrier_score 
                Barrier_score= PopFit(i); % Update Barrier
                Barrier_pos=PopPos(i,:);
            end
            
            if  PopFit(i)>Attacker_score &&  PopFit(i)>Barrier_score &&  PopFit(i)<Chaser_score 
                Chaser_score= PopFit(i); % Update Chaser
                Chaser_pos=PopPos(i,:);
            end
             if  PopFit(i)>Attacker_score &&  PopFit(i)>Barrier_score &&  PopFit(i)>Chaser_score &&  PopFit(i)>Driver_score 
                Driver_score= PopFit(i); % Update Driver
                Driver_pos=PopPos(i,:);
             end
    end


   
   
    
    l=l+1;    
    Convergence_curve(l)=Attacker_score;
    for i=1:nPop
        if PopFit(i)<BestF
           BestF=PopFit(i);
           BestX=PopPos(i,:);            
        end
    end

    bestSolution=Attacker_pos;
    bestFitness=Attacker_score;
end

end


