function [msgshift] = assignment_new(f,K)
%% Assigning messages according to our strategy
% We pass the variables "f" and "K", where K is the size of the subnetwork.
% Furthermore we seperate K1=100 user into groups of K user, such that the
% number of iterations in montecarlo.m automatically increases by the factor of number of these
% existing groups. Each of these groups will have the same massage
% assignment modulo 5.
K1=100;
msg = zeros(K1,2); % msg denotes the matrix of the message assignment, where transmitter 0 is not in the network.
msgshift = zeros(K1,2);% in msgshift we consider transmitter 0.
 for i=1:K1 % First we consider the message assignments where transmitter 0 is not in the network.
     if mod(i,K)==1
        msg(i,:)=[i,i+1];
     end
     if mod(i,K)==0
     msg(i,:)=[i-2,i-1];
     end 
     limit=(min(f*K-2, floor(K/2-1)));
     if limit>0
        n=linspace(1,limit,limit);
        t=(max(2,floor(K/(f*K-1))));
        for j=1:length(n)
            s=floor(1+n(1,j)*t);
            if mod(i,K)==s
            msg(i,:)=[i,i+1];
            end
        end
     end
    limit1=ceil((f-1/2)*K)-1;
    if limit1>0
        n1=linspace(1,limit1,limit1);
        for j=1:length(n1)
            s=floor(2*n1(1,j));
            if mod(i,K)==s
            msg(i,:)=[i,i+1];
            end
        end
    end
 for i=1:K1
     if ((msg(i,1)==0) && (msg(i,2)==0))
         msg(i,:)=[i-1,i];
     end
 end
 if f==0.01
 msg(1,:)=[0,1];
 end
 for s=2:K1+1 %we shift the message assignments by 1 such that transmitter 0 is in the network
     msgshift(s,:)=msg(s-1,:)+1;
 end
end
