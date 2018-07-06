function [DoFMax,H,msg] = montecarlo(n)
p=linspace(0,1,10000);  %probability vector
DoF=zeros(100,length(p));
for i=1:99
    tt=gcd(i,100);%greatest common devider of I and K=100 users
    f=i/100;      % fraction f(p)
    K=100+1; %K=K-1 %K=100 user + transmitter 0 since we considered receiver 0 as inactive. Thus the index "i" of the channel topology is mapped to "i+1" in the simulation
    N=(K-1)/tt+1; %New network size of for the assignment strategy

%% Message Assignments
        msg = zeros(K,2); 

        msg = assignment_new(f,N-1);% messages are assigned according to the assignment strategy in the table

        for j=1:length(p)    % different values for p
            DoF_k = zeros(1, n); % DoF for a particular probability p, with n realizations of H-matrix
 
            for k=1:n     % n different channel realizations
            % Generate random H-matrix where links exist with probability (1-p).
            % Uniformly distributed random variable on [0,1] is greater than p(j) with
            % probability (1-p).
            H = zeros(K,K); %channel matrix, H(i,j) denotes the link from Txj to Rxi  
  
                 for l=1:K %simulating the channel realization for one iteration step                    
                    if not(l==K)
                    H(l+1,l)= rand >= p(j);
                    end
                    H(l,l)= rand >= p(j);
                 end
                 for l=1:K %identifying all inactive transmitter and erasing the links that are connected to those transmitters to decrease the duration of the simulation 
                      if (not(any(msg(:,1)==l))&&not(any(msg(:,2)==l))) 
                         if not(l==K)
                            H(l+1,l)=0;
                         end
                         H(l,l)=0;
                     end
                 end
                H(1,1)=0;% Erasing the first link since we have 100 user + transmitter 0. Receiver 0 is inactive.        

                for l1=2:K  % starting to scan the network with H_21 and identifying all subnetworks
                    for l2=(l1-1):l1
                             if((l1==l2) && (H(l1,l2)==1) && (H(l1,l1-1)==0)) % A direct link exists and a diagonal link that shares the same receiver is erased. Hence the direct link is the first link of the identified subnetwork
                                    counter=l1; %marking the first user of the identified subnetwork
                                stop=0;
                                    for l11=(l1+1):K %finding the first existing link, where the next consequtive link is erased. Since the first existing link is direct, we scan links pairwise i.e. H(l11,l11-1) and H(l11,l11)
                                        for l22=(l11-1):l11 % considering a pair of links which consists of H(l11,l11-1) and H(l11,l11)
                                            if ((l11>l22)&&(H(l22,l22)==1)&&(H(l11,l22)==0))  %diagonal link is erased
                                                if l22==counter % special case if subnetwork is of size 1
                                                    temp=1;
                                                else
                                                    temp=l22-counter+1; %get size of scanned cluster that is greater then 1
                                                end
                                                    DoF_k(k)=DoF_k(k)+compute_DoF_if(counter,temp, msg, 0,1); %Calling the function that calculates the DoF for the identified subnetwork
                                                    stop=1;
                                            end
                                            if ((l11==l22)&&(H(l22,l22-1)==1)&&(H(l11,l22)==0)&&not(l11==K)&&not(l22==K)) %direct link is erased
                                                if l11==counter % speciale case for a cluster of size 1
                                                     temp=1;
                                                else
                                                    temp=l11-counter+1;    %get size of scanned cluster
                                                end
                                                DoF_k(k)=DoF_k(k)+compute_DoF_if(counter,temp, msg, 0,0);
                                                stop=1;
                                            end
                                            if ((l11==K)&&(l22==K)&&(H(l1,l2)==1)) %Special case if all links of the whole network exist
                                                if l11==counter
                                                     temp=1;
                                                else
                                                    temp=l11-counter+1;    
                                                end
                                                stop=1;
                                                DoF_k(k)=DoF_k(k)+compute_DoF_if(counter,temp, msg, 0,1);
                                            end
                                        if stop==1
                                        break;
                                        end
                                        end
                                        if stop==1
                                        break;
                                        end
                                    end
                            end
                             if((l1>l2) && (H(l1,l2)==1) && (H(l2,l2)==0))% Identifying that the first existing link of the subnetwork is a diagonal link
                                counter=l1; 
                                stop=0;
                                    for l22=l2+1:K %Since the first existing link in the subnetwork is diagonal, we start with the folloiwing direct link to scan the consequtive links of the subnetwork
                                        if l22==K  %If the first existing link is H(K,K-1)=1, then the we only have to scan the last direct link H(K,K). If not, then we scan the network pairwise by checking H(l22,l22) and H(l22+1,l22).
                                            a=l22;
                                        else
                                            a=l22+1;
                                        end
                                        for l11=l22:a
                                            if ((l11>l22)&&(H(l22,l22)==1)&&H(l11,l22)==0)  
                                                if l22==counter
                                                    temp=1;
                                                else
                                                    temp=l22-counter+1;
                                                end
                                                DoF_k(k)=DoF_k(k)+compute_DoF_if(counter,temp, msg,1, 1);
                                                stop=1;
                                            end
                                            if ((l11==l22)&&(H(l22,l22-1)==1)&&(H(l11,l22)==0)&&not(l11==K)&&not(l22==K))
                                                if l22==counter
                                                    temp=1;
                                                else
                                                    temp=l11-counter+1;
                                                end
                                                DoF_k(k)=DoF_k(k)+compute_DoF_if(counter,temp, msg,1, 0);
                                                stop=1;
                                            end
                                            if ((l11==K)&&(l22==K)&&(H(l1,l2)==1))
                                                if l22==counter
                                                    temp=1;
                                                else
                                                    temp=l11-counter+1;
                                                end
                                                stop=1;
                                                DoF_k(k)=DoF_k(k)+compute_DoF_if(counter,temp, msg, 1,1);
                                            end
                                     if stop==1
                                        break;
                                     end
                                        end
                                     if stop==1
                                        break;
                                     end
                                    end
                                         
                            end
                    end
                end
                DoF_k(k)=DoF_k(k)/(K-1); %calculating the puDoF of the whole network
            end
%Compute average DoF
              DoF(i,j)=sum(DoF_k)/n; %averaging the puDoF over all iteration points n for one fixed message assignment
        end
end

DoFMax=max(DoF); %storing the maximum puDoF of all possible message assignments for every value of p
end
