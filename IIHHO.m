
%%
function [Rabbit_Energy,Rabbit_Location,CNVG]=IIHHO(N,T,lb,ub,dim,fobj)
% % disp('HHO is now tackling your problem')
% % tic
% initialize the location and Energy of the prey(rabbit)
Rabbit_Location=zeros(1,dim);
Rabbit_Energy=inf;
%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);
CNVG=zeros(1,T);
t=0; % Loop counter
while t<T
    for i=1:size(X,1)
        % Check boundries边界限制
        FU=X(i,:)>ub;
        FL=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness=fobj(X(i,:));
        % Update the location of Rabbit
        if fitness<Rabbit_Energy
            Rabbit_Energy=fitness;
            Rabbit_Location=X(i,:);
        end
    end
     
  
    E1=2*exp(-pi*(t/T));
    for i=1:size(X,1)
        CC=cos((10*pi*t)/T);
        E0=2*rand()-1; %-1<E0<1
        if CC>=0
              Escaping_Energy=E0*E1*CC;  % escaping energy of rabbit
        elseif CC<0
              Escaping_Energy=0;
        end

%% 开始判断动态逃逸能量  abs用作绝对值
        if abs(Escaping_Energy)>=1
           
%% Exploration（探索）:
% Harris' hawks perch randomly based on 2 strategy:
            q=rand();
            rand_Hawk_index = floor(N*rand()+1);
            X_rand=X(rand_Hawk_index,:);
            if q<0.5
                % perch based on other family members
                X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
            elseif q>=0.5
                % perch on a random tall tree (random site inside group's home range)
                % 栖息在一颗随机的高树上（群体范围内的随机点）
                X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand()+lb);
            end
            
        elseif abs(Escaping_Energy)<1
            
%% Exploitation（开发）:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
%% phase 1: surprise pounce (seven kills)
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event
            RR=Rabbit_Location;
%           SS=rand()*(sin((pi/2)*(t/T))+cos((pi/2)*(t/T))-1);
            if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                X(i,:)=((RR)-Escaping_Energy*abs(RR-X(i,:)));
            end
            
            if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege（软围攻）
                Jump_strength=2*(1-rand());
                X(i,:)=((RR-X(i,:))-Escaping_Energy*abs(Jump_strength*RR-X(i,:)));
            end
            
%% phase 2: performing team rapid dives (leapfrog movements)
            if r<0.5 && abs(Escaping_Energy)>=0.5, % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                Jump_strength=2*(1-rand());
                X1=(RR-Escaping_Energy*abs(Jump_strength*RR-X(i,:)));
                
                if fobj(X1)<fobj(X(i,:)) % improved move?
                    X(i,:)=X1;
                else % hawks perform levy-based short rapid dives around the rabbit
                    X2=RR-Escaping_Energy*abs(Jump_strength*RR-X(i,:))+Levy(dim,t,T)*((1-(t/T)).^2);
                    % rand(1,dim)即为S，是大小为1xD的随机向量，Levy飞行函数
                    if (fobj(X2)<fobj(X(i,:))), % improved move?
                        X(i,:)=X2;
                    end
                end
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5, % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                % 渐进快速俯冲的硬围攻，用于兔子没有足够的能量逃跑的情况
                % hawks try to decrease their average location with the rabbit
                % 老鹰尽量缩小自己与兔子的“平均位置”                      
                Jump_strength=2*(1-rand());
                X1=(RR-Escaping_Energy*abs(Jump_strength*RR-mean(X)));      
                if fobj(X1)<fobj(X(i,:)) % improved move?
                    X(i,:)=X1;
                else % Perform levy-based short rapid dives around the rabbit    rand(1,dim).*      *((1-(t/T)).^2)
                    X2=RR-Escaping_Energy*abs(Jump_strength*RR-mean(X))+Levy(dim,t,T)*((1-(t/T)).^2);
                    if (fobj(X2)<fobj(X(i,:))), % improved move?
                        X(i,:)=X2;
                    end
                end
            end
%%
        end
    end
    t=t+1;
    CNVG(t)=Rabbit_Energy;

end
% toc
end

% ___________________________________
function o=Levy(d,t,T)
beta=1.5*(1-(-((t/T)/2)+((t/(2*T)).^2+(t/(3*T)).^3).^(1/2)).^(1/3));
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
if beta>1.2
      o=step;
elseif beta<1.2
      o=step*(rand());
end
end