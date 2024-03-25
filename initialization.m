

% This function initialize the first population of search agents
function Positions=initialization(pop,dim,d,c)

Boundary_no= size(d,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(pop,dim).*(d-c)+c;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=d(i);
        lb_i=c(i);
        Positions(:,i)=rand(pop,1).*(ub_i-lb_i)+lb_i;
    end
end