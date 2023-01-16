function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)
    A = [A;-A];
    b = [bu;-bl];
    options = optimset('MaxIter',1e3,'Display','on');
    [x,info] = quadprog(H,g,A,b,[],[],l,u,xinit,options);
end

