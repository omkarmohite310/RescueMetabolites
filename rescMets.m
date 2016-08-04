function [ResMet] = rescMets(model,Lethal, order)
%To find the rescue metabolites for corresponding lethal reaction set
%Input
% model (the following fields are required - others can be supplied)       
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
%   mets         Metabolite Names

% Lethal         Structure of lethal reaction sets
% Jsl            Single Lethals
% Jdl            Double Lethals
% Jtl            Triple Lethals

% order          Order of lethal sets

% initCobraToolbox;
solWT=optimizeCbModel(model,'max','one');
grWT=solWT.f;
modeldel=model;
% Single lethal reactions and rescue metabolites
if order==1
    ResMet.SL=[];
    for iLeth=1:length(Lethal.Jsl)
        del_Idx= find(strcmp(Lethal.Jsl(iLeth),model.rxns));
        modeldel.lb(del_Idx)=0;
        modeldel.ub(del_Idx)=0;
        solKO_i=optimizeCbModel(modeldel);
        if solKO_i.f<0.01*grWT
           metAll=find(modeldel.S(:,del_Idx));
           met=[];
           if modeldel.rev(del_Idx)==0
               for iMet=1:length(metAll)
                   if modeldel.S(metAll(iMet),del_Idx)>0
                       met=[met;modeldel.mets(metAll(iMet))];
                   end
               end
           else    
               met=modeldel.mets(find(modeldel.S(:,del_Idx)));
           end    
           for jMet=1:length(met)
               modelExc=addExchangeRxn(modeldel,met(jMet));
                  solRM_ij=optimizeCbModel(modelExc,'max','one');
                   if(solRM_ij.f>0.01*grWT)
                       ResMet.SL=[ResMet.SL;model.rxns(del_Idx),met(jMet),solRM_ij.f];
                   end
           end
        end
        modeldel.lb(del_Idx)=model.lb(del_Idx);
        modeldel.ub(del_Idx)=model.ub(del_Idx);
    end
end

% Double and single lethal reactions and rescue metabolites
if order==2
    ResMet.SL=[];
    ResMet.DL=[];
    for iLeth=1:length(Lethal.Jsl)
        del_Idx= find(strcmp(Lethal.Jsl(iLeth),model.rxns));
        modeldel.lb(del_Idx)=0;
        modeldel.ub(del_Idx)=0;
        solKO_i=optimizeCbModel(modeldel);
        if solKO_i.f<0.01*grWT
           metAll=find(modeldel.S(:,del_Idx));
           met=[];
           if modeldel.rev(del_Idx)==0
               for iMet=1:length(metAll)
                   if modeldel.S(metAll(iMet),del_Idx)>0
                       met=[met;modeldel.mets(metAll(iMet))];
                   end
               end
           else    
               met=modeldel.mets(find(modeldel.S(:,del_Idx)));
           end 
           for jMet=1:length(met)
               modelExc=addExchangeRxn(modeldel,met(jMet));
                  solRM_ij=optimizeCbModel(modelExc,'max','one');
                   if(solRM_ij.f>0.01*grWT)
                       ResMet.SL=[ResMet.SL;model.rxns(del_Idx),met(jMet),solRM_ij.f];
                   end
           end
        end
        modeldel.lb(del_Idx)=model.lb(del_Idx);
        modeldel.ub(del_Idx)=model.ub(del_Idx);
    end
    
    for iLeth=1:length(Lethal.Jdl)
        del_Idx1= find(strcmp(Lethal.Jdl(iLeth,1),model.rxns));
        del_Idx2= find(strcmp(Lethal.Jdl(iLeth,2),model.rxns));
        modeldel.lb(del_Idx1)=0;
        modeldel.ub(del_Idx1)=0;
        modeldel.lb(del_Idx2)=0;
        modeldel.ub(del_Idx2)=0;
        solKO_i=optimizeCbModel(modeldel);
        if solKO_i.f<0.01*grWT
           metAll1=find(modeldel.S(:,del_Idx1));
           met1=[];
           if modeldel.rev(del_Idx1)==0
               for iMet=1:length(metAll1)
                   if modeldel.S(metAll1(iMet),del_Idx1)>0
                       met1=[met1;modeldel.mets(metAll1(iMet))];
                   end
               end
           else    
               met1=modeldel.mets(find(modeldel.S(:,del_Idx1)));
           end 
           metAll2=find(modeldel.S(:,del_Idx2));
           met2=[];
           if modeldel.rev(del_Idx2)==0
               for iMet=1:length(metAll2)
                   if modeldel.S(metAll2(iMet),del_Idx2)>0
                       met2=[met2;modeldel.mets(metAll2(iMet))];
                   end
               end
           else    
               met2=modeldel.mets(find(modeldel.S(:,del_Idx2)));
           end 
           met=unique([met1;met2]);
           
           for jMet=1:length(met)
               modelExc=addExchangeRxn(modeldel,met(jMet));
                  solRM_ij=optimizeCbModel(modelExc,'max','one');
                   if(solRM_ij.f>0.01*grWT)
                       ResMet.DL=[ResMet.DL;model.rxns(del_Idx1),model.rxns(del_Idx2),met(jMet),solRM_ij.f];
                   end
           end
        end
        modeldel.lb(del_Idx1)=model.lb(del_Idx1);
        modeldel.ub(del_Idx1)=model.ub(del_Idx1);
        modeldel.lb(del_Idx2)=model.lb(del_Idx2);
        modeldel.ub(del_Idx2)=model.ub(del_Idx2);
    end
end    

% Triple lethal reactions and rescue metabolites
if order==3
    ResMet.SL=[];
    ResMet.DL=[];
    ResMet.TL=[];
    for iLeth=1:length(Lethal.Jsl)
        del_Idx= find(strcmp(Lethal.Jsl(iLeth),model.rxns));
        modeldel.lb(del_Idx)=0;
        modeldel.ub(del_Idx)=0;
        solKO_i=optimizeCbModel(modeldel);
        if solKO_i.f<0.01*grWT
           metAll=find(modeldel.S(:,del_Idx));
           met=[];
           if modeldel.rev(del_Idx)==0
               for iMet=1:length(metAll)
                   if modeldel.S(metAll(iMet),del_Idx)>0
                       met=[met;modeldel.mets(metAll(iMet))];
                   end
               end
           else    
               met=modeldel.mets(find(modeldel.S(:,del_Idx)));
           end 
           for jMet=1:length(met)
               modelExc=addExchangeRxn(modeldel,met(jMet));
                  solRM_ij=optimizeCbModel(modelExc,'max','one');
                   if(solRM_ij.f>0.01*grWT)
                       ResMet.SL=[ResMet.SL;model.rxns(del_Idx),met(jMet),solRM_ij.f];
                   end
           end
        end
        modeldel.lb(del_Idx)=model.lb(del_Idx);
        modeldel.ub(del_Idx)=model.ub(del_Idx);
    end
    
    for iLeth=1:length(Lethal.Jdl)
        del_Idx1= find(strcmp(Lethal.Jdl(iLeth,1),model.rxns));
        del_Idx2= find(strcmp(Lethal.Jdl(iLeth,2),model.rxns));
        modeldel.lb(del_Idx1)=0;
        modeldel.ub(del_Idx1)=0;
        modeldel.lb(del_Idx2)=0;
        modeldel.ub(del_Idx2)=0;
        solKO_i=optimizeCbModel(modeldel);
        if solKO_i.f<0.01*grWT
           metAll1=find(modeldel.S(:,del_Idx1));
           met1=[];
           if modeldel.rev(del_Idx1)==0
               for iMet=1:length(metAll1)
                   if modeldel.S(metAll1(iMet),del_Idx1)>0
                       met1=[met1;modeldel.mets(metAll1(iMet))];
                   end
               end
           else    
               met1=modeldel.mets(find(modeldel.S(:,del_Idx1)));
           end 
           metAll2=find(modeldel.S(:,del_Idx2));
           met2=[];
           if modeldel.rev(del_Idx2)==0
               for iMet=1:length(metAll2)
                   if modeldel.S(metAll2(iMet),del_Idx2)>0
                       met2=[met2;modeldel.mets(metAll2(iMet))];
                   end
               end
           else    
               met2=modeldel.mets(find(modeldel.S(:,del_Idx2)));
           end 
           met=unique([met1;met2]);
           
           for jMet=1:length(met)
               modelExc=addExchangeRxn(modeldel,met(jMet));
                  solRM_ij=optimizeCbModel(modelExc,'max','one');
                   if(solRM_ij.f>0.01*grWT)
                       ResMet.DL=[ResMet.DL;model.rxns(del_Idx1),model.rxns(del_Idx2),met(jMet),solRM_ij.f];
                   end
           end
        end
        modeldel.lb(del_Idx1)=model.lb(del_Idx1);
        modeldel.ub(del_Idx1)=model.ub(del_Idx1);
        modeldel.lb(del_Idx2)=model.lb(del_Idx2);
        modeldel.ub(del_Idx2)=model.ub(del_Idx2);
    end
    
    for iLeth=1:length(Lethal.Jtl)
        del_Idx1= find(strcmp(Lethal.Jtl(iLeth,1),model.rxns));
        del_Idx2= find(strcmp(Lethal.Jtl(iLeth,2),model.rxns));
        del_Idx3= find(strcmp(Lethal.Jtl(iLeth,3),model.rxns));
        modeldel.lb(del_Idx1)=0;
        modeldel.ub(del_Idx1)=0;
        modeldel.lb(del_Idx2)=0;
        modeldel.ub(del_Idx2)=0;
        modeldel.lb(del_Idx3)=0;
        modeldel.ub(del_Idx3)=0;
        solKO_i=optimizeCbModel(modeldel);
        if solKO_i.f<0.01*grWT
           metAll1=find(modeldel.S(:,del_Idx1));
           met1=[];
           if modeldel.rev(del_Idx1)==0
               for iMet=1:length(metAll1)
                   if modeldel.S(metAll1(iMet),del_Idx1)>0
                       met1=[met1;modeldel.mets(metAll1(iMet))];
                   end
               end
           else    
               met1=modeldel.mets(find(modeldel.S(:,del_Idx1)));
           end 
           metAll2=find(modeldel.S(:,del_Idx2));
           met2=[];
           if modeldel.rev(del_Idx2)==0
               for iMet=1:length(metAll2)
                   if modeldel.S(metAll2(iMet),del_Idx2)>0
                       met2=[met2;modeldel.mets(metAll2(iMet))];
                   end
               end
           else    
               met2=modeldel.mets(find(modeldel.S(:,del_Idx2)));
           end 
           metAll3=find(modeldel.S(:,del_Idx3));
           met3=[];
           if modeldel.rev(del_Idx3)==0
               for iMet=1:length(metAll3)
                   if modeldel.S(metAll3(iMet),del_Idx3)>0
                       met3=[met3;modeldel.mets(metAll3(iMet))];
                   end
               end
           else    
               met3=modeldel.mets(find(modeldel.S(:,del_Idx3)));
           end 
           met=unique([met1;met2;met3]);
           
           for jMet=1:length(met)
               modelExc=addExchangeRxn(modeldel,met(jMet));
                  solRM_ij=optimizeCbModel(modelExc,'max','one');
                   if(solRM_ij.f>0.01*grWT)
                       ResMet.TL=[ResMet.TL;model.rxns(del_Idx1),model.rxns(del_Idx2),model.rxns(del_Idx3),met(jMet),solRM_ij.f];
                   end
           end
        end
        modeldel.lb(del_Idx1)=model.lb(del_Idx1);
        modeldel.ub(del_Idx1)=model.ub(del_Idx1);
        modeldel.lb(del_Idx2)=model.lb(del_Idx2);
        modeldel.ub(del_Idx2)=model.ub(del_Idx2);
        modeldel.lb(del_Idx3)=model.lb(del_Idx3);
        modeldel.ub(del_Idx3)=model.ub(del_Idx3);
    end
end 
           
end