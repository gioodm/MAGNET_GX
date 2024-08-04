initCobraToolbox(false)
changeCobraSolver('glpk','all');

%load the model
load('Recon3D_301.mat');

load('GIMME_modelControlFemale');
load('GIMME_modelControlMale');
load('GIMME_modelDCMFemale');
load('GIMME_modelDCMMale');

GIMME_modelControlFemale = changeObjective(GIMME_modelControlFemale, 'DM_atp_c_');
GIMME_modelControlMale = changeObjective(GIMME_modelControlMale, 'DM_atp_c_');
GIMME_modelDCMFemale = changeObjective(GIMME_modelDCMFemale, 'DM_atp_c_');
GIMME_modelDCMMale = changeObjective(GIMME_modelDCMMale, 'DM_atp_c_');

modelfva_controlF = GIMME_modelControlMale; %load model Control
modelfva_controlM = GIMME_modelControlFemale; 
modelfva_DCMF = GIMME_modelDCMFemale; %load model DCM
modelfva_DCMM = GIMME_modelDCMMale;

modelfva_controlF.subSystems=string(modelfva_controlF.subSystems);
modelfva_controlM.subSystems=string(modelfva_controlM.subSystems);
modelfva_DCMF.subSystems=string(modelfva_DCMF.subSystems);
modelfva_DCMM.subSystems=string(modelfva_DCMM.subSystems);

% Selecting several reactions of the model that we want to analyse with FVA
subSystem_fao = string(Recon3D.subSystems(169)); % fatty acid oxidation, 17 rxns
subSystem_gg = string(Recon3D.subSystems(233)); % glycolysis/gluconeogenesis, no rxns
subSystem_gm = string(Recon3D.subSystems(156)); % glutamate metabolism, no rxns
subSystem_cac = string(Recon3D.subSystems(203)); % citric acid cycle, no rxns
subSystem_pm = string(Recon3D.subSystems(239)); % pyruvate metabolism, 5 rxns
subSystem_op = string(Recon3D.subSystems(2755)); % oxidative phosphorylation, 1 rxns

rxnsListCF = modelfva_controlF.rxns(ismember(modelfva_controlF.subSystems,subSystem_pm));
rxnsListCM = modelfva_controlM.rxns(ismember(modelfva_controlM.subSystems,subSystem_pm));
rxnsListDCMF = modelfva_DCMF.rxns(ismember(modelfva_DCMF.subSystems,subSystem_pm));
rxnsListDCMM = modelfva_DCMM.rxns(ismember(modelfva_DCMM.subSystems,subSystem_pm));

% % %remove the genes only present in one set or the other
% % %rxn only in Control
% rxnsListdiffpresCnotD= setdiff(rxnsListCF, rxnsListDCMF);
% % %rxn only in DCM
% rxnsListdiffnotprescontrol= setdiff(rxnsListDCMF, rxnsListCF);
% %remove the unique rxns
% % %normally we should have 2 sets of identical rxn, using one or the other
% % %CF and DF should be the same 
% rxnsListCF2= setdiff(rxnsListCF, rxnsListdiffpresCnotD);
% rxnsListDF2= setdiff(rxnsListDCMF, rxnsListdiffnotprescontrol);

allrxnsList = intersect(intersect(intersect(rxnsListDCMM, rxnsListDCMF), rxnsListCM), rxnsListCF);

% Run FVA analysis for the model with the constraints that simulates aerobic conditions:
[minFlux1, maxFlux1] = fluxVariability(modelfva_controlF, 90, 'max', allrxnsList);
[minFlux2, maxFlux2] = fluxVariability(modelfva_controlM, 90, 'max', allrxnsList);
[minFlux3, maxFlux3] = fluxVariability(modelfva_DCMF, 90, 'max', allrxnsList);
[minFlux4, maxFlux4] = fluxVariability(modelfva_DCMM, 90, 'max', allrxnsList);

diffFluxesmin = zeros(length(minFlux1),5);
j=1;

for i = 1:length(minFlux1)-1
    if round(minFlux1(i)) == round(minFlux2(i)) & round(minFlux3(i)) ~= round(minFlux4(i))
       diffFluxesmin(j, :) = [i, minFlux1(i), minFlux2(i), minFlux3(i), minFlux4(i)];
       j = j+1;
    elseif round(minFlux1(i)) ~= round(minFlux2(i)) & round(minFlux3(i)) == round(minFlux4(i))
        diffFluxesmin(j, :) = [i, minFlux1(i), minFlux2(i), minFlux3(i), minFlux4(i)];
        j = j+1;
    elseif round(minFlux1(i)) ~= round(minFlux2(i)) & round(minFlux1(i)) ~= round(minFlux3(i))
        diffFluxesmin(j, :) = [i, minFlux1(i), minFlux2(i), minFlux3(i), minFlux4(i)];
        j = j+1;
    end
end

diffFluxesmax = zeros(length(maxFlux1),5);
j=1;

for i = 1:length(maxFlux1)-1
    if round(maxFlux1(i)) == round(maxFlux2(i)) & round(maxFlux3(i)) ~= round(maxFlux4(i))
       diffFluxesmax(j, :) = [i, maxFlux1(i), maxFlux2(i), maxFlux3(i), maxFlux4(i)];
       j = j+1;
    elseif round(maxFlux1(i)) ~= round(maxFlux2(i)) & round(maxFlux3(i)) == round(maxFlux4(i))
        diffFluxesmax(j, :) = [i, maxFlux1(i), maxFlux2(i), maxFlux3(i), maxFlux4(i)];
        j = j+1;
    elseif round(maxFlux1(i)) ~= round(maxFlux2(i)) & round(maxFlux1(i)) ~= round(maxFlux3(i))
        diffFluxesmax(j, :) = [i, maxFlux1(i), maxFlux2(i), maxFlux3(i), maxFlux4(i)];
        j = j+1;
    end
end

rxns = union(diffFluxesmin(:,1), diffFluxesmax(:,1));
rxns(1,:) = [];
allrxnsListsub = allrxnsList(rxns, :);
indexes = findRxnIDs(Recon3D, allrxnsListsub);
rxnNames = Recon3D.rxnNames(indexes);

minFlux1sub = minFlux1(rxns, :);
maxFlux1sub = maxFlux1(rxns, :);
minFlux2sub = minFlux2(rxns, :);
maxFlux2sub = maxFlux2(rxns, :);
minFlux3sub = minFlux3(rxns, :);
maxFlux3sub = maxFlux3(rxns, :);
minFlux4sub = minFlux4(rxns, :);
maxFlux4sub = maxFlux4(rxns, :);

%making the graph
maxf = table(maxFlux1sub, maxFlux2sub, maxFlux3sub, maxFlux4sub);
minf = table(minFlux1sub, minFlux2sub, minFlux3sub, minFlux4sub);

maxfxs = table2cell(maxf);
minfxs = table2cell(minf);
figure
plot1 = bar(cell2mat(maxfxs(1:end, :)));
hold on
plot2 = bar(cell2mat(minfxs(1:end, :)));
hold off
xticklabels(rxnsList);
set(gca, 'XTickLabelRotation', -80);
yticks([-1000 -800 -600 -400 -200 0 200 400 600 800 1000])
xlabel('Reactions from the models')
ylabel('Fluxes')
legend({'ControlF', 'ControlM'}, 'Location', 'southwest')
title('Variations in fluxes in Control and DCM conditions for oxidative phosphorylation)-iMat')

