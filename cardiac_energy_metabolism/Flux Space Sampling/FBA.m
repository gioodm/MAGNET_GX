initCobraToolbox(false)
changeCobraSolver('gurobi','all');

GIMME_modelControlFemale = changeObjective(GIMME_modelControlFemale, 'DM_atp_c_');
GIMME_modelDCMFemale = changeObjective(GIMME_modelDCMFemale, 'DM_atp_c_');
GIMME_modelControlMale = changeObjective(GIMME_modelControlMale, 'DM_atp_c_');
GIMME_modelDCMMale = changeObjective(GIMME_modelDCMMale, 'DM_atp_c_');

[P,GIMME_modelDCMMale] = chrrParseModel(GIMME_modelDCMMale);

[samples, roundedPolytope] = chrrSampler(GIMME_modelDCMMale, 100000, 100);

save('samplesDCMMale100.mat', 'samples')
save('roundedPolytopeDCMMale100.mat', 'roundedPolytope')
