close all
clear all

GSM = GlobalSkyModel([],'haslam');
GSM = GSM.generate(408);
GSM.view(1,true)


