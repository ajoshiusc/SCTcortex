clc;clear all;close all;

restoredefaultpath;
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab\dev'));
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab\src'));

subname1='E:\sipi_data\oasis_full_data\data\\OAS1_0001_MR1\\PROCESSED\\MPRAGE\\T88_111\\OAS1_0001_MR1_aaj';
%subname2='C:\Users\ajoshi\Downloads\oasis_data\OAS1_0061_MR2\OAS1_0061_MR2_mpr_n4_anon_sbj_111.RAS';

NRingsNbr = 16
sub_SCT_hist (subname1, NRingsNbr);
%sub_SCT (subname1, 5);

