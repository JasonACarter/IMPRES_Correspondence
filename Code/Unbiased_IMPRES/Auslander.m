%Modified version of Auslander et al. "hilclimbing.m"
%https://github.com/noamaus/IMPRES-codes/blob/master/MAIN_CODES/Feature_selection/hilClibming.m

%For ease of implementation, we have changed the input/output of this script:

%Input
%   NEUB.mat, CPall.mat, CFF.mat, and RATS.mat are all data files downloaded from Auslander et al
%   "classifyImmuneCOMP3.m" and "SelectFeaturs.m" are IMPRES algorithm scripts used exactly as done by Auslander et al
%   see https://github.com/noamaus/IMPRES-codes/blob/master/MAIN_CODES/Feature_selection

%   To allow easy output, we have created the "G1.mat" and "G2.mat" files. These contain the gene names for the 756 considered pairs (i.e. export gene names rather than list indices).

%Output
%   Runs 500 iteration and generates new feature sets. Again, the actual IMPRES algorithm has not been changed and features are scored/selected as described by Auslander et al.
%   Final feature set is written to "Genes.txt" file in this folder

%Running this script using the 500 NBL training sets provided by Auslander et al. exactly reproduces the original IMPRES feature set
%As these training sets were supposed to have been generated randomly, this script creates truly random training sets for each iteration. That is, this script runs the IMPRES algorithm as it was originally described

%Array job submission
%   This script runs all steps of the IMPRES algorithm. As it involves 500 hill-climbing iterations and is not optimized, it has a lengthy runtime
%   To run the full algorithm X times, we submitted an array job to the CSHL High Performance Computing Cluster. This script was therefore run independently across X nodes.
%   see "impres_qsub.sh"

%   Alternatively, one feature set can be generated using "run('./Auslander.m')" in MATLAB

load('CFF.mat')
load('CPall.mat')
load('RATS.mat')
load('NEUB.mat')
load('G1.mat')
load('G2.mat')
P8 = find(NEUB.high_risk==0&NEUB.age<=1.5&NEUB.progression==0);
N8 = find(NEUB.high_risk==1&NEUB.age<=1.5&NEUB.progression==1);

reps = 500;
FEATS = zeros(756,reps);

rng('shuffle')

for k = 1:reps

    TSP=randsample(P8,13,false);
    TSN=randsample(N8,13,false);
    TP=randsample(setdiff(P8,TSP),3,false);
    TN=randsample(setdiff(N8,TSN),3,false);

    curAUC = 0.1;
    mAUC = 0;
    round = 1;
    csel=[];

    while round<=15

        for i = 1:length(CFF)
            c2 = unique([csel; CFF(i)]);

            [AUCTR(i)] = classifyImmuneCOMP3(NEUB,CPall,c2,TSP,TSN,0,RATS);
        end
        ni = find(AUCTR==max(AUCTR));
        curAUC = max(AUCTR);
        if mAUC<curAUC
            csel = [csel;CFF(ni)];
            round = round+1;
            mAUC = curAUC;
        else
            break
        end
    end
    [AUCTS(k)] = classifyImmuneCOMP3(NEUB,CPall,csel,TP,TN,0,RATS);
    FEATS(csel,k) = 1;
    STSN(:,k) = TSN;
    STSP(:,k) = TSP;
    STN(:,k) = TN;
    STP(:,k) = TP;
end

RES.STSP=STSP;
RES.STSN=STSN;
RES.STP=STP;
RES.STN=STN;
RES.AUCTS = AUCTS;
RES.FEATS=FEATS;

[FEATS] = SelectFeaturs(RES,CFF);

genes=[];
for i=1:length(FEATS)
genes=[genes,g1{FEATS(i)},'_',g2{FEATS(i)},'","'];
end

fid = fopen('Genes.txt','a');
fprintf(fid, '%s\n', string(genes));
