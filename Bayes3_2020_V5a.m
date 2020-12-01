%Simulator for Jenkins, Samuelson, Penny & Spencer (2020) DF Model
%Simulates the 'suspicious coincidence' effect and its reversal
%
%Model is set up with a,b,c versions for simulating (a) the 'simultaneous'
%presentation condition (Exp1 in the paper), (b) the 'sequential'
%presentataion condition (Exp 2 in the paper), and (c) the generalization
%experiment with 6 sequential exemplars all in the same position (Exp 3 in
%the paper). 
%
%The simulator takes several hours to run through each experiment. To
%optimise this, you can open 3 instances of matlab and run the experiments
%simultaneously. 
%
%The simulator has some plotting options. Currently, it is set to do an
%'endplot' which plots the results and the AIC-BIC comparison. This option
%also saves two files: (1) a .mat file with the percentages for each condition
%as well as the standard error over runs (see data names with 'S' at the
%end); and (2) a space-delimited .prn file with the percentages along with
%the AIC/BIC/log-likelihood values for each model (DFT, Bayes, Dummy).
%
%Note: if you want to look at plots of field activities for each trial, set
%'subjects = 1' and then set 'testplot = 1' below. This will generate 12
%plots (one for each trial). We also included a plot of the 'os_to_ls'
%projection in the paper. You can inspect a variant of this by setting
%'os_to_ls_plot = 1'. This latter option isn't optimized as, in truth, you
%need to set the details of the stimuli / spatial positions to precisely
%replicate our figure from the paper. But this gets close to the spirit of
%the plot for interested users.

%close all;
clear;
rand('state', sum(100*clock));
randn('state', sum(100*clock));

%*******Main Running Mode and Logistical Parameters**********
subjects=        60;      %Number of subjects (trials x 12)
simultaneous=    1;      %Simultaneous presentation = 1, sequential = 0; Exp 3 = 2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endplot = 1; %plot of quantitative results including AIC/BIC comparison
testplot = 0; %plot of results from each subject; only turn on if subjects = 1
os_to_ls_plot = 0; %plot of projection from object to label along MDS dimension

%misc options
randomspace = 1; %randomize spatial positions as in experiment
shufflesupers = 3;

%n_field must be odd
n_field=63;
padding=10; %padding on each side of the actual useful portion of the field.  No items can appear here.  It just stops activity from wrapping around.

%iterating to graph stuff:
os_to_ls_store=zeros(4,(n_field + padding));
peakHistory = zeros(12, 2); %storage space for results.  spots for 8 potential peaks, subject, round, trial, and trial type (single/sub/basic/super) in that order
superCatHistory = zeros(9,3);

%set up the conditions
trialsB = 4;
numRounds = 3;
numLabels = 4;
labelSpacing = (n_field - 60) / numLabels; %numLabels is PER ROUND

%model specific stuff:
center=40; %where the center exemplar shows up in space
spread=10; %spatial spread of exemplars
D2U=1.0; %mapping from dimensions in experiment to units in field -- 1:1 here

init_time=50;
label_event_time=125; %duration of each object presentation
h_time=0;
delay=1;
test_event_time=125;

%RESTING LEVEL
h_os=-4;
h_ls=-31;
h_vos=-6.30;
h_vls=-11.55;

%AMARI DYNAMICS
tau_excite=20.;
tau_inhib=2.5;
beta_os=5;
beta_ls=5;
beta_os_to_ls=0.8;

%OBJECT SHAPE (SPACE SIMILARITY) PARAMETERS
osos_sigma=3;
osossim_sigma=3;
w_osos=6.0;

vos_sigma=15;
vossim_sigma=15;
w_vos=40.0;

osv_sigma=15;
osvsim_sigma=15;
w_osv=40.0;
w_osv_const=0.06;

%LABEL SIMILARITY PARAMETERS
lsls_sigma=1;
lslssim_sigma=1;
w_lsls=4;

vls_sigma=6;
vlssim_sigma=6;
w_vls=22.5140;

lsv_sigma=60;
lsvsim_sigma=60;
w_lsv=29.1037;
w_lsv_const=.004;

%PRESHAPE DYNAMICS: LABEL FIELD
lspre_sigma=2;
w_lspre= 4.0065;
taupre_build=4500;
taupre_decay=50000;

%SIMILARITY DIMENSION COUPLING
os_to_ls_sigma=1;
w_os_to_ls= 0.6;

ls_to_os_sigma=10;
w_ls_to_os=0.01;

%STOCHASTICS
noise_sigma= 4;
noise_strength= 3;

%Similarity Structure Preshape
ct_pre_sim=[1];
pre_sim_pos=[0];
pre_sim_strength=[0];
pre_sim_width=[30];

%LABEL RIDGE INPUTS: LABELING AND TEST EVENTS
ct_label_ridge=[1];
label_ridge_strength=8;
label_ridge_width=5;

%TEST VARIABLES
test_peak_threshold = 0;
obj_space_width = 3.25;
obj_space_strength=11.1;
space_input_widen=1;
ct_obj_space_input_test=24;
obj_space_strength_test=16;
obj_space_width_test=1.75;

%SET PROPERTIES
sub_dist=1;
basic_dist=6;
super_dist=18;

for runtimes = 1:subjects
    
    runtimes
    
    %Construct superordinate sets
    if shufflesupers < 2
        allSuperCats = padding+[1 2 3 6 8 18 22 3 2 10 8 14 16 18 22;1 2 3 6 8 18 22 3 2 10 8 14 16 18 22;1 2 3 6 8 18 22 3 2 10 8 14 16 18 22];
    elseif shufflesupers == 2
        allSuperCats(1,:) = padding + [20 20 20 81 73 215 183 20 38 73 103 196 149 166 241]; %animals
        allSuperCats(2,:) = padding + [133 133 133 105 165 218 241 139 125 20 99 172 212 189 146];%vehicles
        allSuperCats(3,:) = padding + [79 79 79 42 60 176 161 79 79 20 58 240 177 106 241];%vegetables
    else
        allSuperCats(1,:) = padding + treeSolve(1,1,n_field,1); %animals
        allSuperCats(2,:) = padding + treeSolve(1,2,n_field,1); %vehicles
        allSuperCats(3,:) = padding + treeSolve(1,3,n_field,1); %vegetables
    end
    
    %store test item positions for accurate plotting later.
    superCatHistoryT = zeros(9,3);
    superCatHistoryT(1:8,1) = allSuperCats(1,8:15);
    superCatHistoryT(1:8,2) = allSuperCats(2,8:15);
    superCatHistoryT(1:8,3) = allSuperCats(3,8:15);
    superCatHistoryT(9,:) = [runtimes runtimes runtimes];
    if runtimes == 1
        superCatHistory = superCatHistoryT;
    else
        superCatHistory = [superCatHistory superCatHistoryT];
    end
    
    %Construct randomized label list, such that labels do not overlap, though
    labelList = [padding + 30:labelSpacing:padding + n_field-30];
    labelList = randswap(labelList);
    
    toa = -1;
    for currentRound = 1:numRounds
        subSimPos = mean(allSuperCats(currentRound,8:9));
        workingSet = allSuperCats(currentRound,:);
        toa = [2 3 4];
        toa = randswap(toa);
        toa = [1 toa];
        
        for (trial = 1:trialsB)
            
            label_ridge_pos = labelList(trial);
            
            %SPATIAL SHUFFLE - In behavioral studies, this was done every trial
            if randomspace < 2
                spatialOptions = [padding + 20:((n_field-40)/(7)):(padding + 20 + 7*(n_field-40)/(7))]; %in order from 20 to whatever
                spatialOptions = round(spatialOptions);
            else
                spatialOptions = [20 41 83 124 165 207 41 58 92 131 170 211 83 92 117 149 185 223 124 131 149 175 207 241]; %all possible euclidean distances from the corner of the array.  Basically collapsing along the diagonal, which is the longest axis, for minimum distortion.
                spatialOptions = randswap(spatialOptions); %shuffle
                spatialOptions = spatialOptions(1:8); %cut off all but 8 spots.
            end
            workingSet(16:23) = spatialOptions;
            %END SPATIAL SHUFFLE
            
            %get an extra set of stimuli for 6 subordinate examples
            %case for experiment 3
            if (simultaneous == 2 && toa(trial) == 2)
                subs = [workingSet(1:3) workingSet(8:9) workingSet(16:17)];
                bass = [workingSet(4:5) workingSet(10:11) workingSet(18:19)];
                sups = [workingSet(6:7) workingSet(12:15) workingSet(20:23)];
                
                subs = randswap(subs,'full');
                bass = randswap(bass,'full');
                sups = randswap(sups,'full');
                
                workingSet2 = [subs(1:3) bass(1:2) sups(1:2) subs(4:5) bass(3:4) sups(3:6) subs(6:7) bass(5:6) sups(7:10)];
            end
            
            %Now shuffle all items of like kinds (subs with subs, bas
            %with bas, sup with sup:
            subs = [workingSet(1:3) workingSet(8:9) workingSet(16:17)];
            bass = [workingSet(4:5) workingSet(10:11) workingSet(18:19)];
            sups = [workingSet(6:7) workingSet(12:15) workingSet(20:23)];
            
            subs = randswap(subs,'full');
            bass = randswap(bass,'full');
            sups = randswap(sups,'full');
            
            workingSet = [subs(1:3) bass(1:2) sups(1:2) subs(4:5) bass(3:4) sups(3:6) subs(6:7) bass(5:6) sups(7:10)];
            
            %ONE SUBORDINATE EXAMPLE
            if toa(trial) == 1
                ct_obj_space_input=[1];
                space_pos=[center]; %used to be +padding, and 50,30,70.  getting rid of that and using low numbers means this happens IN the padding area.
                obj_pos=[workingSet(1)];
                obj_space_strength1=[obj_space_strength];
                obj_space_strength2=[obj_space_strength];
                obj_space_strength3=[obj_space_strength];
                obj_space_strength4=[obj_space_strength];
                obj_space_strength5=[obj_space_strength];
                obj_space_strength6=[obj_space_strength];
                space_pos_test= [workingSet(16),workingSet(17),workingSet(18),workingSet(19),workingSet(20),workingSet(21),workingSet(22),workingSet(23)];
                obj_pos_test=  [workingSet(8),workingSet(12),workingSet(10),workingSet(13),workingSet(9),workingSet(14),workingSet(15),workingSet(11)];
                %BASIC LEVEL CONDITION
            elseif toa(trial) == 3
                ct_obj_space_input=[3];
                space_pos=[center,center-spread,center+spread];
                obj_pos=[workingSet(1),workingSet(4),workingSet(5)];
                obj_space_strength1=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength2=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength3=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength4=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength5=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength6=[obj_space_strength,obj_space_strength,obj_space_strength];
                space_pos_test=[workingSet(16),workingSet(17),workingSet(18),workingSet(19),workingSet(20),workingSet(21),workingSet(22),workingSet(23)];
                obj_pos_test=  [workingSet(8),workingSet(12),workingSet(10),workingSet(13),workingSet(9),workingSet(14),workingSet(15),workingSet(11)];
                %SUPERORDINATE LEVEL CONDITION
            elseif toa(trial) == 4
                ct_obj_space_input=[3];
                space_pos=[center,center-spread,center+spread];
                obj_pos=[workingSet(1),workingSet(6),workingSet(7)];
                obj_space_strength1=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength2=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength3=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength4=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength5=[obj_space_strength,obj_space_strength,obj_space_strength];
                obj_space_strength6=[obj_space_strength,obj_space_strength,obj_space_strength];
                space_pos_test=[workingSet(16),workingSet(17),workingSet(18),workingSet(19),workingSet(20),workingSet(21),workingSet(22),workingSet(23)];
                obj_pos_test=  [workingSet(8),workingSet(12),workingSet(10),workingSet(13),workingSet(9),workingSet(14),workingSet(15),workingSet(11)];
            end
            
            if (simultaneous == 1)
                %THREE SUBORDINATE EXAMPLE
                if toa(trial) == 2
                    ct_obj_space_input=[3];
                    space_pos=[center,center-spread,center+spread];
                    obj_pos=[workingSet(1),workingSet(2),workingSet(3)];
                    obj_space_strength1=[obj_space_strength,obj_space_strength,obj_space_strength];
                    obj_space_strength2=[obj_space_strength,obj_space_strength,obj_space_strength];
                    obj_space_strength3=[obj_space_strength,obj_space_strength,obj_space_strength];
                    obj_space_strength4=[obj_space_strength,obj_space_strength,obj_space_strength];
                    obj_space_strength5=[obj_space_strength,obj_space_strength,obj_space_strength];
                    obj_space_strength6=[obj_space_strength,obj_space_strength,obj_space_strength];
                    space_pos_test=[workingSet(16),workingSet(17),workingSet(18),workingSet(19),workingSet(20),workingSet(21),workingSet(22),workingSet(23)];
                    obj_pos_test=  [workingSet(8),workingSet(12),workingSet(10),workingSet(13),workingSet(9),workingSet(14),workingSet(15),workingSet(11)];
                end
            elseif (simultaneous == 0) %sequential
                %THREE SUBORDINATE EXAMPLE
                if toa(trial) == 2
                    ct_obj_space_input=[3];
                    space_pos=[center,center-spread,center+spread];
                    obj_pos=[workingSet(1),workingSet(2),workingSet(3)];
                    obj_space_strength1=[obj_space_strength,0,0];
                    obj_space_strength2=[0,obj_space_strength,0];
                    obj_space_strength3=[0,0,obj_space_strength];
                    obj_space_strength4=[obj_space_strength,0,0];
                    obj_space_strength5=[0,obj_space_strength,0];
                    obj_space_strength6=[0,0,obj_space_strength];
                    space_pos_test=[workingSet(16),workingSet(17),workingSet(18),workingSet(19),workingSet(20),workingSet(21),workingSet(22),workingSet(23)];
                    obj_pos_test=  [workingSet(8),workingSet(12),workingSet(10),workingSet(13),workingSet(9),workingSet(14),workingSet(15),workingSet(11)];
                end
            elseif (simultaneous == 2) %exp 3
                %THREE SUBORDINATE EXAMPLE
                if toa(trial) == 2
                    ct_obj_space_input=[6];
                    space_pos=[center, center, center, center, center, center];
                    obj_pos=[workingSet(1),workingSet(2),workingSet(3),workingSet2(1),workingSet2(2),workingSet2(3)];
                    obj_space_strength1=[obj_space_strength,0,0,0,0,0];
                    obj_space_strength2=[0,obj_space_strength,0,0,0,0];
                    obj_space_strength3=[0,0,obj_space_strength,0,0,0];
                    obj_space_strength4=[0,0,0,obj_space_strength,0,0];
                    obj_space_strength5=[0,0,0,0,obj_space_strength,0];
                    obj_space_strength6=[0,0,0,0,0,obj_space_strength];
                    space_pos_test=[workingSet(16),workingSet(17),workingSet(18),workingSet(19),workingSet(20),workingSet(21),workingSet(22),workingSet(23)];
                    obj_pos_test=  [workingSet(8),workingSet(12),workingSet(10),workingSet(13),workingSet(9),workingSet(14),workingSet(15),workingSet(11)];
                end
            end
            
            if toa(trial) == 2
                test_event_time2 = test_event_time;
                label_event_time2 = label_event_time;
            else
                test_event_time2 = 1;
                label_event_time2 = 1;
            end
            
            %exp 3 has 6 stimuli * 8 tests = 48 instead of 3 stim * 8 tests = 24 events
            %keep initial 6 events for 54 total
            %adding new data structures so 3sub trials treated differently.
            %that way, all experiments comparable...
            n_time = init_time+...
                label_event_time+label_event_time2+label_event_time+label_event_time2+label_event_time+label_event_time2+ ...
                h_time+delay+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2+ ...
                test_event_time+test_event_time2+test_event_time+test_event_time2;
            
            n_events=57;
            
            event_time = cumsum([init_time ...
                label_event_time label_event_time2 label_event_time label_event_time2 label_event_time label_event_time2 ...
                h_time delay ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2]);
            
            %added full OS field resets
            %reset after 6 iterations
            reset_os_time = [event_time(7) event_time(15) event_time(21) event_time(27) ...
                event_time(33) event_time(39) event_time(45) ...
                event_time(51) event_time(57)];
            
            reset_ls_time=[event_time(7) event_time(15) event_time(21) event_time(27) ...
                event_time(33) event_time(39) event_time(45) ...
                event_time(51) event_time(57)]; %add on all such times as an array here.
            
            probe_ct=13;
            
            probe_time=cumsum([init_time ... %ctps are based on this.  1 = after the first label event
                label_event_time+label_event_time2 label_event_time+label_event_time2 label_event_time+label_event_time2 ...
                h_time+delay 3*test_event_time+3*test_event_time2 3*test_event_time+3*test_event_time2 ...
                3*test_event_time+3*test_event_time2 3*test_event_time+3*test_event_time2 ...
                3*test_event_time+3*test_event_time2 3*test_event_time+3*test_event_time2 ...
                3*test_event_time+3*test_event_time2 3*test_event_time+3*test_event_time2]);
            
            probe_ct_ilsos=48;
            
            probe_time_ilsos=cumsum([init_time+ ... %ctps are based on this.  1 = after the first label event
                label_event_time++label_event_time2+label_event_time++label_event_time2+label_event_time++label_event_time2+...
                h_time+delay+ ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2 ...
                test_event_time test_event_time2 test_event_time test_event_time2]);
            
            % input kernels
            pre_sim_input_t=zeros(1,padding+n_field);
            for ct=1:ct_pre_sim
                pre_sim_input_t=pre_sim_input_t+gauss(padding+n_field-1,(fix(pre_sim_pos(ct)*D2U))-1,pre_sim_strength(ct),pre_sim_width(ct),0);
            end
            pre_sim_input=(ones(padding+n_field,1)*pre_sim_input_t)';
            
            %create object_space inputs for each presentation cycle
            obj_space_input1=zeros(padding+n_field,padding+n_field);
            for ct=1:ct_obj_space_input
                obj_space_input1=obj_space_input1 + obj_space_strength1(ct)*gauss(padding+n_field-1,(fix(space_pos(ct)*D2U))-1,1.0,(obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2)*space_input_widen,0)'...
                    *gauss(padding+n_field-1,(fix(obj_pos(ct)*D2U))-1,1.0,obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2,0);
            end
            obj_space_input2=zeros(padding+n_field,padding+n_field);
            for ct=1:ct_obj_space_input
                obj_space_input2=obj_space_input2 + obj_space_strength2(ct)*gauss(padding+n_field-1,(fix(space_pos(ct)*D2U))-1,1.0,(obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2)*space_input_widen,0)'...
                    *gauss(padding+n_field-1,(fix(obj_pos(ct)*D2U))-1,1.0,obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2,0);
            end
            obj_space_input3=zeros(padding+n_field,padding+n_field);
            for ct=1:ct_obj_space_input
                obj_space_input3=obj_space_input3 + obj_space_strength3(ct)*gauss(padding+n_field-1,(fix(space_pos(ct)*D2U))-1,1.0,(obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2)*space_input_widen,0)'...
                    *gauss(padding+n_field-1,(fix(obj_pos(ct)*D2U))-1,1.0,obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2,0);
            end
            obj_space_input4=zeros(padding+n_field,padding+n_field);
            for ct=1:ct_obj_space_input
                obj_space_input4=obj_space_input4 + obj_space_strength4(ct)*gauss(padding+n_field-1,(fix(space_pos(ct)*D2U))-1,1.0,(obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2)*space_input_widen,0)'...
                    *gauss(padding+n_field-1,(fix(obj_pos(ct)*D2U))-1,1.0,obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2,0);
            end
            obj_space_input5=zeros(padding+n_field,padding+n_field);
            for ct=1:ct_obj_space_input
                obj_space_input5=obj_space_input5 + obj_space_strength5(ct)*gauss(padding+n_field-1,(fix(space_pos(ct)*D2U))-1,1.0,(obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2)*space_input_widen,0)'...
                    *gauss(padding+n_field-1,(fix(obj_pos(ct)*D2U))-1,1.0,obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2,0);
            end
            obj_space_input6=zeros(padding+n_field,padding+n_field);
            for ct=1:ct_obj_space_input
                obj_space_input6=obj_space_input6 + obj_space_strength6(ct)*gauss(padding+n_field-1,(fix(space_pos(ct)*D2U))-1,1.0,(obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2)*space_input_widen,0)'...
                    *gauss(padding+n_field-1,(fix(obj_pos(ct)*D2U))-1,1.0,obj_space_width+(toa(trial)==3)+(toa(trial)==4)*2,0);
            end
            
            label_ridge_input_t=zeros(1,padding+n_field);
            for ct=1:ct_label_ridge
                label_ridge_input_t=label_ridge_input_t+gauss(padding+n_field-1,(fix(label_ridge_pos(ct)*D2U))-1,label_ridge_strength(ct),label_ridge_width(ct),0);
            end
            label_ridge_input=(ones(padding+n_field,1)*label_ridge_input_t)';
            
            %ridges used instead here:
            test_ridge_input_t = zeros(1,padding+n_field,ct_obj_space_input_test);
            test_ridge_input = zeros(padding+n_field,padding+n_field,ct_obj_space_input_test);
            for ct=1:ct_obj_space_input_test
                dt = ct;
                if ct>8
                    dt = dt - 8;
                end
                if ct > 16
                    dt = dt - 8;
                end
                test_ridge_input_t(1,:,dt) = gauss(padding+n_field-1,(fix(obj_pos_test(dt)*D2U))-1,obj_space_strength_test,obj_space_width_test,0);
                for i=1:padding+n_field
                    test_ridge_input(i,:,dt) = test_ridge_input_t(1,:,dt);
                end
            end
            
            % Input to Space-Similarity Field
            % NOTE: x=space, y=similarity
            s_stimulus=zeros(n_events,padding+n_field,padding+n_field);
            s_stimulus(1,:,:)=pre_sim_input;
            s_stimulus(2,:,:)=pre_sim_input + obj_space_input1;
            s_stimulus(3,:,:)=pre_sim_input + obj_space_input4;
            s_stimulus(4,:,:)=pre_sim_input + obj_space_input2;
            s_stimulus(5,:,:)=pre_sim_input + obj_space_input5;
            s_stimulus(6,:,:)=pre_sim_input + obj_space_input3;
            s_stimulus(7,:,:)=pre_sim_input + obj_space_input6;
            s_stimulus(8,:,:)=pre_sim_input;
            s_stimulus(9,:,:)=pre_sim_input;
            s_stimulus(10,:,:)=pre_sim_input + obj_space_input1;
            s_stimulus(11,:,:)=pre_sim_input + obj_space_input4;
            s_stimulus(12,:,:)=pre_sim_input + obj_space_input2;
            s_stimulus(13,:,:)=pre_sim_input + obj_space_input5;
            s_stimulus(14,:,:)=pre_sim_input + obj_space_input3;
            s_stimulus(15,:,:)=pre_sim_input + obj_space_input6;
            s_stimulus(16,:,:)=pre_sim_input + obj_space_input1;
            s_stimulus(17,:,:)=pre_sim_input + obj_space_input4;
            s_stimulus(18,:,:)=pre_sim_input + obj_space_input2;
            s_stimulus(19,:,:)=pre_sim_input + obj_space_input5;
            s_stimulus(20,:,:)=pre_sim_input + obj_space_input3;
            s_stimulus(21,:,:)=pre_sim_input + obj_space_input6;
            s_stimulus(22,:,:)=pre_sim_input + obj_space_input1;
            s_stimulus(23,:,:)=pre_sim_input + obj_space_input4;
            s_stimulus(24,:,:)=pre_sim_input + obj_space_input2;
            s_stimulus(25,:,:)=pre_sim_input + obj_space_input5;
            s_stimulus(26,:,:)=pre_sim_input + obj_space_input3;
            s_stimulus(27,:,:)=pre_sim_input + obj_space_input6;
            s_stimulus(28,:,:)=pre_sim_input + obj_space_input1;
            s_stimulus(29,:,:)=pre_sim_input + obj_space_input4;
            s_stimulus(30,:,:)=pre_sim_input + obj_space_input2;
            s_stimulus(31,:,:)=pre_sim_input + obj_space_input5;
            s_stimulus(32,:,:)=pre_sim_input + obj_space_input3;
            s_stimulus(33,:,:)=pre_sim_input + obj_space_input6;
            s_stimulus(34,:,:)=pre_sim_input + obj_space_input1;
            s_stimulus(35,:,:)=pre_sim_input + obj_space_input4;
            s_stimulus(36,:,:)=pre_sim_input + obj_space_input2;
            s_stimulus(37,:,:)=pre_sim_input + obj_space_input5;
            s_stimulus(38,:,:)=pre_sim_input + obj_space_input3;
            s_stimulus(39,:,:)=pre_sim_input + obj_space_input6;
            s_stimulus(40,:,:)=pre_sim_input + obj_space_input1;
            s_stimulus(41,:,:)=pre_sim_input + obj_space_input4;
            s_stimulus(42,:,:)=pre_sim_input + obj_space_input2;
            s_stimulus(43,:,:)=pre_sim_input + obj_space_input5;
            s_stimulus(44,:,:)=pre_sim_input + obj_space_input3;
            s_stimulus(45,:,:)=pre_sim_input + obj_space_input6;
            s_stimulus(46,:,:)=pre_sim_input + obj_space_input1;
            s_stimulus(47,:,:)=pre_sim_input + obj_space_input4;
            s_stimulus(48,:,:)=pre_sim_input + obj_space_input2;
            s_stimulus(49,:,:)=pre_sim_input + obj_space_input5;
            s_stimulus(50,:,:)=pre_sim_input + obj_space_input3;
            s_stimulus(51,:,:)=pre_sim_input + obj_space_input6;
            s_stimulus(52,:,:)=pre_sim_input + obj_space_input1;
            s_stimulus(53,:,:)=pre_sim_input + obj_space_input4;
            s_stimulus(54,:,:)=pre_sim_input + obj_space_input2;
            s_stimulus(55,:,:)=pre_sim_input + obj_space_input5;
            s_stimulus(56,:,:)=pre_sim_input + obj_space_input3;
            s_stimulus(57,:,:)=pre_sim_input + obj_space_input6;
            
            % Input to label-Similarity Field
            l_stimulus=zeros(n_events,padding+n_field,padding+n_field);
            l_stimulus(1,:,:)=pre_sim_input;
            l_stimulus(2,:,:)=pre_sim_input + label_ridge_input;
            l_stimulus(3,:,:)=pre_sim_input + label_ridge_input;
            l_stimulus(4,:,:)=pre_sim_input + label_ridge_input;
            l_stimulus(5,:,:)=pre_sim_input + label_ridge_input;
            l_stimulus(6,:,:)=pre_sim_input + label_ridge_input;
            l_stimulus(7,:,:)=pre_sim_input + label_ridge_input;
            l_stimulus(8,:,:)=pre_sim_input;
            l_stimulus(9,:,:)=pre_sim_input;
            l_stimulus(10,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,1);
            l_stimulus(11,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,1);
            l_stimulus(12,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,1);
            l_stimulus(13,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,1);
            l_stimulus(14,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,1);
            l_stimulus(15,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,1);
            l_stimulus(16,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,2);
            l_stimulus(17,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,2);
            l_stimulus(18,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,2);
            l_stimulus(19,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,2);
            l_stimulus(20,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,2);
            l_stimulus(21,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,2);
            l_stimulus(22,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,3);
            l_stimulus(23,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,3);
            l_stimulus(24,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,3);
            l_stimulus(25,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,3);
            l_stimulus(26,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,3);
            l_stimulus(27,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,3);
            l_stimulus(28,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,4);
            l_stimulus(29,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,4);
            l_stimulus(30,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,4);
            l_stimulus(31,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,4);
            l_stimulus(32,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,4);
            l_stimulus(33,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,4);
            l_stimulus(34,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,5);
            l_stimulus(35,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,5);
            l_stimulus(36,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,5);
            l_stimulus(37,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,5);
            l_stimulus(38,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,5);
            l_stimulus(39,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,5);
            l_stimulus(40,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,6);
            l_stimulus(41,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,6);
            l_stimulus(42,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,6);
            l_stimulus(43,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,6);
            l_stimulus(44,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,6);
            l_stimulus(45,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,6);
            l_stimulus(46,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,7);
            l_stimulus(47,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,7);
            l_stimulus(48,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,7);
            l_stimulus(49,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,7);
            l_stimulus(50,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,7);
            l_stimulus(51,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,7);
            l_stimulus(52,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,8);
            l_stimulus(53,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,8);
            l_stimulus(54,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,8);
            l_stimulus(55,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,8);
            l_stimulus(56,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,8);
            l_stimulus(57,:,:)=pre_sim_input + label_ridge_input + test_ridge_input(:,:,8);
            
            % kernels
            kernel_width_multiplier=3; %DO NOT SET TO SMALLER THAN 3
            
            %%%% preparing the kernels
            osos_kernel_width=floor(min(kernel_width_multiplier*osos_sigma,(padding+n_field-1)/2));
            osossim_kernel_width=floor(min(kernel_width_multiplier*osossim_sigma,(padding+n_field-1)/2));
            osos_int_kernel=gaussNorm(-osos_kernel_width:osos_kernel_width,0,osos_sigma);
            osossim_int_kernel=gaussNorm(-osossim_kernel_width:osossim_kernel_width,0,osossim_sigma);
            osos_ext_index=[padding+n_field-osos_kernel_width+1:padding+n_field, 1:padding+n_field, 1: osos_kernel_width];
            osossim_ext_index=[padding+n_field-osossim_kernel_width+1:padding+n_field, 1:padding+n_field, 1: osossim_kernel_width];
            
            osv_kernel_width=floor(min(kernel_width_multiplier*osv_sigma,(padding+n_field-1)/2));
            osvsim_kernel_width=floor(min(kernel_width_multiplier*osvsim_sigma,(padding+n_field-1)/2));
            osv_int_kernel=gaussNorm(-osv_kernel_width:osv_kernel_width,0,osv_sigma);
            osvsim_int_kernel=gaussNorm(-osvsim_kernel_width:osvsim_kernel_width,0,osvsim_sigma);
            osv_ext_index=[padding+n_field-osv_kernel_width+1:padding+n_field, 1:padding+n_field, 1: osv_kernel_width];
            osvsim_ext_index=[padding+n_field-osvsim_kernel_width+1:padding+n_field, 1:padding+n_field, 1: osvsim_kernel_width];
            
            vos_kernel_width=floor(min(kernel_width_multiplier*vos_sigma,(padding+n_field-1)/2));
            vossim_kernel_width=floor(min(kernel_width_multiplier*vossim_sigma,(padding+n_field-1)/2));
            vos_int_kernel=gaussNorm(-vos_kernel_width:vos_kernel_width,0,vos_sigma);
            vossim_int_kernel=gaussNorm(-vossim_kernel_width:vossim_kernel_width,0,vossim_sigma);
            vos_ext_index=[padding+n_field-vos_kernel_width+1:padding+n_field, 1:padding+n_field, 1: vos_kernel_width];
            vossim_ext_index=[padding+n_field-vossim_kernel_width+1:padding+n_field, 1:padding+n_field, 1: vossim_kernel_width];
            
            lsls_kernel_width=floor(min(kernel_width_multiplier*lsls_sigma,(padding+n_field-1)/2));
            lslssim_kernel_width=floor(min(kernel_width_multiplier*lslssim_sigma,(padding+n_field-1)/2));
            lsls_int_kernel=gaussNorm(-lsls_kernel_width:lsls_kernel_width,0,lsls_sigma);
            lslssim_int_kernel=gaussNorm(-lslssim_kernel_width:lslssim_kernel_width,0,lslssim_sigma);
            lsls_ext_index=[padding+n_field-lsls_kernel_width+1:padding+n_field, 1:padding+n_field, 1: lsls_kernel_width];
            lslssim_ext_index=[padding+n_field-lslssim_kernel_width+1:padding+n_field, 1:padding+n_field, 1: lslssim_kernel_width];
            
            lsv_kernel_width=floor(min(kernel_width_multiplier*lsv_sigma,(padding+n_field-1)/2));
            lsvsim_kernel_width=floor(min(kernel_width_multiplier*lsvsim_sigma,(padding+n_field-1)/2));
            lsv_int_kernel=gaussNorm(-lsv_kernel_width:lsv_kernel_width,0,lsv_sigma);
            lsvsim_int_kernel=gaussNorm(-lsvsim_kernel_width:lsvsim_kernel_width,0,lsvsim_sigma);
            lsv_ext_index=[padding+n_field-lsv_kernel_width+1:padding+n_field, 1:padding+n_field, 1: lsv_kernel_width];
            lsvsim_ext_index=[padding+n_field-lsvsim_kernel_width+1:padding+n_field, 1:padding+n_field, 1: lsvsim_kernel_width];
            
            vls_kernel_width=floor(min(kernel_width_multiplier*vls_sigma,(padding+n_field-1)/2));
            vlssim_kernel_width=floor(min(kernel_width_multiplier*vlssim_sigma,(padding+n_field-1)/2));
            vls_int_kernel=gaussNorm(-vls_kernel_width:vls_kernel_width,0,vls_sigma);
            vlssim_int_kernel=gaussNorm(-vlssim_kernel_width:vlssim_kernel_width,0,vlssim_sigma);
            vls_ext_index=[padding+n_field-vls_kernel_width+1:padding+n_field, 1:padding+n_field, 1: vls_kernel_width];
            vlssim_ext_index=[padding+n_field-vlssim_kernel_width+1:padding+n_field, 1:padding+n_field, 1: vlssim_kernel_width];
            
            %Coupling Across 2D Fields
            os_to_ls_kernel_width=floor(min(kernel_width_multiplier*os_to_ls_sigma,(padding+n_field-1)/2));
            os_to_ls_int_kernel=gaussNorm(-os_to_ls_kernel_width:os_to_ls_kernel_width,0,os_to_ls_sigma);
            os_to_ls_ext_index=[padding+n_field-os_to_ls_kernel_width+1:padding+n_field, 1:padding+n_field, 1: os_to_ls_kernel_width];
            
            ls_to_os_kernel_width=floor(min(kernel_width_multiplier*ls_to_os_sigma,(padding+n_field-1)/2));
            ls_to_os_int_kernel=gaussNorm(-ls_to_os_kernel_width:ls_to_os_kernel_width,0,ls_to_os_sigma);
            ls_to_os_ext_index=[padding+n_field-ls_to_os_kernel_width+1:padding+n_field, 1:padding+n_field, 1: ls_to_os_kernel_width];
            
            %LTM
            lspre_kernel_width=floor(min(kernel_width_multiplier*lspre_sigma,(padding+n_field-1)/2));
            lspre_int_kernel=gaussNorm(-lspre_kernel_width:lspre_kernel_width,0,lspre_sigma);
            lspre_ext_index=[padding+n_field-lspre_kernel_width+1:padding+n_field, 1:padding+n_field, 1: lspre_kernel_width];
            
            %Noise
            noise_kernel_width=floor(min(kernel_width_multiplier*noise_sigma,(padding+n_field-1)/2));
            noise_kernel=gaussNorm(-noise_kernel_width:noise_kernel_width,0,noise_sigma);
            noise_size=padding+n_field+2*noise_kernel_width;
            
            %%%%% initialization, initial conditions for this one trial
            ls_field =zeros(padding+n_field,padding+n_field) + h_ls;
            vls_field=zeros(padding+n_field,padding+n_field) + h_vls;
            os_field =zeros(padding+n_field,padding+n_field) + h_os;
            vos_field=zeros(padding+n_field,padding+n_field) + h_vos;
            lspre_field =zeros(padding+n_field,padding+n_field);
            %(perceptual field and STM field should always reset)
            
            os_to_ls_storage =zeros(probe_ct_ilsos,padding+n_field,padding+n_field);
            os_field_storage =zeros(probe_ct,padding+n_field,padding+n_field);
            vos_field_storage=zeros(probe_ct,padding+n_field,padding+n_field);
            ls_field_storage =zeros(probe_ct,padding+n_field,padding+n_field);
            vls_field_storage=zeros(probe_ct,padding+n_field,padding+n_field);
            lspre_field_storage =zeros(probe_ct,padding+n_field,padding+n_field);
            ctp=1;
            ctpilsos = 1;
            
            tic;
            
            peaks = zeros(12, 2);
            peaks(9,:) = runtimes;
            peaks(10,:) = currentRound;
            peaks(11,:) = trial;
            peaks(12,:) = toa(trial);
            bestPeak = -100;
            bestPos = -10;
            
            for time=2:n_time
                
                [maxe,event_ind]=max([1-max(diff(time<=event_time)), diff(time<=event_time)]);
                sS= squeeze(s_stimulus(event_ind, :, :));
                lS= squeeze(l_stimulus(event_ind, :, :));
                
                if (time == probe_time(ctp))
                    os_field_storage(ctp,:,:)=os_field;
                    ls_field_storage(ctp,:,:)=ls_field;
                    vos_field_storage(ctp,:,:)=vos_field;
                    vls_field_storage(ctp,:,:)=vls_field;
                    lspre_field_storage(ctp,:,:)=lspre_field;
                end
                
                if (time == probe_time_ilsos(ctpilsos))
                    os_to_ls_storage(ctpilsos,:,:)=input_os_to_ls;
                end
                
                if ismember(time,reset_os_time)
                    os_field =zeros(padding+n_field,padding+n_field) + h_os ;
                    vos_field=zeros(padding+n_field,padding+n_field) + h_vos ;
                end
                
                if ismember(time,reset_ls_time+1)
                    bestPeak = -100;
                    bestPos = -10;
                end
                [V,I] = max(max(ls_field));
                
                if (V > bestPeak)
                    bestPeak = V;
                    bestPos = I;
                end
                
                if ismember(time,reset_ls_time)
                    %save any sufficiently high peaks and clear label field.
                    ctp;
                    if (bestPeak > test_peak_threshold) && (ctp ~= 4) %There is a peak, and it is not the ls_field reset BEFORE testing...
                        peaks(ctp-5,1) = 1; %12 used to be nPeaks, but this makes it compatible to tack on new data in the save file.  first 8 spots are for any peaks above thresh.  Last 3 spots are for Subject, Round, Trial, Trial type
                        peaks(ctp-5,2) = bestPos; %store the position in similarity dimension of the max peak.
                    end
                    ls_field =zeros(padding+n_field,padding+n_field) + h_ls ;
                    vls_field=zeros(padding+n_field,padding+n_field) + h_vls ;
                end
                
                f_os=sigmoid(os_field,beta_os,0);
                f_os_to_ls=sigmoid(os_field,beta_os_to_ls,0);
                f_vos=sigmoid(vos_field,beta_os,0);
                f_ls=sigmoid(ls_field,beta_ls,0);
                f_vls=sigmoid(vls_field,beta_ls,0);
                f_vos_total = sum(sum(f_vos));
                f_vls_total = sum(sum(f_vls));
                f_ossim=sum(f_os);
                f_os_to_ls_sim=sum(f_os_to_ls);
                f_lssim=sum(f_ls);
                
                %Compute Convolutions
                osos_conv=conv2(1, osossim_int_kernel, f_os(:, osossim_ext_index),'valid');
                osos_conv=conv2(osos_int_kernel, 1, osos_conv(osos_ext_index, :),'valid');
                osvos_conv=conv2(1, osvsim_int_kernel, f_vos(:, osvsim_ext_index),'valid');
                osvos_conv=conv2(osv_int_kernel, 1, osvos_conv(osv_ext_index, :),'valid');
                vosos_conv=conv2(1,vossim_int_kernel,f_os(:,vossim_ext_index),'valid');
                vosos_conv=conv2(vos_int_kernel, 1, vosos_conv(vos_ext_index,:),'valid');
                
                lsls_conv=conv2(1,lslssim_int_kernel,f_ls(:,lslssim_ext_index),'valid');
                lsls_conv=conv2(lsls_int_kernel, 1, lsls_conv(lsls_ext_index,:),'valid');
                lsvls_conv=conv2(1,lsvsim_int_kernel,f_ls(:,lsvsim_ext_index),'valid');
                lsvls_conv=conv2(lsv_int_kernel, 1, lsvls_conv(lsv_ext_index,:),'valid');
                vlsls_conv=conv2(1,vlssim_int_kernel,f_ls(:,vlssim_ext_index),'valid');
                vlsls_conv=conv2(vls_int_kernel, 1, vlsls_conv(vls_ext_index,:),'valid');
                
                os_to_ls_conv=conv2(1,os_to_ls_int_kernel,f_os_to_ls_sim(os_to_ls_ext_index),'valid');
                input_os_to_ls=ones(padding+n_field,1)*(w_os_to_ls*os_to_ls_conv);
                ls_to_os_conv=conv2(1,ls_to_os_int_kernel,f_lssim(ls_to_os_ext_index),'valid');
                input_ls_to_os=ones(padding+n_field,1)*(w_ls_to_os*ls_to_os_conv);
                lspre_conv=conv2(1,lspre_int_kernel,lspre_field(:,lspre_ext_index),'valid');
                lspre_conv=conv2(lspre_int_kernel,1,lspre_conv(lspre_ext_index,:),'valid');
                
                os_noise=conv2(noise_kernel,noise_kernel',randn(noise_size,noise_size),'valid');
                vos_noise=conv2(noise_kernel,noise_kernel',randn(noise_size,noise_size),'valid');
                ls_noise=conv2(noise_kernel,noise_kernel',randn(noise_size,noise_size),'valid');
                vls_noise=conv2(noise_kernel,noise_kernel',randn(noise_size,noise_size),'valid');
                
                % equations
                os_field = os_field + (- os_field + h_os + sS + w_osos*osos_conv - w_osv*osvos_conv - w_osv_const*f_vos_total + input_ls_to_os)/tau_excite + noise_strength*os_noise;
                vos_field = vos_field + (- vos_field + h_vos + w_vos*vosos_conv)/tau_inhib + noise_strength*vos_noise;
                
                ls_field = ls_field + (- ls_field + h_ls  + lS + w_lsls*lsls_conv - w_lsv*lsvls_conv - w_lsv_const*f_vls_total + input_os_to_ls + w_lspre*lspre_conv)/tau_excite + noise_strength*ls_noise;
                vls_field = vls_field + (- vls_field + h_vls + w_vls*vlsls_conv)/tau_inhib + noise_strength*vls_noise;
                
                theta_pre=ls_field>0;
                lspre_field = lspre_field + (theta_pre).*(-lspre_field + f_ls)/taupre_build + (1-theta_pre).*(-lspre_field)/taupre_decay;
                
                % store data every now and then
                if (time == probe_time(ctp))
                    ctp = ctp+1;
                    if (ctp > probe_ct) ctp=probe_ct; end
                end
                
                if (time == probe_time_ilsos(ctpilsos))
                    ctpilsos = ctpilsos + 1;
                end
                
            end %end time loop
            
            toc;
            n_timesteps=n_time;
            time=[1:n_time];
            
            temp_os_ls = zeros(1,n_field+padding);
            for i = 1:probe_ct_ilsos
                temp_os_ls = temp_os_ls + max(squeeze(os_to_ls_storage(i,:,:)));
            end
            temp_os_ls = temp_os_ls / probe_ct_ilsos;
            os_to_ls_store(1,:) = temp_os_ls;
            
            peakHistory = [peakHistory peaks]; %appending new data
            
            if testplot == 1
                EndName = 'Test';
                hFig = figure('Name',EndName);
                set(hFig, 'Position', [20 50 1900 900])
                subplot(3,8,1)
                hold on;
                mesh(squeeze(os_field_storage(6,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(space)');
                zlabel('s(x,y)');
                title('OS Test 1');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,2)
                hold on;
                mesh(squeeze(os_field_storage(7,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(space)');
                zlabel('s(x,y)');
                title('OS Test 2');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                
                subplot(3,8,3)
                hold on;
                mesh(squeeze(os_field_storage(8,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(space)');
                zlabel('s(x,y)');
                title('OS Test 3');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,4)
                hold on;
                mesh(squeeze(os_field_storage(9,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(space)');
                zlabel('s(x,y)');
                title('OS Test 4');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,5)
                hold on;
                mesh(squeeze(os_field_storage(10,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(space)');
                zlabel('s(x,y)');
                title('OS Test 5');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,6)
                hold on;
                mesh(squeeze(os_field_storage(11,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(space)');
                zlabel('s(x,y)');
                title('OS Test 6');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,7)
                hold on;
                mesh(squeeze(os_field_storage(12,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(space)');
                zlabel('s(x,y)');
                title('OS Test 7');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,8)
                hold on;
                mesh(squeeze(os_field_storage(13,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(space)');
                zlabel('s(x,y)');
                title('OS Test 8');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                
                subplot(3,8,9)
                hold on;
                mesh(squeeze(ls_field_storage(6,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LS Test 1');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,10)
                hold on;
                mesh(squeeze(ls_field_storage(7,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LS Test 2');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,11)
                hold on;
                mesh(squeeze(ls_field_storage(8,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LS Test 3');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,12)
                hold on;
                mesh(squeeze(ls_field_storage(9,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LS Test 4');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,13)
                hold on;
                mesh(squeeze(ls_field_storage(10,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LS Test 5');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,14)
                hold on;
                mesh(squeeze(ls_field_storage(11,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(space)');
                zlabel('s(x,y)');
                title('LS Test 6');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,15)
                hold on;
                mesh(squeeze(ls_field_storage(12,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LS Test 7');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,16)
                hold on;
                mesh(squeeze(ls_field_storage(13,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LS Test 8');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                
                subplot(3,8,17)
                hold on;
                mesh(squeeze(lspre_field_storage(6,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LSpre Test 1');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,18)
                hold on;
                mesh(squeeze(lspre_field_storage(7,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LSpre Test 2');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,19)
                hold on;
                mesh(squeeze(lspre_field_storage(8,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LSpre Test 3');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,20)
                hold on;
                mesh(squeeze(lspre_field_storage(9,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LSpre Test 4');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,21)
                hold on;
                mesh(squeeze(lspre_field_storage(10,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LSpre Test 5');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,22)
                hold on;
                mesh(squeeze(lspre_field_storage(11,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LSpre Test 6');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,23)
                hold on;
                mesh(squeeze(lspre_field_storage(12,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LSpre Test 7');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
                subplot(3,8,24)
                hold on;
                mesh(squeeze(lspre_field_storage(13,:,:)));
                plot(obj_pos(1),space_pos(1),'kv');
                xlabel('(sim)');
                ylabel('(label)');
                zlabel('s(x,y)');
                title('LSpre Test 8');
                view(2);
                axis([0 n_field+padding 0 n_field+padding]);
                axis square
                hold off;
            end
            
        end
    end
end

%plotting results numerically for whole run
peakHistory(:,1)=[];
peakHistory(:,1)=[];
peakHistory;
if endplot==1
    graphMatrix = plotBayesForAlphaGamma_V2(peakHistory,0,simultaneous);
end

disp('herpderp')

if (os_to_ls_plot == 1)
    figure
    plot(os_to_ls_store(1,:));
    hold;
    plot(os_to_ls_store(2,:),'g');
    plot(os_to_ls_store(3,:),'r');
    plot(os_to_ls_store(4,:),'k');
    hold off;
end


