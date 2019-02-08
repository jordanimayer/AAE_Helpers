%%%%%
% Jordan Mayer
% AAE 590: Space Flight Operations
% HW 01
% January 22, 2019
%
% Create pareto graph using cost and reliability data for various mission
% architecture options for Lunar Lava Tube Explorer mission. See parts
% 1 and 2 of assignment for creation of morphological matrix and generation
% of cost/reliability estimates.
%%%%%

close all;  % close any open figures

% set up cost and reliability matrices (see parts 1 and 2 for creation)
% use 1i (imaginary) to represent nonexistent fields
costs = [0 -5 -5; 0 7 7; 0 3 1i; 0 15 1i; -10 0 0; 0 15 15; 0 5 2; ... 
        0 5 10; -1 0 5; -2 0 100];
rels = [0.99 0.75 0.75; 0.995 0.99 0.98; 0.99 0.995 1i; 0.90 0.85 1i; ...
        0.99 0.98 0.98; 0.995 0.85 0.85; 0.995 0.85 0.88; ...
        0.995 0.99 0.97; 0.98 0.99 0.85; 0.995 0.995 0.90];
    
figure();  % open new figure for plot
hold on;  % overlay plots

% label based on mobility and 
% science payload
nmvl = 'NM, VL';
mvl = 'M, VL';
nmvlir = 'NM, VLIR';
mvlir = 'M, VLIR';
nmvlsr = 'NM, VLSR';
mvlsr = 'M, VLSR';
markSize = 12;

legend_ref = [187 0.8346; 202 0.7130; 192 0.8304; 207 0.7094; ...
              197 0.8137; 212 0.6951; 200 0.8346];
  % precalculated costs for certain mission architectures
  % used only to create legend (to avoid having to process all data points)
  
% create legend
plot(legend_ref(1,2), legend_ref(1,1), '+g');
plot(legend_ref(2,2), legend_ref(2,1), 'og');
plot(legend_ref(3,2), legend_ref(3,1), '+r');
plot(legend_ref(4,2), legend_ref(4,1), 'or');
plot(legend_ref(5,2), legend_ref(5,1), '+k');
plot(legend_ref(6,2), legend_ref(6,1), 'ok');
plot(legend_ref(7,2), legend_ref(7,1), 'vb', 'MarkerSize', markSize);

legend({nmvl, mvl, nmvlir, mvlir, nmvlsr, mvlsr, 'Candidate'});

% analyze every possible combination of options and obtain total
% mission cost and reliability

for opt1 = 1:3
    for opt2 = 1:3
        for opt3 = 1:2
            for opt4 = 1:2
                for opt5 = 1:3
                    for opt6 = 1:3
                        for opt7 = 1:3
                            for opt8 = 1:3
                                for opt9 = 1:3
                                    for opt10 = 1:3
                                        
                                        % ignore infeasible combinations
                                        % (e.g. sample return without
                                        %  sample collection payload)
                                        if (opt7 == 1 && opt9 == 3) || ...
                                           (opt10 == 3 && (opt8 ~= 3 || ...
                                                           opt9 == 1)) || ...
                                           (opt8 == 3 && opt10 ~= 3)
                                            continue;
                                        end
                                        
                                        % estimate reference mission costs
                                        % $200M
                                        mission_cost = ...
                                            200 + ...
                                            costs(1, opt1) + ...
                                            costs(2, opt2) + ...
                                            costs(3, opt3) + ...
                                            costs(4, opt4) + ...
                                            costs(5, opt5) + ...
                                            costs(6, opt6) + ...
                                            costs(7, opt7) + ...
                                            costs(8, opt8) + ...
                                            costs(9, opt9) + ...
                                            costs(10, opt10);
                                        mission_rel = ...
                                            rels(1, opt1) * ...
                                            rels(2, opt2) * ...
                                            rels(3, opt3) * ...
                                            rels(4, opt4) * ...
                                            rels(5, opt5) * ...
                                            rels(6, opt6) * ...
                                            rels(7, opt7) * ...
                                            rels(8, opt8) * ...
                                            rels(9, opt9) * ...
                                            rels(10, opt10);
                                        
                                        if opt8 == 1
                                            if opt6 == 1
                                                style = '+g';
                                                if (imag(legend_ref(1,1)) ~= 0)
                                                    legend_ref(1,1) = mission_cost;
                                                    legend_ref(1,2) = mission_rel;
                                                end
                                            else
                                                style = 'og';
                                                if (imag(legend_ref(2,1)) ~= 0)
                                                    legend_ref(2,1) = mission_cost;
                                                    legend_ref(2,2) = mission_rel;
                                                end
                                            end
                                        elseif opt8 == 2
                                            if opt6 == 1
                                                style = '+r';
                                                if (imag(legend_ref(3,1)) ~= 0)
                                                    legend_ref(3,1) = mission_cost;
                                                    legend_ref(3,2) = mission_rel;
                                                end
                                            else
                                                style = 'or';
                                                if (imag(legend_ref(4,1)) ~= 0)
                                                    legend_ref(4,1) = mission_cost;
                                                    legend_ref(4,2) = mission_rel;
                                                end
                                            end
                                        else
                                            if opt6 == 1
                                                style = '+k';
                                                if (imag(legend_ref(5,1)) ~= 0)
                                                    legend_ref(5,1) = mission_cost;
                                                    legend_ref(5,2) = mission_rel;
                                                end
                                            else
                                                style = 'ok';
                                                if (imag(legend_ref(6,1)) ~= 0)
                                                    legend_ref(6,1) = mission_cost;
                                                    legend_ref(6,2) = mission_rel;
                                                end
                                            end
                                        end
                                        
                                        plot(mission_rel, mission_cost, ...
                                             style);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% format plot
title('Jordan Mayer, AAE 590-SFO, HW 01: Pareto Evaluation of Lava Tube Explorer Concepts');
xlabel('Reliability'); ylabel('Cost ($M)'); grid on;

% calculate cost and reliability for candidate missions
candidate1 = [1 1 1 1 2 1 1 1 2 2];
  % reference mission
candidate2 = [1 1 1 1 2 2 1 1 2 2];
  % mini-rover, enter/exit with crane, VL, X-band
candidate3 = [1 1 1 1 2 2 2 2 1 1];
  % mini-rover, enter with rockets, VLIR, no exit, Ka-band
candidate4 = [1 2 1 1 3 3 3 1 1 1];
  % mini-rover, both rovers treaded, enter with cushioning, no exit, X-band
candidate5 = [1 1 1 1 2 1 1 3 2 3];
  % crane, sample return
cands = [candidate1; candidate2; candidate3; candidate4; candidate5];
% 
% styles = ['*r'; '*g'; '*b'; '*c'; '*k'];
% 
% figure(); hold on;

for k=1:5
    mission_cost = 200 + ...
                   costs(1, cands(k,1)) + ...
                   costs(2, cands(k,2)) + ...
                   costs(3, cands(k,3)) + ...
                   costs(4, cands(k,4)) + ...
                   costs(5, cands(k,5)) + ...
                   costs(6, cands(k,6)) + ...
                   costs(7, cands(k,7)) + ...
                   costs(8, cands(k,8)) + ...
                   costs(9, cands(k,9)) + ...
                   costs(10, cands(k,10));
     mission_rel = rels(1, cands(k,1)) * ...
                   rels(2, cands(k,2)) * ...
                   rels(3, cands(k,3)) * ...
                   rels(4, cands(k,4)) * ...
                   rels(5, cands(k,5)) * ...
                   rels(6, cands(k,6)) * ...
                   rels(7, cands(k,7)) * ...
                   rels(8, cands(k,8)) * ...
                   rels(9, cands(k,9)) * ...
                   rels(10, cands(k,10));
     
     plot(mission_rel, mission_cost, 'vb', 'MarkerSize', markSize);
end