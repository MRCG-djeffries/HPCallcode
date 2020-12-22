function leout=life_expectancy(z)
% from GBD 2017
% http://ghdx.healthdata.org/record/global-burden-disease-study-2017-gbd-2017-reference-life-table
y=[0 87.885872
1	87.007248
5	83.035378
10	78.050774
15	73.069237
20	68.110138
25	63.157372
30	58.207291
35	53.27124
40	48.368408
45	43.49641
50	38.703121
55	33.98209
60	29.31563
65	24.73456
70	20.32095
75	16.09445
80	12.18093
85	8.7796783
90	6.0613198
95	3.8977709
100	2.2286451
105	1.6117361
110	1.363304];
% yy=interp1(y(:,1),y(:,2),20:1/2:90);
% figure;plot(20:1/2:90,yy)
% grid on
% xlabel('age')
% ylabel('GBD LE')
% z=[22.5 27.5 32.5 40 50 60 70 80 90];
leout=round(interp1(y(:,1),y(:,2),z));


