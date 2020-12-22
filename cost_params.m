function [costanual,costoneoff,namesof_anualcost,namesof_oneoff]=cost_params()
                
namesof_anualcost={'c_F0','c_F1','c_F2','c_F3','c_F4','c_DC','c_HCC'};
namesof_oneoff={'c_TF4','c_THCC','c_LT','HCV_diagnosis'};
costanual=[446.7 446.7 446.7 690.85 935.05 15202.43 10759.75] + 34.22;
costoneoff=[565.67 970.20 145565.00 936.63];
