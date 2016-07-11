param N;
param M;
set IJ dimen 3;
set I dimen 2;
param Q{IJ};
param L{I};
param LB{1..M} default -1e20;
param UB{1..M} default +1e20;
set EQ  := {m in 1..M:abs(LB[m]-UB[m])<1e-10};
set LEQ := {m in 1..M:LB[m]<=-1e20 && UB[m]< +1e20};
set GEQ := {m in 1..M:LB[m]> -1e20 && UB[m]>=+1e20};
set RNG := {m in 1..M:abs(LB[m]-UB[m])>1e-10 && LB[m]> -1e20 && UB[m]<+1e20};
var x{0..N-1};
minimize obj:sum{(0,i,j) in IJ}(Q[0,i,j]*x[i]*x[j])+sum{(0,i) in I}(L[0,i]*x[i]);
subject to ctr_eq {m in  EQ}:sum{(m,i,j) in IJ}(Q[m,i,j]*x[i]*x[j])+sum{(m,i) in I}(L[m,i]*x[i])=LB[m];
subject to ctr_geq{m in GEQ}:sum{(m,i,j) in IJ}(Q[m,i,j]*x[i]*x[j])+sum{(m,i) in I}(L[m,i]*x[i])>=LB[m];
subject to ctr_leq{m in LEQ}:sum{(m,i,j) in IJ}(Q[m,i,j]*x[i]*x[j])+sum{(m,i) in I}(L[m,i]*x[i])<=UB[m];
subject to ctr_rng{m in RNG}:LB[m]<=sum{(m,i,j) in IJ}(Q[m,i,j]*x[i]*x[j])+sum{(m,i) in I}(L[m,i]*x[i])<=UB[m];