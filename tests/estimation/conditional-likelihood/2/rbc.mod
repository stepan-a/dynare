var k A c l i y B;

varexo e_a e_b;

parameters alp bet tet tau delt rho_a rho_b ;

alp = 0.4;
bet = 0.99;
tet = 0.357;
tau =  50 ;
delt = 0.02;
rho_a = 0.95;
rho_b = 0.80;

model;

c = ((1 - alp)*tet/(1-tet))*A*(1-l)*((k(-1)/l)^alp) ;
y = A*(k(-1)^alp)*(l^(1-alp)) ;
i = y-c ;
k = (1-delt)*k(-1) + i ;
log(A) = rho_a*log(A(-1)) + e_a ;
log(B) = rho_b*log(B(-1)) + e_b ;
(((c^(tet))*((1-l)^(1-tet)))^(1-tau))/c - bet*B(1)/B*((((c(+1)^(tet))*((1-l(+1))^(1-tet)))^(1-tau))/c(+1))*(1 -delt+alp*(A(1)*(k^alp)*(l(1)^(1-alp)))/k)=0 ;

end;

shocks;
var e_a; stderr 0.0035;
var e_b; stderr 0.0025;
end;

steady_state_model;
kss = -(alp-1)*(alp^(1/(1-alp)))*(bet^(1/(1-alp)))*((bet*(delt-1)+1)^(alp/(alp-1)))*tet;
k = kss/(-alp*delt*bet+delt*bet+alp*tet*bet-bet-alp*tet+1);
lss = (alp-1)*(bet*(delt-1)+1)*tet;
l = lss/(alp*tet+bet*((alp-1)*delt-alp*tet+1)-1) ;
y = (k^alp)*(l^(1-alp)) ;
i = delt*k ;
c = y - i ;
A = 1;
B = 1;
end;

%if ~exist('rbcdata.m', 'file')
steady;
stoch_simul(nograph, irf=0, periods=1000, order=2);
datatomfile('rbcdata');
%end;

estimated_params;
alp, uniform_pdf,,, 0.0001, 0.99;
bet, uniform_pdf,,, 0.0001, 0.99999;
tet, uniform_pdf,,, 0.0001, .999;
tau, uniform_pdf,,, 0.0001, 100;
delt, uniform_pdf,,, 0.0001, 0.05;
rho_a, uniform_pdf,,, 0.0001, 0.9999;
rho_b, uniform_pdf,,, 0.0001, 0.9999;
stderr e_a, uniform_pdf,,, 0.00001, 0.1;
stderr e_b, uniform_pdf,,, 0.00001, 0.1;
end;

estimated_params_init;
  alp, 0.4;
  bet, 0.99;
  tet, 0.357;
  tau, 50;
  delt, 0.02;
  rho_a, 0.95;
  rho_b, 0.80;
  stderr e_a, .035;
  stderr e_b, .025;
end;


varobs y c ;

data(file='./rbcdata.m');

estimation(nograph, conditional_likelihood, order=2, mode_compute=4, mh_replic=0,mode_check);
