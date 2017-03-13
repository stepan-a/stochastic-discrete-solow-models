/* © 2017 Stéphane Adjemian <stephane.adjemian@univ-lemans.fr>
 *
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with the file.  If not, see <http://www.gnu.org/licenses/>.
 */

var Efficiency                  // $A$
    EfficiencyGrowth            // $X$
    Population                  // $L$
    PopulationGrowth            // $N$
    Output                      // $Y$
    PhysicalCapitalStock ;      // $K$

varexo e_x   // $\varepsilon_x$
       e_n ; // $\varepsilon_n$

parameters alpha                               // $\alpha$
	   delta                               // $\delta$
	   s                                   // $s$
           rho_x                               // $\rho_x$
           rho_n                               // $\rho_n$
           EfficiencyGrowth_ss                 // $X^{\star}$
           PopulationGrowth_ss ;               // $N^{\star}$

alpha = .33;
delta = .02;
s     = .20;
rho_x = .90;
rho_n = .95;
EfficiencyGrowth_ss = 1.00; // Do not change this calibration
PopulationGrowth_ss = 1.00; // Do not change this calibration

alpha = .33;
delta = .02;
s     = .20;
rho_x = .90;
rho_n = .95;
EfficiencyGrowth_ss = 1.02;
PopulationGrowth_ss = 1.02;

model;
    Efficiency = EfficiencyGrowth*Efficiency(-1);
    EfficiencyGrowth/EfficiencyGrowth_ss = (EfficiencyGrowth(-1)/EfficiencyGrowth_ss)^(rho_x)*exp(e_x);
    Population = PopulationGrowth*Population(-1);
    PopulationGrowth/PopulationGrowth_ss = (PopulationGrowth(-1)/PopulationGrowth_ss)^(rho_n)*exp(e_n);
    Output = PhysicalCapitalStock(-1)^alpha*(Efficiency*Population)^(1-alpha);
    PhysicalCapitalStock = (1-delta)*PhysicalCapitalStock(-1) + s*Output;
end;

histval;
    Efficiency(0) = 1;
    EfficiencyGrowth(0) = 1.02;
    Population(0) = 1;
    PopulationGrowth(0) = 1.02;
    PhysicalCapitalStock(0) = 1;
end;

shocks;
    var e_x = 0.005;
    var e_n = 0.001;
end;

oo_ = simul_backward_nonlinear_model([], 5000, options_, M_, oo_);

plot(log(oo_.endo_simul(5,:)));
title('log(Output)');
axis tight
box on
