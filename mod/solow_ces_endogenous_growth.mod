/* © Stéphane Adjemian 2017 <stephane.adjemian@univ-lemans.fr>
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
    Population                  // $L$
    Output                      // $Y$
    PhysicalCapitalStock ;      // $K$

varexo e_a   // $\varepsilon_a$
       e_l ; // $\varepsilon_l$

parameters alpha                               // $\alpha$
           epsilon                             // $\varepsilon$
	   delta                               // $\delta$
	   s                                   // $s$
           rho_a                               // $\rho_a$
           rho_l                               // $\rho_l$
           Efficiency_ss                       // $A^{\star}$
           Population_ss ;                     // $L^{\star}$

alpha   = 0.33;
epsilon = 2.00;
delta   = 0.02;
s       = 0.20;
rho_a   = 0.90;
rho_l   = 0.95;

Efficiency_ss = 1.0;
Population_ss = 1.0;

if s>delta*alpha^(-epsilon/(epsilon-1))
    disp('The model does not admit a steady state => Endogenous growth.')
end


model;
    Efficiency/Efficiency_ss = (Efficiency(-1)/Efficiency_ss)^rho_a*exp(e_a);
    Population/Population_ss = (Population(-1)/Population_ss)^rho_l*exp(e_l);
    Output = (alpha*PhysicalCapitalStock(-1)^((epsilon-1)/epsilon)+(1-alpha)*(Efficiency*Population)^((epsilon-1)/epsilon))^(epsilon/(epsilon-1));
    PhysicalCapitalStock = (1-delta)*PhysicalCapitalStock(-1) + s*Output;
end;

histval;
    Efficiency(0) = 1;
    Population(0) = 1;
    PhysicalCapitalStock(0) = .5;
end;

shocks;
    var e_a = 0.005;
    var e_l = 0.001;
end;

oo_ = simul_backward_nonlinear_model([], 5000, options_, M_, oo_);
