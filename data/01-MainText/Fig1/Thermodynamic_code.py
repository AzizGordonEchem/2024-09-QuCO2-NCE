# -*- coding: utf-8 -*-

from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize, OptimizeResult
from scipy.special import expit

K1_values = np.logspace(4, 24, num=120)

for K1 in K1_values:

    @dataclass
    class PhysicalConstants:
        """Physical constants for aqueous CO2 equilibria with AQDS"""

        # Equlibrium constant for water dissociation in units of M^-2
        # H+ + OH- <-> H20
        # [H+] [OH-] = Kw
        Kw: float = pow(10.0, -14)

        # Henry's constant for CO2 in water, in units of M/bar
        # CO2(g) <-> CO2(aq)
        # [CO2(aq)] = KH * [pCO2(g)]
        KH: float = 0.035

        # Equilibrium constant for formation of C02 adduct, in units of M
        # AQ-- + 2 CO2(aq) <-> AQ(C02)2--
        # [AQ(C02)2--] = K1 * [AQ--] [CO2(aq)]^2
        K1: float = K1

        # Equilibrium constant for first protonation of AQDS
        # AQ-- + H+ <-> AQH-
        # [AQH-] = K2 * [AQ--] [H+]
        K2: float = pow(10.0, 13)

        # Equilibrium constant for second protonation of AQDS
        # AQH- + H+ <-> AQH2
        # [AQH2] = K3 * [AQH-] [H+]
        K3: float = pow(10.0, 11)

        # Equilibrium constant for protonation of bicarbonate into water and CO2
        # CO2 + H20 <-> HCO3- + H+
        # [HCO3-] [H+] = K4 * [CO2]
        K4: float = 1.1E-6

        # Equilibrium constant for deprotonation of bicarbonate to carbonate
        # HCO3- <-> CO3-- + H+ 
        # [CO3--] [H+] = K5 * [HCO3-]
        K5: float = 4.1E-10

    # Create a single default instance of PhysicalConstants
    pc_dflt: PhysicalConstants = PhysicalConstants()

    # *************************************************************************************************
    @dataclass
    class SpeciesConcentration:
        """State of the system: concentrations of all nine aqueous species in one dataclass"""
        # Concentration of H+ in units of M
        H: float
        # Concentration of OH- in units of M
        OH: float
        # Concentration of CO2(aq) in units of M
        CO2: float
        # Concentration of HCO3- (bicarbonate) in units of M
        HCO3: float
        # Concentration of CO3-- (carbonate) in units of M
        CO3: float
        # Concentration of AQ-- (deprotonated AQDS) in units of M
        AQ: float
        # Concentration of AQH- (protonated AQDS) in units of M
        AQH: float
        # Concentration of AQH2 (doubly protonated AQDS) in units of M
        AQH2: float
        # Concentration of AQ(CO2)2-- (CO2 adduct of AQDS) in units of M
        AQCO2_2: float
        # Amount of potassium ions from AQDS salt, in units of M
        K: float
        # Concerved quantity: total concentration of AQDS in units of M (AQ + AQH + AQH2 + AQCO2_2) = [K]
        AQ_total: float

        def aq_sum(self) -> float:
            """Sum of AQDS over four species; should equal the conserved total """
            return self.AQ + self.AQH + self.AQH2 + self.AQCO2_2
        
        def aq_inbal_abs(self) -> float:
            """Difference between simulated sum of AQDS species and conserved total, in absolute terms (M)"""
            return self.aq_sum() - self.AQ_total
        
        def aq_inbal_rel(self) -> float:
            """Difference between sum of AQDS and conserved total, in relative terms (dimensionless)"""
            return self.aq_inbal_abs() / self.AQ_total
        
        def charge_pos(self) -> float:
            """Total positive charge in the solution; in M of electrons equivalent"""
            # Review all nine species and their charges
            # H       = +1
            # OH      = -1
            # CO2     =  0
            # HC03    = -1
            # CO3     = -2
            # AQ      = -2
            # AQH     = -1
            # AQH2    =  0
            # AQCO2_2 = -2
            # K       = +1
            return (self.H + self.K)

        def charge_neg(self) -> float:
            """Total negative charge in the solution; in M of electrons equivalent"""
            return self.OH + self.HCO3 + self.AQH + 2.0 * (self.AQ + self.CO3 + self.AQCO2_2)

        def charge_net(self) -> float:
            """Total net charge in the solution; in M of electrons equivalent"""
            return self.charge_pos() - self.charge_neg()

        def charge_inbal_abs(self) -> float:
            """Charge inbalance in the solution; in M of electrons equivalent"""
            return abs(self.charge_net())
        
        def charge_inbal_rel(self) -> float:
            """Charge inbalance in the solution; in relative terms (dimensionless)"""
            return (2.0 * self.charge_inbal_abs()) / (self.charge_pos() + self.charge_neg())
        
        def error(self) -> float:
            """A single scalar function summarizing the error from the two inbalance functions"""
            # The relative AQ inbalance
            err_aq: float = self.aq_inbal_rel()
            # The relative charge inbalance
            err_charge: float = self.charge_inbal_rel()
            # Use the square root of the unweighted sum of squares of these two errors (i.e. the Euclidean norm)
            return np.hypot(err_aq, err_charge)
        
        def __str__(self) -> str:
            """One string writing out the nine concentrations"""
            out: str = "Species Concentrations:\n"
            out += '--------------------------------\n'
            out += f'     pH: {-np.log10(self.H):6.3f}\n'
            out += f'      H: {self.H:10.3e}\n'
            out += f'     OH: {self.OH:10.3e}\n'
            out += f'    CO2: {self.CO2:10.3e}\n'
            out += f'   HCO3: {self.HCO3:10.3e}\n'
            out += f'    CO3: {self.CO3:10.3e}\n'
            out += f'     AQ: {self.AQ:10.3e}\n'
            out += f'    AQH: {self.AQH:10.3e}\n'
            out += f'   AQH2: {self.AQH2:10.3e}\n'
            out += f'AQCO2_2: {self.AQCO2_2:10.3e}\n'
            out += f'      K: {self.K:10.3e}\n'
            out += f'\nTotal AQ and Charge:\n'
            out += '--------------------------------\n'
            out += f' AQ_sum: {self.aq_sum():10.6e}\n'
            out += f' charge: {self.charge_net():10.6e}\n'
            out += f'\nRelative Errors:\n'
            out += '--------------------------------\n'
            out += f' aq_sum: {self.aq_inbal_rel():10.3e}\n'
            out += f' charge: {self.charge_inbal_rel():10.3e}\n'
            out += f' error:  {self.error():10.3e}\n'
            return out

    # *************************************************************************************************
    def species_conc(H: float, AQH2: float, AQ_total: float, pCO2: float, pc: PhysicalConstants = pc_dflt) \
            -> SpeciesConcentration:
        """
        Populate the concentration of all nine species given the minimal set of inputs
        INPUTS:
            H:              The concentration of protons in M
            AQH2:           The concentration of AQH2 in M
            AQ_tota;:       The total concentration of all four AQ species in M; this is a conserved quantity
            pCO2:           The partial pressure of CO2 in bar
            pc:             The PhysicalConstants dataclass
        OUTPUTS:
            sc:             The SpeciesConcentration dataclass
        """ 
        # Compute hydroxide concentration from the water equilibrium
        # [H+] [OH-] = Kw
        OH: float = pc.Kw / H
        
        # Compute the concentration of CO2(aq) from the Henry's law equilibrium
        # [CO2(aq)] = KH * [CO2(g)]
        CO2: float = pc.KH * pCO2

        # Compute the concentration of bicarbonate from the water / CO2 equilibrium (deprotonation of bicarbonate)
        # [HCO3-] [H+] = K4 * [CO2]
        HCO3: float = (pc.K4 * CO2) / H

        # Compute the concentration of carbonate from the equilibrium for protonation of carbonate
        # [CO3--] [H+] = K5 * [HCO3-]
        CO3: float = (pc.K5 * HCO3) / H

        # Compute the concentration of AQH- from the equilibrium for the second protonation of AQ--
        # [AQH2] = K3 * [AQH-] [H+]
        AQH: float = AQH2 / (pc.K3 * H)

        # Compute the concentration of AQ from the equilibrium for the first protonation of AQH-
        # [AQH-] = K2 * [AQ--] [H+]
        AQ: float = AQH / (pc.K2 * H)

        # Compute the concentration of AQ(CO2)2-- from the equilibrium for the formation of the CO2 adduct
        # [AQ(C02)2--] = K1 * [AQ--] [CO2(aq)]^2
        AQCO2_2: float = pc.K1 * AQ * CO2 * CO2

        # The amount of potassium matches the total of AQ species *2
        K: float = 2.0 * AQ_total

        # Return the SpeciesConcentration dataclass
        return SpeciesConcentration(H=H, OH=OH, CO2=CO2, HCO3=HCO3, CO3=CO3, 
                                    AQ=AQ, AQH=AQH, AQH2=AQH2, AQCO2_2=AQCO2_2, K=K, AQ_total=AQ_total)

    # *************************************************************************************************
    def solve_CO2_conc(AQ_total: float, pCO2: float):
        """
        Solve for the species of all nine species given experimental conditions
        INPUTS:
            AQ_total:   The total concentration of all four AQ species in M; this is a conserved quantity
            pCO2:       Partial pressure of CO2 in bar
        """

        # Create an objective function for the optimizer
        def fun(x0: np.ndarray) -> float:
            """x0 is a vector of two values: log[H] and log[AQ]"""
            # Unpack x0
            x: float = x0[0]
            y: float = x0[1]
            # H is exponential of first control variable
            H: float = np.exp(x)
            # Relative abundance of AQH2 as fraction of total is sigmoid function of y
            AQH2: float = expit(y) * AQ_total
            # Compute the concentration of all nine species from H and AQH2
            sc: SpeciesConcentration = species_conc(H=H, AQH2=AQH2, AQ_total=AQ_total, pCO2=pCO2)
            # Return the error in these concentrations
            return sc.error()

        # Initial guess for the proton concentration
        H: float = 1.0E-7
        x: float = np.log(H)
        # Initial guess for the AQH2 fraction of all AQ species is 0.5 with logit 0
        y: float = 0.0
        # Wrap initial guess into one array
        x0: np.ndarray = np.array([x, y])

        # Set optimizer options
        tol: float = np.finfo(dtype=np.float64).eps
        options: dict = {'maxiter': 60000, 'disp': False}
        # Call the optimizer
        res: OptimizeResult = minimize(fun=fun, x0=x0, method='Nelder-Mead', tol=tol*50, options=options)

        # Unpack the result
        H: float = np.exp(res.x[0])
        AQH2: float = expit(res.x[1]) * AQ_total
        # Compute the concentration of all nine species from H and AQ
        sc: SpeciesConcentration = species_conc(H=H, AQH2=AQH2, AQ_total=AQ_total, pCO2=pCO2)

        # Return the SpeciesConcentration dataclass and the OptimizeResult
        return sc, res

    # *************************************************************************************************
    def main():
        """Entry for console program"""

        # Desired total concentration of AQDS in M
        AQ_total: float = 0.1

        # Desired partial pressure of CO2 in bar
        pCO2: float = 0.1

        # Solve for the species concentrations
        sc: SpeciesConcentration
        res: OptimizeResult
        sc, res = solve_CO2_conc(AQ_total=AQ_total, pCO2=pCO2)

        # Print the results
        print(sc)
        print("\nOptimization Results:")
        print('--------------------------------')
        print(res)
        
        filename = f"concentrationofcarbontt_{K1:.6e}.txt"
        with open(filename, 'w') as f:
            f.write("K1\CO2\HCO3\CO3\AQCO2_2\H\n")
            f.write("{0:.64f} " " {1:.64f} " " {2:.64f}" " {3:.64f}" " {4:.64f}" " {5:.64f} \n".format(K1,sc.CO2,sc.HCO3,sc.CO3,sc.AQCO2_2,sc.H))

        
  

    # *************************************************************************************************
    if __name__ == '__main__':
        main()
