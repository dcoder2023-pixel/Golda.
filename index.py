import numpy as np
from typing import Union, Tuple

class AcidDissociationCalculator:
    """Comprehensive calculator for acid dissociation constants"""
    
    def __init__(self):
        pass
    
    def ka_from_ph(self, pH: float, C: float, approximation: bool = True) -> float:
        """
        Calculate Ka from pH and initial concentration
        
        Parameters:
        -----------
        pH : float
            pH of the solution
        C : float
            Initial concentration of acid (M)
        approximation : bool
            If True, uses [H⁺]²/C approximation for very weak acids
            If False, solves exact quadratic equation
            
        Returns:
        --------
        float : Ka value
        """
        [H] = 10 ** (-pH)
        
        if approximation and (H / C < 0.05):  # 5% rule
            # Approximation for very weak acids
            Ka = (H ** 2) / C
        else:
            # Exact solution
            Ka = (H ** 2) / (C - H)
        
        return Ka
    
    def ph_from_ka(self, Ka: float, C: float) -> float:
        """
        Calculate pH from Ka and concentration
        
        Parameters:
        -----------
        Ka : float
            Acid dissociation constant
        C : float
            Initial concentration of acid (M)
            
        Returns:
        --------
        float : pH of solution
        """
        # Solve quadratic: x² + Ka*x - Ka*C = 0
        # where x = [H⁺]
        discriminant = Ka**2 + 4*Ka*C
        H = (-Ka + np.sqrt(discriminant)) / 2
        
        pH = -np.log10(H)
        return pH
    
    def pka_from_ka(self, Ka: float) -> float:
        """Calculate pKa from Ka"""
        return -np.log10(Ka)
    
    def ka_from_pka(self, pKa: float) -> float:
        """Calculate Ka from pKa"""
        return 10 ** (-pKa)
    
    def degree_of_dissociation(self, Ka: float, C: float) -> float:
        """
        Calculate degree of dissociation (α)
        α = [H⁺]/C
        """
        pH = self.ph_from_ka(Ka, C)
        H = 10 ** (-pH)
        return H / C
    
    def buffer_ph(self, pKa: float, acid_conc: float, salt_conc: float) -> float:
        """
        Calculate pH of buffer using Henderson-Hasselbalch equation
        
        pH = pKa + log([A⁻]/[HA])
        """
        pH = pKa + np.log10(salt_conc / acid_conc)
        return pH

# Example usage
calc = AcidDissociationCalculator()

# Example 1: Calculate Ka from known pH
ka_value = calc.ka_from_ph(2.89, 0.100)
print(f"Ka from pH: {ka_value:.2e}")
print(f"pKa: {calc.pka_from_ka(ka_value):.2f}")

# Example 2: Calculate pH from known Ka
pH_value = calc.ph_from_ka(1.8e-5, 0.10)  # Acetic acid
print(f"pH from Ka: {pH_value:.2f}")

# Example 3: Buffer calculation
buffer_pH = calc.buffer_ph(4.76, 0.10, 0.10)  # Acetic acid/acetate
print(f"Buffer pH: {buffer_pH:.2f}")






