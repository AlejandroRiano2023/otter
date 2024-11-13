#pragma once

// Including the "ADKernel" base class here so we can extend it
#include "AuxKernel.h"

class CLCalc2 : public AuxKernel
{
public:
    static InputParameters validParams();
    CLCalc2(const InputParameters & parameters);

protected:
  /// ADKernel objects must override precomputeQpResidual
  virtual RealVectorValue computeValue() override;

private:
    const VariableValue & _C_L;                 // Concentración C_L
    const VariableGradient & _grad_C_L;         // Gradient C_L
    const VariableValue & _der_C_L;             // Derivate of C_L in time
    const VariableValue & _d;                   // Variable d
    const VariableValue & _hidrostatic_stress;  // Tensión hidrostática
    const VariableGradient & _grad_hidrostatic_stress;
    const VariableValue & _N_T;                 // Variable N_T
    const VariableValue & _effective_plastic_strain; // Deformación plástica efectiva
    const VariableValue & _der_effective_plastic_strain;
    const Real _D_L;                           // Coeficiente de difusión D_L (micrometros cuadrados por segundo)
    const Real _V_H;                           // Volumen o velocidad V_H (micrometros cúbicos por mol)
    const Real _R;                             // Constante de gas R
    const Real _Theta;                         // Temperatura Theta
    const Real _Beta_L;
    const Real _N_L;
    const Real _W_B;
};
