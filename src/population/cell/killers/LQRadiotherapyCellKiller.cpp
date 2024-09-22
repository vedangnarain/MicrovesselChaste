/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "LQRadiotherapyCellKiller.hpp"
#include "CancerCellMutationState.hpp"
#include "StalkCellMutationState.hpp"
#include "BaseUnits.hpp"
#include "Secomb04Parameters.hpp"
#include "GenericParameters.hpp"
#include "Owen11Parameters.hpp"
// #include "QuiescentCancerCellMutationState.hpp"
#include <iostream>  // for debugging only

template<unsigned DIM>
LQRadiotherapyCellKiller<DIM>::LQRadiotherapyCellKiller(AbstractCellPopulation<DIM>* pCellPopulation) :
        AbstractCellKiller<DIM>(pCellPopulation),
        cancerousLinearRadiosensitivity(0.3 * unit::per_gray),
        cancerousQuadraticRadiosensitivity(0.03 * unit::per_gray_squared),
        normalLinearRadiosensitivity(0.15 * unit::per_gray),
        normalQuadraticRadiosensitivity(0.05 * unit::per_gray_squared),
        mDose(2.0 * unit::gray),
        mRadiationTimes(),
        mOerAlphaMax(1.75),
        mOerAlphaMin(1.0),
        mOerBetaMax(3.25),
        mOerBetaMin(1.0),
        mKOer(0.0*unit::mole_per_metre_cubed),
        mAlphaMax(0.3 * unit::per_gray),
        mBetaMax(0.03 * unit::per_gray_squared),
        mUseOer(false),
        mUseLewinModel(false),
        mUseScottModel(false),
        mUseWoutersModel(false)
{
    mKOer = 0.004499758549541243*unit::mole_per_metre_cubed;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::AddTimeOfRadiation(QTime time)
{
    mRadiationTimes.push_back(time);
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetDoseInjected(QAbsorbedDose d)
{
    mDose = d;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetTimeOfRadiation(std::vector<QTime > t)
{
    mRadiationTimes = t;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetOerAlphaMax(double value)
{
    mOerAlphaMax = value;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetOerAlphaMin(double value)
{
    mOerAlphaMin = value;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetOerBetaMax(double value)
{
    mOerBetaMax = value;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetOerBetaMin(double value)
{
    mOerBetaMin = value;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetOerConstant(QConcentration value)
{
    mKOer = value;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetAlphaMax(QPerAbsorbedDose value)
{
    mAlphaMax = value;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetBetaMax(QPerAbsorbedDoseSquared value)
{
    mBetaMax = value;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::UseOer(bool useOer)
{
    mUseOer = useOer;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::UseLewinModel(bool useLewinModel)
{
    mUseLewinModel = useLewinModel;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::UseScottModel(bool useScottModel)
{
    mUseScottModel = useScottModel;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::UseWoutersModel(bool useWoutersModel)
{
    mUseWoutersModel = useWoutersModel;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetCancerousRadiosensitivity(QPerAbsorbedDose alpha,
                                                                 QPerAbsorbedDoseSquared beta)
{
    cancerousLinearRadiosensitivity = alpha;
    cancerousQuadraticRadiosensitivity = beta;
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::SetNormalRadiosensitivity(QPerAbsorbedDose alpha,
                                                              QPerAbsorbedDoseSquared beta)
{
    normalLinearRadiosensitivity = alpha;
    normalQuadraticRadiosensitivity = beta;
}


template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
    if (!pCell->GetMutationState()->IsType<StalkCellMutationState>())
    {
        for (unsigned idx = 0; idx < mRadiationTimes.size(); idx++)
        {
            if (SimulationTime::Instance()->GetTime()*BaseUnits::Instance()->GetReferenceTimeScale() == mRadiationTimes[idx])
            {
                // Model radiation hit
                double death_probability = 0.0;
                if (mUseOer)
                {
                    // Get the cell's oxygen concentration
                    QConcentration oxygen = pCell->GetCellData()->GetItem("oxygen")*BaseUnits::Instance()->GetReferenceConcentrationScale();  // mol/metre_cubed
                    
                    // If we want to use a constant OER for hypoxic cells
                    if (mUseLewinModel)
                    {
                        // Set the oxygen solubility at STP and the radioresistant threshold
                        QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("LQRadiotherapyCellKiller") * GenericParameters::mpGasConcentrationAtStp->GetValue("LQRadiotherapyCellKiller");  // 1/Pa * mol/metre_cubed 
                        QPressure radioresistant_threshold = Owen11Parameters::mpOxygenPartialPressureAtHalfMaxCycleRateNormal->GetValue("LQRadiotherapyCellKiller");  // Pa
                        double oer_constant; 

                        // For debugging
                        // std::cout << "oxygen (mol/m3): " << oxygen << std::endl;
                        // std::cout << "mpOxygenVolumetricSolubility (1/Pa)" << Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("LQRadiotherapyCellKiller") << std::endl;
                        // std::cout << "mpGasConcentrationAtStp (mol/metre_cubed_mmHg): " << GenericParameters::mpGasConcentrationAtStp->GetValue("LQRadiotherapyCellKiller") << std::endl;
                        // std::cout << "oxygen (Pa): " << oxygen/oxygen_solubility_at_stp << std::endl;
                        // std::cout << "radioresistant_threshold (Pa): " << radioresistant_threshold << std::endl;

                        // If the cell is a hypoxic cancer cell, use the OER constant
                        if (oxygen/oxygen_solubility_at_stp<=radioresistant_threshold)
                        {
                            oer_constant = 3.0;
                        }
                        else
                        {
                            oer_constant = 1.0;
                        }
                        QPerAbsorbedDose alpha = mAlphaMax / oer_constant;
                        QPerAbsorbedDoseSquared beta = mBetaMax / (oer_constant * oer_constant);
                        death_probability = 1.0 - exp(-alpha * mDose - beta * mDose * mDose);
                    }
                    else if (mUseScottModel)
                    {
                        // Scott model
                        double oer_alpha = (mOerAlphaMax - mOerAlphaMin) * mKOer / (oxygen + mKOer) + mOerAlphaMin;
                        double oer_beta = (mOerBetaMax - mOerBetaMin) * mKOer / (oxygen + mKOer) + mOerBetaMin;
                        QPerAbsorbedDose alpha = mAlphaMax / oer_alpha;
                        QPerAbsorbedDoseSquared beta = mBetaMax / (oer_beta * oer_beta);
                        death_probability = 1.0 - exp(-alpha * mDose - beta * mDose * mDose);
                    }
                    else if (mUseWoutersModel)
                    {
                        // Wouters and Brown (1997)
                        double oer_alpha = (oxygen*mOerAlphaMax + mKOer)/(oxygen + mKOer);
                        double oer_beta = (oxygen*mOerBetaMax + mKOer)/(oxygen + mKOer);                        
                        double alpha_h = mAlphaMax/ mOerAlphaMax;
                        double beta_h = mBetaMax/ (mOerBetaMax * mOerBetaMax);
                        QPerAbsorbedDose alpha = alpha_h * oer_alpha;
                        QPerAbsorbedDoseSquared beta = beta_h * oer_beta * oer_beta;
                        death_probability = 1.0 - exp(-alpha * mDose - beta * mDose * mDose);
                    }
                }
                else
                {
                    death_probability = 1.0 - exp(-cancerousLinearRadiosensitivity * mDose
                                    - cancerousQuadraticRadiosensitivity * mDose * mDose);
                }
                if (!pCell->HasApoptosisBegun() && RandomNumberGenerator::Instance()->ranf() < death_probability)
                {
                    pCell->StartApoptosis();
                }
            }
        }
    }
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
            cell_iter != this->mpCellPopulation->End(); ++cell_iter)
    {
        CheckAndLabelSingleCellForApoptosis(*cell_iter);
    }
}

template<unsigned DIM>
void LQRadiotherapyCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile;

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class LQRadiotherapyCellKiller<1> ;
template class LQRadiotherapyCellKiller<2> ;
template class LQRadiotherapyCellKiller<3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LQRadiotherapyCellKiller)
