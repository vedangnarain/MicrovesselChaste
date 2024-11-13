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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vedang's Notes (vedang.narain@msdtc.ox.ac.uk)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

This script contains two tests to generate an oxygen field from a source in the middle of a domain, with the consumption represented by (a) a uniform sink term (b) a sink term for individual cells (population uniformly distributed). It looks like the consumption rate, despite being set to the same value in simulations (a) and (b), produces different oxygen distributions. The lowest difference is obtained assuming a grid/cell size of ~1 um.

18/9/24

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialisation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TESTOXYGENCONSUMPTIONBUG_HPP_
#define TESTOXYGENCONSUMPTIONBUG_HPP_
#define _BACKWARD_BACKWARD_WARNING_H 1  // Cut out the VTK deprecated warning

// Essential functionality
#include <boost/lexical_cast.hpp>
#include <cxxtest/TestSuite.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>  // needed for exp function
#include <sstream>
#include <string>
// #include <time.h>
// #include <vector>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"

// Geometry tools
// #include "Part.hpp"

// Dimensional analysis
#include "BaseUnits.hpp"
#include "UnitCollection.hpp"
#include "GenericParameters.hpp"
#include "ParameterCollection.hpp"
#include "Owen11Parameters.hpp"
#include "Secomb04Parameters.hpp"
// #include "Vertex.hpp"

// Grids and PDEs
#include "CellBasedDiscreteSource.hpp"
#include "CellStateDependentDiscreteSource.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
// #include "SimpleLinearEllipticFiniteElementSolver.hpp"
// #include "DiscreteContinuumBoundaryCondition.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "RegularGrid.hpp"
#include "GridCalculator.hpp"
// #include "CellwiseSourceEllipticPde.hpp"
// #include "ConstBoundaryCondition.hpp"
// #include "EllipticGrowingDomainPdeModifier.hpp"
#include "VesselBasedDiscreteSource.hpp"

// Vessel networks
#include "Vessel.hpp"
#include "VesselNetwork.hpp"
// #include "VesselNetworkReader.hpp"
#include "VesselNetworkGenerator.hpp"
#include "VesselNetworkGeometryCalculator.hpp"
#include "VesselNetworkPropertyManager.hpp"
#include "VesselNode.hpp"
#include "VesselSegment.hpp"

// Flow
#include "NodeFlowProperties.hpp"
#include "VesselImpedanceCalculator.hpp"
#include "FlowSolver.hpp"
// #include "WallShearStressCalculator.hpp"
// #include "MechanicalStimulusCalculator.hpp"
// #include "MetabolicStimulusCalculator.hpp"
// #include "ShrinkingStimulusCalculator.hpp"
// #include "StructuralAdaptationSolver.hpp"
#include "ViscosityCalculator.hpp"

// Haematocrit
// #include "AlarconHaematocritSolver.hpp"
#include "BetteridgeHaematocritSolver.hpp"
#include "ConstantHaematocritSolver.hpp"
// #include "LinnengarHaematocritSolver.hpp"
// #include "YangHaematocritSolver.hpp"
#include "PriesHaematocritSolver.hpp"
#include "PriesWithMemoryHaematocritSolver.hpp"

// Cells
#include "CellsGenerator.hpp"
// #include "CellBasedSimulationArchiver.hpp"
#include "CellLabel.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
// #include "UniformCellCycleModel.hpp"
// #include "HoneycombMeshGenerator.hpp"
// #include "CellProliferativeTypesCountWriter.hpp"
// #include "CellPopulationAreaWriter.hpp"
// #include "CellPopulationElementWriter.hpp"
// #include "RandomCellKiller.hpp"
#include "LQRadiotherapyCellKiller.hpp"
#include "CancerCellMutationState.hpp"
// #include "StalkCellMutationState.hpp"
#include "QuiescentCancerCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Owen11CellPopulationGenerator.hpp"
#include "Owen2011TrackingModifier.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ApoptoticCellKiller.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
// #include "TransitCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"

// Forces
#include "GeneralisedLinearSpringForce.hpp"

// Angiogenesis
#include "Owen2011SproutingRule.hpp"
#include "Owen2011MigrationRule.hpp"
#include "AngiogenesisSolver.hpp"

// Vessel regression solver
// #include "WallShearStressBasedRegressionSolver.hpp"

// General solver to collect all the flows
#include "MicrovesselSolver.hpp"
#include "MicrovesselSimulationModifier.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"

// Visualisation
// #include "MicrovesselVtkScene.hpp"
// #include "VtkSceneMicrovesselModifier.hpp"

// Keep this last
#include "PetscAndVtkSetupAndFinalize.hpp"

// Includes standard library
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tests
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set key vessel parameters
QLength domain_span(1000.0*unit::microns);

// Set up the domain parameters
QLength domain_x = domain_span;  // this should extend to the x-position of outlet node
QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
QLength mid_domain_y = domain_y*0.5;

// Set nodes based on an equilateral network
std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y);
std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(domain_x, mid_domain_y);

// Make segments 
std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);

// Make vessels
std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);

// Set grid spacing
QLength grid_spacing = 20.0_um;

// Set haematocrit
double initial_haematocrit = 0.3;

class TestHaematocritSolvers : public CxxTest::TestSuite
{

public:

  // Make a single line across a PDE grid that acts as a Dirichlet BC in 2D 
  void xTestSingleLineSource2D()
  {    
      // Add the vessels to a vessel network
      std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
      p_network->AddVessel(p_vessel_1);
      
      // No pruning
      std::ostringstream strs;
      strs << std::fixed << std::setprecision( 1 );
      strs << "TestDebuggingNetworks/SingleLineSource/";
      std::string str_directory_name = strs.str();
      auto p_output_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);         
      // p_network->RemoveVessel(p_vessel_7,true);

      // Set up the reference length for the simulation
      QLength reference_length(1_um);
      BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

      // Set up the grid for the finite difference solver
      auto p_grid = RegularGrid<2>::Create();
      p_grid->SetSpacing(grid_spacing);
      c_vector<unsigned, 3> dimensions;
      dimensions[0] = unsigned((domain_x)/(grid_spacing)+1); // num x
      dimensions[1] = unsigned((domain_y)/(grid_spacing)+1); // num_y
      dimensions[2] = 1;
      p_grid->SetDimensions(dimensions);

      // Choose the PDE
      std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();

      // Set the diffusivity and decay terms
      p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
      p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

      // Set the O2 concentration value
      QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
      QConcentration vessel_oxygen_concentration = (initial_haematocrit/0.45)*oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");
      // std::cout << Owen11Parameters::mpOxygenDiffusivity->GetValue("User") << std::endl; 
      // std::cout << -1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User") << std::endl; 
      // std::cout << vessel_oxygen_concentration << std::endl; 

      // Set the boundary condition to be the network acting as a Dirichlet source
      std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
      p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
      p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
      p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

      // Set up the finite difference solver for oxygen (which handles everything)
      SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
      solver.SetGrid(p_grid);
      solver.SetPde(p_oxygen_pde);
      solver.AddBoundaryCondition(p_vessel_boundary_condition);
      solver.SetVesselNetwork(p_network);
      solver.SetLabel("oxygen");
      solver.SetFileName("oxygen_solution_0");

      // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
      solver.SetFileHandler(p_output_file_handler);
      solver.SetWriteSolution(true);
      solver.Solve();
      p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "line_network.vtp");

      // Print the median oxygenation
      std::vector<double> solution = solver.GetSolution();
      std::sort(solution.begin(), solution.end());

      double median_oxygen;
      size_t size = solution.size();
      if (size % 2 == 0) {
          // Even number of elements
          median_oxygen = (solution[size / 2 - 1] + solution[size / 2]) / 2.0;
      } else {
          // Odd number of elements
          median_oxygen = solution[size / 2];
      }
      std::cout << "Median oxygen: " << median_oxygen << std::endl;

      // Dump our parameter collection to an xml file and, importantly, clear it for the next test
      ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
      ParameterCollection::Instance()->Destroy();
      BaseUnits::Instance()->Destroy();
  }

  // Make a single line across a PDE grid with flow and h-splitting in 2D 
  void TestSingleFlowSource2D()
  {    
      // Add the vessels to a vessel network
      std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
      p_network->AddVessel(p_vessel_1);
    
      unsigned h_solver=3;

      // Specify which nodes are the inlets and outlets
      p_network->GetNode(0)->GetFlowProperties()->SetIsInputNode(true);
      p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
      p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetIsOutputNode(true);
      p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));
      
      // Set segment radii values
      QLength vessel_radius(7.5_um);
      VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);
      
      // Set segment viscosity values
      QDynamicViscosity viscosity = 1.96*1.e-3*unit::poiseuille;
      auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
      p_viscosity_calculator->SetPlasmaViscosity(viscosity);
      p_viscosity_calculator->SetVesselNetwork(p_network);
      p_viscosity_calculator->Calculate();
      // VesselNetworkPropertyManager<2>::SetSegmentViscosity(p_network, viscosity);

      // Set up the impedance calculator
      VesselImpedanceCalculator<2> impedance_calculator = VesselImpedanceCalculator<2>();
      impedance_calculator.SetVesselNetwork(p_network);
      impedance_calculator.Calculate();

      // Set up the flow solver
      FlowSolver<2> flow_solver = FlowSolver<2>();
      flow_solver.SetVesselNetwork(p_network);
      flow_solver.SetUp();
      flow_solver.Solve();
      std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
      if (h_solver==1)
      {
          std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
          auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
          p_haematocrit_solver->SetVesselNetwork(p_network);
          p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
          p_abstract_haematocrit_solver = p_haematocrit_solver;       
      }
      else if (h_solver==2)
      {
          std::cout << "Now using PriesHaematocritSolver..." << std::endl;
          auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
          // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
          p_haematocrit_solver->SetVesselNetwork(p_network);
          p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
          p_abstract_haematocrit_solver = p_haematocrit_solver;                    
      }
      else if (h_solver==3)
      {
          std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
          auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
          // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
          p_haematocrit_solver->SetVesselNetwork(p_network);
          p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
          p_abstract_haematocrit_solver = p_haematocrit_solver;      
      }
      else if (h_solver==4)
      {
          std::cout << "Now using FungHaematocritSolver..." << std::endl;
          auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
          // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
          p_haematocrit_solver->SetVesselNetwork(p_network);
          p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
          p_abstract_haematocrit_solver = p_haematocrit_solver;      
      }

      // No pruning
      std::ostringstream strs;
      strs << std::fixed << std::setprecision( 1 );
      strs << "TestDebuggingNetworks/SingleFlowSource/";
      std::string str_directory_name = strs.str();
      auto p_output_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);         
      // p_network->RemoveVessel(p_vessel_7,true);

      // Set up the reference length for the simulation
      QLength reference_length(1_um);
      BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

      // Set up the grid for the finite difference solver
      auto p_grid = RegularGrid<2>::Create();
      p_grid->SetSpacing(grid_spacing);
      c_vector<unsigned, 3> dimensions;
      dimensions[0] = unsigned((domain_x)/(grid_spacing)+1); // num x
      dimensions[1] = unsigned((domain_y)/(grid_spacing)+1); // num_y
      dimensions[2] = 1;
      p_grid->SetDimensions(dimensions);

      // Choose the PDE
      std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();

      // Set the diffusivity and decay terms
      p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
      p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

      // Set up the discrete source
      std::shared_ptr<VesselBasedDiscreteSource<2> > p_vessel_source = VesselBasedDiscreteSource<2>::Create();
      QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
          GenericParameters::mpGasConcentrationAtStp->GetValue("User");
      QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
          Owen11Parameters::mpReferencePartialPressure->GetValue("User");
      p_vessel_source->SetReferenceConcentration(vessel_oxygen_concentration);
      p_vessel_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
      p_vessel_source->SetVesselPermeability(1.0*Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
      p_oxygen_pde->AddDiscreteSource(p_vessel_source);

      // Set up the finite difference solver for oxygen (which handles everything)
      auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
      p_oxygen_solver->SetPde(p_oxygen_pde);
      p_oxygen_solver->SetLabel("oxygen");
      p_oxygen_solver->SetGrid(p_grid);

      // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
      // double initial_haematocrit = 0.3;
      unsigned max_iter = 1000; 
      double tolerance2 = 1.e-10;
      std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
      std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
      for(unsigned idx=0;idx<max_iter;idx++)
      {
          // Run the solvers
          impedance_calculator.Calculate();
          flow_solver.SetUp();
          flow_solver.Solve();
          p_abstract_haematocrit_solver->Calculate();
          p_viscosity_calculator->Calculate();

          // Get the residual
          double max_difference = 0.0;
          double h_for_max = 0.0;
          double prev_for_max = 0.0;
          for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
          {
              double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
              double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
              if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
              {
                  max_difference = difference;
                  h_for_max = current_haematocrit;
                  prev_for_max = previous_haematocrit[jdx];
              }
              previous_haematocrit[jdx] = current_haematocrit;
          }
          std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

          // Print the final or intermediary convergence results
          if(max_difference<=tolerance2)  
          {
              std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
              break;
          }
          else
          {
              if(idx%1==0)
              {
                  std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                  // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                  // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                  // p_network->Write(output_file);
              }
          }
      }
      
      // Run the simulation 
      SimulationTime::Instance()->SetStartTime(0.0);
      SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
      auto p_microvessel_solver = MicrovesselSolver<2>::Create();
      p_microvessel_solver->SetVesselNetwork(p_network);
      p_microvessel_solver->SetOutputFileHandler(p_output_file_handler);
      p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
      p_microvessel_solver->Run();

      // Print the average oxygenation
      std::vector<double> solution = p_oxygen_solver->GetSolution();
      double average_oxygen = 0.0;
      for(unsigned jdx=0;jdx<solution.size();jdx++)
      {
          average_oxygen += solution[jdx];
      }
      average_oxygen /= double(solution.size());
      std::cout << "Average oxygen: " << average_oxygen << std::endl;
      
      // Print the median oxygenation
      // std::vector<double> solution = solver.GetSolution();
      std::sort(solution.begin(), solution.end());

      double median_oxygen;
      size_t size = solution.size();
      if (size % 2 == 0) {
          // Even number of elements
          median_oxygen = (solution[size / 2 - 1] + solution[size / 2]) / 2.0;
      } else {
          // Odd number of elements
          median_oxygen = solution[size / 2];
      }
      std::cout << "Median oxygen: " << median_oxygen << std::endl;

      // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
      std::string output_file = p_output_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
      p_network->Write(output_file);

      // Dump our parameter collection to an xml file and, importantly, clear it for the next test
      ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
      ParameterCollection::Instance()->Destroy();
      BaseUnits::Instance()->Destroy();
      SimulationTime::Instance()->Destroy();
  }

 };

// Tests with cells
class TestPaper3 : public AbstractCellBasedWithTimingsTestSuite
{

public:

  // Simulate a 2D hexagonal network on a PDE grid with flow and H-splitting (single inlet and outlet) and cells
  void xTestSingleLineSourceWithCACells2D()
  {
      // Set file name
      std::ostringstream strs;
      // strs << std::fixed << std::setprecision( 2 );
      strs << "TestDebuggingNetworks/SingleLineSourceWithCACells/";        

      // Seed the random number generator
      RandomNumberGenerator::Instance()->Reseed(12345);

      // Set up the reference length for the simulation
      QLength reference_length(1_um);
      BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
      QTime reference_time = 3600.0*unit::seconds;
      BaseUnits::Instance()->SetReferenceTimeScale(reference_time);



      auto p_grid = RegularGrid<2>::Create();
      // QLength grid_spacineg = Owen11Parameters::mpLatticeSpacing->GetValue("User");
      p_grid->SetSpacing(grid_spacing);
      c_vector<unsigned, 3> dimensions;


      // dimensions[0] = 51; // num x
      // dimensions[1] = 51; // num_y
      // dimensions[2] = 1;
      // p_grid->SetDimensions(dimensions);



      // Set up the grid for the finite difference solver
      // auto p_grid = RegularGrid<2>::Create();
      // p_grid->SetSpacing(grid_spacing);
      // c_vector<unsigned, 3> dimensions;
      dimensions[0] = unsigned((domain_x)/(grid_spacing))+1; // num x
      dimensions[1] = unsigned((domain_y)/(grid_spacing))+1; // num_y
      dimensions[2] = 1;
      // std::cout << " 0= " << dimensions[0] << " 1= " << dimensions[1]  << " 2= " << dimensions[2] << std::endl;
      p_grid->SetDimensions(dimensions);





      // Add the vessels to a vessel network
      std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
      p_network->AddVessel(p_vessel_1);
      
      // Match this based on cell size
      // QLength grid_spacing = 10.0_um;  // the simulation time gets quite long if you reduce the resolution further
  
      // Set key vessel parameters
      // QLength homogeneous_vessel_radius = 10_um;
      // QLength domain_span(100.0*unit::microns);

      // Set up the domain parameters
      // QLength domain_x = domain_span;  // this should extend to the x-position of outlet node
      // QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
      // QLength mid_domain_y = domain_y*0.5;

      // Set nodes based on an equilateral network
      // std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(50.0_um, mid_domain_y);
      // std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(51.0_um, mid_domain_y);

      // // Make segments 
      // std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);

      // // Make vessels
      // std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);

      // // Add the vessels to a vessel network
      // std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
      // p_network->AddVessel(p_vessel_1);



      // // Run the simulation with different solvers of interest
      // for (unsigned h_solver=1; h_solver<=1; h_solver++)
      // {
          // std::ostringstream strs;
          std::string str_directory_name = strs.str();
          auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);            
          
          // Set radius
          // p_segment = p_network->GetVesselSegments()[0];
          // p_segment->SetRadius(homogeneous_vessel_radius);
          // VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);

          // Save just the network architecture without flow
          // p_network->Write(p_file_handler->GetOutputDirectoryFullPath() + "network_structure.vtp");

          // Save just the grid
          // p_grid->Write(p_file_handler);

          // Set up the cell populations
          std::shared_ptr<GridCalculator<2> > p_grid_calc = GridCalculator<2>::Create();
          p_grid_calc->SetGrid(p_grid);
        
          std::shared_ptr<Owen11CellPopulationGenerator<2> > p_cell_population_genenerator = Owen11CellPopulationGenerator<2>::Create();
          p_cell_population_genenerator->SetGridCalculator(p_grid_calc);
          p_cell_population_genenerator->SetVesselNetwork(p_network);
          QLength tumour_radius(50.0 * unit::microns);
          p_cell_population_genenerator->SetTumourRadius(tumour_radius);
          std::shared_ptr<CaBasedCellPopulation<2> > p_cell_population = p_cell_population_genenerator->Update();

          // Choose the PDE
          std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();

          // Set the diffusivity and decay terms
          p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
          // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
          auto p_cell_oxygen_sink = CellBasedDiscreteSource<2>::Create();  // Set the cell terms
          // p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(13.0);
          p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
          p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

          // Set the O2 concentration value
          QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
          QConcentration vessel_oxygen_concentration = (initial_haematocrit/0.45)*oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");
          // std::cout << Owen11Parameters::mpOxygenDiffusivity->GetValue("User") << std::endl; 
          // std::cout << -1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User") << std::endl; 
          // std::cout << vessel_oxygen_concentration << std::endl; 
          
          // Set the boundary condition to be the network acting as a Dirichlet source
          std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
          p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
          p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
          p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

          // Set up the discrete vessel source
          // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
          // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
          //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
          // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
          //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
          // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
          // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
          // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
          // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

          // Set up the finite difference solver for oxyge
          auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
          p_oxygen_solver->SetPde(p_oxygen_pde);
          p_oxygen_solver->SetLabel("oxygen");
          p_oxygen_solver->SetGrid(p_grid);
          p_oxygen_solver->AddBoundaryCondition(p_vessel_boundary_condition);
          p_oxygen_solver->SetVesselNetwork(p_network);
          // solver.SetFileName("oxygen_solution_0");
          // p_oxygen_solver->SetWriteSolution(true);

          // Set up the VEGF PDE as the O2 one
          auto p_vegf_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
          p_vegf_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpVegfDiffusivity->GetValue("User"));
          // p_vegf_pde->SetIsotropicDiffusionConstant(0.0);
          p_vegf_pde->SetContinuumLinearInUTerm(-Owen11Parameters::mpVegfDecayRate->GetValue("User"));
      
          // Set up a map for different release rates depending on cell type. Also include a threshold intracellular VEGF below which there is no release.
          auto p_normal_and_quiescent_cell_source = CellStateDependentDiscreteSource<2>::Create();  // initialise the source (which combines normal and quiescent cells)
          std::map<unsigned, QConcentrationFlowRate > normal_and_quiescent_cell_rates;  // map to store secretion rates
          std::map<unsigned, QConcentration > normal_and_quiescent_cell_rate_thresholds;  // map to store secretion thresholds
          MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
          MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
          normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for normal cells
          normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 0.27*unit::mole_per_metre_cubed;  // set the threshold for normal cells
          // normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = 0.00000000000000000001;  // set the rate for normal cells
          // normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 1000000_M;  // set the threshold for normal cells
          normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for quiescent cells
          normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 0_M;  // set the threshold for quiescent cells
          // normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = 0.00000000000000000001;  // set the rate for quiescent cells
          // normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 1000000_M;  // set the threshold for quiescent cells
          p_normal_and_quiescent_cell_source->SetStateRateMap(normal_and_quiescent_cell_rates);  // configure the map
          p_normal_and_quiescent_cell_source->SetLabelName("VEGF");  // set the label for the source
          p_normal_and_quiescent_cell_source->SetStateRateThresholdMap(normal_and_quiescent_cell_rate_thresholds);
          // p_vegf_pde->AddDiscreteSource(p_normal_and_quiescent_cell_source);

          // Add a vessel related VEGF sink
          auto p_vessel_vegf_sink = VesselBasedDiscreteSource<2>::Create();
          p_vessel_vegf_sink->SetReferenceConcentration(0.0_M);
          p_vessel_vegf_sink->SetVesselPermeability(Owen11Parameters::mpVesselVegfPermeability->GetValue("User"));
          // p_vessel_vegf_sink->SetVesselPermeability(0.0);
          // p_vegf_pde->AddDiscreteSource(p_vessel_vegf_sink);

          // Set up a finite difference solver as before for the VEGF
          auto p_vegf_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
          p_vegf_solver->SetPde(p_vegf_pde);
          p_vegf_solver->SetLabel("VEGF_Extracellular");  // set the label for the external field
          p_vegf_solver->SetGrid(p_grid);

          // // Set up the viscosity solver
          // auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
          // p_viscosity_calculator->SetPlasmaViscosity(viscosity);
          // p_viscosity_calculator->SetVesselNetwork(p_network);
          // p_viscosity_calculator->Calculate();
          
          // // Set up the impedance solver
          // auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
          // p_impedance_calculator->SetVesselNetwork(p_network);
          // p_impedance_calculator->Calculate();

          // Set up the flow solver 
          // FlowSolver<2> flow_solver;
          // flow_solver.SetVesselNetwork(p_network);
          // // flow_solver.SetUp();
          // flow_solver.SetUseDirectSolver(true);
          // // flow_solver.Solve();

          // Set up flag for broken solver
          // unsigned broken_solver = 0;

          // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
          // unsigned max_iter = 1000;  // 1000 
          // double tolerance2 = 1.e-10;
          // std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
          // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
          // for(unsigned idx=0;idx<max_iter;idx++)
          // {
          //     // Run the solvers
          //     p_viscosity_calculator->Calculate();
          //     p_impedance_calculator->Calculate();
          //     flow_solver.SetUp();
          //     flow_solver.Solve();
          //     p_abstract_haematocrit_solver->Calculate();
          //     p_viscosity_calculator->Calculate();

          //     // Get the residual
          //     double max_difference = 0.0;
          //     double h_for_max = 0.0;
          //     double prev_for_max = 0.0;
          //     for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
          //     {
          //         double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
          //         double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
          //         if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
          //         {
          //             max_difference = difference;
          //             h_for_max = current_haematocrit;
          //             prev_for_max = previous_haematocrit[jdx];
          //         }
          //         previous_haematocrit[jdx] = current_haematocrit;
          //     }
          //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

          //     // Print the final or intermediary convergence results
          //     if(max_difference<=tolerance2)  
          //     {
          //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
          //         broken_solver = 0;
          //         break;
          //     }
          //     else
          //     {
          //         if(idx%1==0)
          //         {
          //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
          //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
          //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
          //             p_network->Write(output_file);
          //         }
          //     }

          //     // If there is no convergence after all the iterations, print the error message.
          //     if(idx==max_iter-1)
          //     {
          //         std::cout << "Problem encountered in " << str_directory_name << std::endl;
          //         error_log << "\n Problem encountered in " << str_directory_name << std::endl; 
          //         broken_solver = 1;
          //         break;
          //     }
          // }

          // // If solver doesn't converge, move on to next one
          // if (broken_solver == 1)
          // {
          //     continue;
          // }

          // Set up the angiogenesis solver (runs on deactivated VEGF function but seems to be necessary for main solver)
          std::shared_ptr<AngiogenesisSolver<2> > p_angiogenesis_solver = AngiogenesisSolver<2>::Create();
          std::shared_ptr<Owen2011SproutingRule<2> > p_sprouting_rule = Owen2011SproutingRule<2>::Create();
          std::shared_ptr<Owen2011MigrationRule<2> > p_migration_rule = Owen2011MigrationRule<2>::Create();
          p_angiogenesis_solver->SetMigrationRule(p_migration_rule);
          p_angiogenesis_solver->SetSproutingRule(p_sprouting_rule);
          p_sprouting_rule->SetDiscreteContinuumSolver(p_vegf_solver);
          p_migration_rule->SetDiscreteContinuumSolver(p_vegf_solver);
          p_angiogenesis_solver->SetVesselGridCalculator(p_grid_calc);
          p_angiogenesis_solver->SetVesselNetwork(p_network);
          
          // Initialise the solver
          auto p_microvessel_solver = MicrovesselSolver<2>::Create();
          p_microvessel_solver->SetVesselNetwork(p_network);
          p_microvessel_solver->SetOutputFrequency(5_h);
          // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
          p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
          p_microvessel_solver->AddDiscreteContinuumSolver(p_vegf_solver);
          p_microvessel_solver->SetAngiogenesisSolver(p_angiogenesis_solver);

          // Specifies which extracellular fields to update based on PDE
          boost::shared_ptr<MicrovesselSimulationModifier<2> > p_microvessel_modifier =
          boost::shared_ptr<MicrovesselSimulationModifier<2> >(new MicrovesselSimulationModifier<2> ());
          p_microvessel_modifier->SetMicrovesselSolver(p_microvessel_solver);
          std::vector<std::string> update_labels;
          update_labels.push_back("oxygen");
          p_microvessel_modifier->SetCellDataUpdateLabels(update_labels);

          // The full simulation is run as a typical Cell Based Chaste simulation
          OnLatticeSimulation<2> simulator(*p_cell_population);
          simulator.AddSimulationModifier(p_microvessel_modifier);

          // Add a killer to remove apoptotic cells (leave it commented so RT is the only killer)
          // boost::shared_ptr<ApoptoticCellKiller<2> > p_apoptotic_cell_killer(new ApoptoticCellKiller<2>(p_cell_population.get())); (leave it commented so RT is the only killer)
          // simulator.AddCellKiller(p_apoptotic_cell_killer);

          // Add a LQ RT killer
          // boost::shared_ptr<LQRadiotherapyCellKiller<2> > p_rt_killer =
          // boost::shared_ptr<LQRadiotherapyCellKiller<2> >(new LQRadiotherapyCellKiller<2> (p_cell_population.get()));
          // p_rt_killer->SetAlphaMax(0.3*unit::per_gray);
          // p_rt_killer->SetBetaMax(0.03*unit::per_gray_squared);
          // p_rt_killer->SetDoseInjected(2.0*unit::gray);
          // p_rt_killer->SetCancerousRadiosensitivity(0.3 * unit::per_gray, 0.03 * unit::per_gray_squared);
          // p_rt_killer->SetNormalRadiosensitivity(0.15 * unit::per_gray, 0.05 * unit::per_gray_squared);
          // p_rt_killer->SetOerAlphaMax(1.75);
          // p_rt_killer->SetOerAlphaMin(1.0);
          // p_rt_killer->SetOerBetaMax(3.25);
          // p_rt_killer->SetOerBetaMin(1.0);
          // // p_rt_killer->SetOerConstant(oxygen_solubility_at_stp * 3.28_Pa);
          // p_rt_killer->SetOerConstant(0.004499758549541243*unit::mole_per_metre_cubed);  // K_OER (3.28 mmHg)
          // p_rt_killer->SetAlphaMax(0.3 * unit::per_gray);
          // p_rt_killer->SetBetaMax(0.03 * unit::per_gray_squared);
          // p_rt_killer->UseOer(true);

          // Set Up the dosage times
          // std::vector<QTime > rt_times;
          // rt_times.push_back(3600.0*24.0*unit::seconds);
          // rt_times.push_back(3600.0*48.0*unit::seconds);
          // rt_times.push_back(3600.0*72.0*unit::seconds);
          // p_rt_killer->SetTimeOfRadiation(rt_times);
          // p_rt_killer->AddTimeOfRadiation(3600.0*96.0*unit::seconds);
          // simulator.AddCellKiller(p_rt_killer);

          // If needed
          // boost::shared_ptr<RandomCellKiller<2> > p_cell_killer(new RandomCellKiller<2>(p_cell_population.get(), 0.01));
          // simulator.AddCellKiller(p_cell_killer);
          // for (typename AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population.get()->Begin(); cell_iter != p_cell_population.get()->End(); ++cell_iter)
          // {
          //     // CheckAndLabelSingleCellForApoptosis(*cell_iter);
          //     std::cout << "hi cell" << (*cell_iter)->GetCellData()->GetItem("oxygen") << std::endl;

          //     //  std::cout <<  BaseUnits::Instance()->GetReferenceConcentrationScale() << std::endl;


          //     //     double oxygen_concentration;

          //     // std::cout << "hi cell" << oxygen_concentration << std::endl;
          // }

          // Add another modifier for updating cell cycle quantities.
          boost::shared_ptr<Owen2011TrackingModifier<2> > p_owen11_tracking_modifier(new Owen2011TrackingModifier<2>);
          simulator.AddSimulationModifier(p_owen11_tracking_modifier);

          // Print the average oxygenation
          // std::vector<double> solution = p_oxygen_solver->GetSolution();
          // double average_oxygen = 0.0;
          // for(unsigned jdx=0;jdx<solution.size();jdx++)
          // {
          //     average_oxygen += solution[jdx];
          // }
          // average_oxygen /= double(solution.size());
          // std::cout << "Average oxygen: " << average_oxygen << std::endl;

          // Set up the simulation time and run it
          simulator.SetOutputDirectory(str_directory_name);
          // simulator.SetSamplingTimestepMultiple(10);  // get sample every x steps (multiple*dt = time between s)
          simulator.SetSamplingTimestepMultiple(1);
          // simulator.SetDt(1);  // simulation time step in hours
          simulator.SetDt(1.0);
          simulator.SetEndTime(1.0);  // end time in hours
          // simulator.SetEndTime(24);  // end time in hours
          // simulator.SetEndTime(2.0);
          simulator.Solve();

          // // Print the median oxygenation
          std::vector<double> solution = p_oxygen_solver->GetSolution();
          std::sort(solution.begin(), solution.end());

          double median_oxygen;
          size_t size = solution.size();
          if (size % 2 == 0) {
              // Even number of elements
              median_oxygen = (solution[size / 2 - 1] + solution[size / 2]) / 2.0;
          } else {
              // Odd number of elements
              median_oxygen = solution[size / 2];
          }
          std::cout << "Median oxygen: " << median_oxygen << std::endl;

          // Write the output network file (visualise with Paraview: set Filters->Alphabetical->Tube)
          std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
          p_network->Write(output_file);

          // Dump our parameter collection to an xml file and, importantly, clear it for the next test
          ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
          ParameterCollection::Instance()->Destroy();
          BaseUnits::Instance()->Destroy();
          SimulationTime::Instance()->Destroy();
      // }
      // Print the error log
      // std::string error_message = error_log.str();
      // std::cout << error_message << std::endl; 
  }

  // Simulate a 2D hexagonal network on a PDE grid with flow and H-splitting (single inlet and outlet) and cells
  void xTestSingleLineSourceWithNodeBasedCells2D()
  {
    // Set file name
    std::ostringstream strs;
    // strs << std::fixed << std::setprecision( 2 );
    strs << "TestDebuggingNetworks/SingleLineSourceWithNodeBasedCells/";        

    // Seed the random number generator
    RandomNumberGenerator::Instance()->Reseed(12345);

    // Set up the reference length for the simulation
    QLength reference_length(1_um);
    BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
    QTime reference_time = 3600.0*unit::seconds;
    BaseUnits::Instance()->SetReferenceTimeScale(reference_time);

    // Change the resolution for the PDE fields
    QLength new_grid_spacing = 10.0_um;

    // Set up the grid
    auto p_grid = RegularGrid<2>::Create();
    // QLength grid_spacineg = Owen11Parameters::mpLatticeSpacing->GetValue("User");
    p_grid->SetSpacing(new_grid_spacing);
    c_vector<unsigned, 3> dimensions;
    dimensions[0] = unsigned((domain_x)/(new_grid_spacing))+1; // num x
    dimensions[1] = unsigned((domain_y)/(new_grid_spacing))+1; // num_y
    dimensions[2] = 1;
    // std::cout << " 0= " << dimensions[0] << " 1= " << dimensions[1]  << " 2= " << dimensions[2] << std::endl;
    p_grid->SetDimensions(dimensions);

    // Add the vessels to a vessel network
    std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
    p_network->AddVessel(p_vessel_1);
  


    // std::ostringstream strs;
    std::string str_directory_name = strs.str();
    auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);            

    // Set up the cell populations
    std::shared_ptr<GridCalculator<2> > p_grid_calc = GridCalculator<2>::Create();
    p_grid_calc->SetGrid(p_grid);
  

    const int DIM = 2;
    // Set up the cell nodes
    std::vector<Node<2>*> nodes;
    unsigned nodeNum=0;
    unsigned numTumourNodes = 0;
    std::vector<unsigned> location_indices;

    // For staggered rows
    unsigned rownum = 1;
    double offset;
    double cell_size = 10.0;
    
    // Mark the stroma nodes
    for (double x=0; x<=100.0; x=x+cell_size)//x++)
    {
        offset = 0;
        // if (rownum % 2)
        // {
        //     offset = cell_size/2.0;
        // }
        for (double y=0+offset; y<=100.0; y=y+cell_size)//y++)
        {
            if (DIM == 2) 
            {
                nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
                nodeNum++;
                numTumourNodes++;
            }
        }
        rownum++;
    }
    // std::cout << "unsigned(domain_x) = " << unsigned(domain_x) << std::endl;

    // std::cout << "numTumourNodes = " << numTumourNodes << std::endl;


    // Pick the cell mutation and proliferation state
    std::vector<CellPtr> cells;
    MAKE_PTR(WildTypeCellMutationState, p_wild_state);
    MAKE_PTR(StemCellProliferativeType, p_stem_state);

    // Create the cancer and stroma labels
    unsigned cancerLabelColour = 9;
    unsigned stromaLabelColour = 1;
    MAKE_PTR_ARGS(CellLabel, p_stroma_label, (stromaLabelColour));
    MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));


    // Print the cell types and their colours
    std::cout << "WildTypeCellMutationState has Legacy Cell type " <<  p_wild_state->GetColour() << std::endl;
    std::cout << "StemCellProliferativeType has Legacy Cell type " <<  p_stem_state->GetColour() << std::endl;
    std::cout << "p_stroma_label has Legacy Cell type " <<  p_stroma_label->GetColour() << std::endl;
    std::cout << "p_cancer_label has Legacy Cell type " <<  p_cancer_label->GetColour() << std::endl;
    // std::cout << "Cancer state has Legacy Cell type " <<  p_cancer_state->GetColour() << std::endl;

    // Seed the tumour cells
    for (unsigned int i=0; i<numTumourNodes; i++)  // iterate over the number of tumour nodes plus one to add a cell and location index for the origin
    {
        SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
        p_model->SetDimension(2);
        CellPtr p_cell(new Cell(p_wild_state, p_model));
        p_cell->SetCellProliferativeType(p_stem_state);
        p_cell->AddCellProperty(p_cancer_label);
        double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                        (  p_model->GetStemCellG1Duration()
                          + p_model->GetSG2MDuration() );
        p_cell->SetBirthTime(birth_time);
        cells.push_back(p_cell);
        location_indices.push_back(i);
    }
    // std::cout << "location_indices.size = " << location_indices.size() << std::endl;

    // Set up a mesh for the cells 
    NodesOnlyMesh<2> mesh;                       
    mesh.ConstructNodesWithoutMesh(nodes, 0.5*cell_size);                       
    NodeBasedCellPopulation<2> cell_population(mesh, cells, location_indices);                       
    cell_population.SetAbsoluteMovementThreshold(DBL_MAX);  //Set big movement threshold

    // std::shared_ptr<Owen11CellPopulationGenerator<2> > p_cell_population_genenerator = Owen11CellPopulationGenerator<2>::Create();
    // p_cell_population_genenerator->SetGridCalculator(p_grid_calc);
    // p_cell_population_genenerator->SetVesselNetwork(p_network);
    // QLength tumour_radius(50.0 * unit::microns);
    // p_cell_population_genenerator->SetTumourRadius(tumour_radius);
    // std::shared_ptr<CaBasedCellPopulation<2> > p_cell_population = p_cell_population_genenerator->Update();

    // Choose the PDE
    std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();

    // Set the diffusivity and decay terms
    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
    auto p_cell_oxygen_sink = CellBasedDiscreteSource<2>::Create();  // Set the cell terms
    // p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(13.0);
    p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
    p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

    // Set the O2 concentration value
    QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
    QConcentration vessel_oxygen_concentration = (initial_haematocrit/0.45)*oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");
    // std::cout << Owen11Parameters::mpOxygenDiffusivity->GetValue("User") << std::endl; 
    // std::cout << -1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User") << std::endl; 
    // std::cout << vessel_oxygen_concentration << std::endl; 

    // Set the boundary condition to be the network acting as a Dirichlet source
    std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
    p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
    p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
    p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

    // Set up the finite difference solver for oxyge
    auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
    p_oxygen_solver->SetPde(p_oxygen_pde);
    p_oxygen_solver->SetLabel("oxygen");
    p_oxygen_solver->SetGrid(p_grid);
    p_oxygen_solver->AddBoundaryCondition(p_vessel_boundary_condition);
    p_oxygen_solver->SetVesselNetwork(p_network);
    // solver.SetFileName("oxygen_solution_0");
    // p_oxygen_solver->SetWriteSolution(true);

    // Set up the VEGF PDE as the O2 one
    auto p_vegf_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
    p_vegf_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpVegfDiffusivity->GetValue("User"));
    // p_vegf_pde->SetIsotropicDiffusionConstant(0.0);
    p_vegf_pde->SetContinuumLinearInUTerm(-Owen11Parameters::mpVegfDecayRate->GetValue("User"));

    // Set up a map for different release rates depending on cell type. Also include a threshold intracellular VEGF below which there is no release.
    auto p_normal_and_quiescent_cell_source = CellStateDependentDiscreteSource<2>::Create();  // initialise the source (which combines normal and quiescent cells)
    std::map<unsigned, QConcentrationFlowRate > normal_and_quiescent_cell_rates;  // map to store secretion rates
    std::map<unsigned, QConcentration > normal_and_quiescent_cell_rate_thresholds;  // map to store secretion thresholds
    MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
    MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
    normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for normal cells
    normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 0.27*unit::mole_per_metre_cubed;  // set the threshold for normal cells
    // normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = 0.00000000000000000001;  // set the rate for normal cells
    // normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 1000000_M;  // set the threshold for normal cells
    normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for quiescent cells
    normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 0_M;  // set the threshold for quiescent cells
    // normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = 0.00000000000000000001;  // set the rate for quiescent cells
    // normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 1000000_M;  // set the threshold for quiescent cells
    p_normal_and_quiescent_cell_source->SetStateRateMap(normal_and_quiescent_cell_rates);  // configure the map
    p_normal_and_quiescent_cell_source->SetLabelName("VEGF");  // set the label for the source
    p_normal_and_quiescent_cell_source->SetStateRateThresholdMap(normal_and_quiescent_cell_rate_thresholds);
    // p_vegf_pde->AddDiscreteSource(p_normal_and_quiescent_cell_source);

    // Add a vessel related VEGF sink
    auto p_vessel_vegf_sink = VesselBasedDiscreteSource<2>::Create();
    p_vessel_vegf_sink->SetReferenceConcentration(0.0_M);
    p_vessel_vegf_sink->SetVesselPermeability(Owen11Parameters::mpVesselVegfPermeability->GetValue("User"));
    // p_vessel_vegf_sink->SetVesselPermeability(0.0);
    // p_vegf_pde->AddDiscreteSource(p_vessel_vegf_sink);

    // Set up a finite difference solver as before for the VEGF
    auto p_vegf_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
    p_vegf_solver->SetPde(p_vegf_pde);
    p_vegf_solver->SetLabel("VEGF_Extracellular");  // set the label for the external field
    p_vegf_solver->SetGrid(p_grid);

    // Set up the angiogenesis solver (runs on deactivated VEGF function but seems to be necessary for main solver)
    std::shared_ptr<AngiogenesisSolver<2> > p_angiogenesis_solver = AngiogenesisSolver<2>::Create();
    std::shared_ptr<Owen2011SproutingRule<2> > p_sprouting_rule = Owen2011SproutingRule<2>::Create();
    std::shared_ptr<Owen2011MigrationRule<2> > p_migration_rule = Owen2011MigrationRule<2>::Create();
    p_angiogenesis_solver->SetMigrationRule(p_migration_rule);
    p_angiogenesis_solver->SetSproutingRule(p_sprouting_rule);
    p_sprouting_rule->SetDiscreteContinuumSolver(p_vegf_solver);
    p_migration_rule->SetDiscreteContinuumSolver(p_vegf_solver);
    p_angiogenesis_solver->SetVesselGridCalculator(p_grid_calc);
    p_angiogenesis_solver->SetVesselNetwork(p_network);
    
    // Initialise the solver
    auto p_microvessel_solver = MicrovesselSolver<2>::Create();
    p_microvessel_solver->SetVesselNetwork(p_network);
    p_microvessel_solver->SetOutputFrequency(5_h);
    // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
    p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
    p_microvessel_solver->AddDiscreteContinuumSolver(p_vegf_solver);
    p_microvessel_solver->SetAngiogenesisSolver(p_angiogenesis_solver);

    // Specifies which extracellular fields to update based on PDE
    boost::shared_ptr<MicrovesselSimulationModifier<2> > p_microvessel_modifier =
    boost::shared_ptr<MicrovesselSimulationModifier<2> >(new MicrovesselSimulationModifier<2> ());
    p_microvessel_modifier->SetMicrovesselSolver(p_microvessel_solver);
    std::vector<std::string> update_labels;
    update_labels.push_back("oxygen");
    p_microvessel_modifier->SetCellDataUpdateLabels(update_labels);

    // The full simulation is run as a typical Cell Based Chaste simulation
    OffLatticeSimulation<2> simulator(cell_population);
    simulator.AddSimulationModifier(p_microvessel_modifier);

    // Add a killer to remove apoptotic cells (leave it commented so RT is the only killer)
    // boost::shared_ptr<ApoptoticCellKiller<2> > p_apoptotic_cell_killer(new ApoptoticCellKiller<2>(p_cell_population.get())); (leave it commented so RT is the only killer)
    // simulator.AddCellKiller(p_apoptotic_cell_killer);

    // Add a LQ RT killer
    // boost::shared_ptr<LQRadiotherapyCellKiller<2> > p_rt_killer =
    // boost::shared_ptr<LQRadiotherapyCellKiller<2> >(new LQRadiotherapyCellKiller<2> (p_cell_population.get()));
    // p_rt_killer->SetAlphaMax(0.3*unit::per_gray);
    // p_rt_killer->SetBetaMax(0.03*unit::per_gray_squared);
    // p_rt_killer->SetDoseInjected(2.0*unit::gray);
    // p_rt_killer->SetCancerousRadiosensitivity(0.3 * unit::per_gray, 0.03 * unit::per_gray_squared);
    // p_rt_killer->SetNormalRadiosensitivity(0.15 * unit::per_gray, 0.05 * unit::per_gray_squared);
    // p_rt_killer->SetOerAlphaMax(1.75);
    // p_rt_killer->SetOerAlphaMin(1.0);
    // p_rt_killer->SetOerBetaMax(3.25);
    // p_rt_killer->SetOerBetaMin(1.0);
    // // p_rt_killer->SetOerConstant(oxygen_solubility_at_stp * 3.28_Pa);
    // p_rt_killer->SetOerConstant(0.004499758549541243*unit::mole_per_metre_cubed);  // K_OER (3.28 mmHg)
    // p_rt_killer->SetAlphaMax(0.3 * unit::per_gray);
    // p_rt_killer->SetBetaMax(0.03 * unit::per_gray_squared);
    // p_rt_killer->UseOer(true);

    // Set Up the dosage times
    // std::vector<QTime > rt_times;
    // rt_times.push_back(3600.0*24.0*unit::seconds);
    // rt_times.push_back(3600.0*48.0*unit::seconds);
    // rt_times.push_back(3600.0*72.0*unit::seconds);
    // p_rt_killer->SetTimeOfRadiation(rt_times);
    // p_rt_killer->AddTimeOfRadiation(3600.0*96.0*unit::seconds);
    // simulator.AddCellKiller(p_rt_killer);

    // If needed
    // boost::shared_ptr<RandomCellKiller<2> > p_cell_killer(new RandomCellKiller<2>(p_cell_population.get(), 0.01));
    // simulator.AddCellKiller(p_cell_killer);
    // for (typename AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population.get()->Begin(); cell_iter != p_cell_population.get()->End(); ++cell_iter)
    // {
    //     // CheckAndLabelSingleCellForApoptosis(*cell_iter);
    //     std::cout << "hi cell" << (*cell_iter)->GetCellData()->GetItem("oxygen") << std::endl;

    //     //  std::cout <<  BaseUnits::Instance()->GetReferenceConcentrationScale() << std::endl;


    //     //     double oxygen_concentration;

    //     // std::cout << "hi cell" << oxygen_concentration << std::endl;
    // }

    // Add another modifier for updating cell cycle quantities.
    // boost::shared_ptr<Owen2011TrackingModifier<2> > p_owen11_tracking_modifier(new Owen2011TrackingModifier<2>);
    // simulator.AddSimulationModifier(p_owen11_tracking_modifier);

    // Print the average oxygenation
    // std::vector<double> solution = p_oxygen_solver->GetSolution();
    // double average_oxygen = 0.0;
    // for(unsigned jdx=0;jdx<solution.size();jdx++)
    // {
    //     average_oxygen += solution[jdx];
    // }
    // average_oxygen /= double(solution.size());
    // std::cout << "Average oxygen: " << average_oxygen << std::endl;

    // Set up the simulation time and run it
    simulator.SetOutputDirectory(str_directory_name);
    // simulator.SetSamplingTimestepMultiple(10);  // get sample every x steps (multiple*dt = time between s)
    simulator.SetSamplingTimestepMultiple(1);
    // simulator.SetDt(1);  // simulation time step in hours
    simulator.SetDt(1.0);
    simulator.SetEndTime(1.0);  // end time in hours
    // simulator.SetEndTime(24);  // end time in hours
    // simulator.SetEndTime(2.0);
    simulator.Solve();

    // // Print the median oxygenation
    std::vector<double> solution = p_oxygen_solver->GetSolution();
    std::sort(solution.begin(), solution.end());

    double median_oxygen;
    size_t size = solution.size();
    if (size % 2 == 0) {
        // Even number of elements
        median_oxygen = (solution[size / 2 - 1] + solution[size / 2]) / 2.0;
    } else {
        // Odd number of elements
        median_oxygen = solution[size / 2];
    }
    std::cout << "Median oxygen: " << median_oxygen << std::endl;

    // Write the output network file (visualise with Paraview: set Filters->Alphabetical->Tube)
    std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
    p_network->Write(output_file);

    // Dump our parameter collection to an xml file and, importantly, clear it for the next test
    ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
    ParameterCollection::Instance()->Destroy();
    BaseUnits::Instance()->Destroy();
    SimulationTime::Instance()->Destroy();

    // Print the error log
    // std::string error_message = error_log.str();
    // std::cout << error_message << std::endl; 
  }


};

#endif /*TESTOXYGENCONSUMPTIONBUG_HPP_*/   