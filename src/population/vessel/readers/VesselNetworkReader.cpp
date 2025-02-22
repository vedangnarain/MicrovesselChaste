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

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkCleanPolyData.h>
#include <vtkSplineFilter.h>
#include <vtkVersion.h>
#include "Exception.hpp"
#include "VesselNetworkReader.hpp"

template<unsigned DIM>
VesselNetworkReader<DIM>::VesselNetworkReader()
    : mFileName(),
      mRadiusLabel("Node Radius"),
      mRadiusConversionFactor(1.0),
      mMergeCoincidentPoints(false),
      mTargetSegmentLength(0.0_m),
      mReferenceLength(1_um)
{

}

template<unsigned DIM>
VesselNetworkReader<DIM>::~VesselNetworkReader()
{
}

template <unsigned DIM>
std::shared_ptr<VesselNetworkReader<DIM> > VesselNetworkReader<DIM>::Create()
{
    return std::make_shared<VesselNetworkReader<DIM> >();

}

template <unsigned DIM>
void VesselNetworkReader<DIM>::SetRadiusArrayName(const std::string& rRadius)
{
    mRadiusLabel = rRadius;
}

template <unsigned DIM>
void VesselNetworkReader<DIM>::SetMergeCoincidentPoints(bool mergePoints)
{
	mMergeCoincidentPoints = mergePoints;
}

template <unsigned DIM>
void VesselNetworkReader<DIM>::SetTargetSegmentLength(QLength targetSegmentLength)
{
	mTargetSegmentLength = targetSegmentLength;
}

template <unsigned DIM>
void VesselNetworkReader<DIM>::SetReferenceLengthScale(QLength rReferenceLength)
{
    mReferenceLength = rReferenceLength;
}

template<unsigned DIM>
std::shared_ptr<VesselNetwork<DIM> > VesselNetworkReader<DIM>::Read()
{
    if(mFileName.empty())
    {
        EXCEPTION("File name not set in vessel network reader");
    }

    // Create an empty vessel network
    std::shared_ptr<VesselNetwork<DIM> > p_network = VesselNetwork<DIM>::Create();

    // Create a VTK PolyData object based on the contents of the input VTK file
    vtkSmartPointer<vtkXMLPolyDataReader> p_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    p_reader->SetFileName(mFileName.c_str());
    p_reader->Update();

    vtkSmartPointer<vtkPolyData> p_polydata = p_reader->GetOutput();

    if (mMergeCoincidentPoints) {
        vtkSmartPointer<vtkCleanPolyData> p_cleaner = vtkSmartPointer<
                vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
        p_cleaner->SetInput(p_polydata);
#else
        p_cleaner->SetInputData(p_polydata);
#endif
        p_cleaner->Update();
        p_polydata = p_cleaner->GetOutput();
    }

    if (mTargetSegmentLength != 0_m) {
        vtkSmartPointer<vtkSplineFilter> p_spline_filter = vtkSmartPointer<
                vtkSplineFilter>::New();
#if VTK_MAJOR_VERSION <= 5
        p_spline_filter->SetInput(p_polydata);
#else
        p_spline_filter->SetInputData(p_polydata);
#endif
        p_spline_filter->SetSubdivideToLength();
        p_spline_filter->SetLength(mTargetSegmentLength / mReferenceLength);
        p_spline_filter->Update();
        p_polydata = p_spline_filter->GetOutput();
    }

    // Create the nodes
    vtkSmartPointer<vtkPointData> p_point_data = vtkSmartPointer<vtkPointData>::New();
    std::vector<std::shared_ptr<VesselNode<DIM> > > nodes;
    for (vtkIdType i = 0; i < p_polydata->GetNumberOfPoints(); i++)
    {
        double point_coords[3];
        p_polydata->GetPoint(i, point_coords);
        nodes.push_back(VesselNode<DIM>::Create(Vertex<DIM>(point_coords, mReferenceLength)));
    }

    // Extract radii corresponding to each node from the VTK Polydata and store them in a list.
    p_point_data = p_polydata->GetPointData();
    std::vector<double> radii;
    for (vtkIdType i = 0; i < p_point_data->GetNumberOfArrays(); i++)
    {
        std::string array_name = p_point_data->GetArrayName(i);
        if (array_name.compare(mRadiusLabel) == 0)
        {
            for (vtkIdType j = 0; j < p_point_data->GetArray(i)->GetNumberOfTuples(); j++)
            {
                radii.push_back(p_point_data->GetArray(i)->GetTuple1(j));
            }
        }
    }

    // Extract vessels from the VTK Polydata. This is done by iterating over a VTK CellArray object which
    // returns a 'pointList' vtkIdList object. This object contains the point IDs of the nodes which make up
    // the vessel.
    vtkSmartPointer<vtkCellArray> pCellArray = vtkSmartPointer<vtkCellArray>::New();
    pCellArray = p_polydata->GetLines();
    for (int i = 0; i < p_polydata->GetNumberOfLines(); i++)
    {
        // Make a new vessel
        vtkIdType num_segments;
        vtkIdType* pSegmentList;
        pCellArray->GetNextCell(num_segments, pSegmentList);
        std::vector<std::shared_ptr<VesselSegment<DIM> > > segments;
        // Add segments to the vessels in order
        for (int j = 1; j < num_segments; j++)
        {
            std::shared_ptr<VesselSegment<DIM> > p_segment = VesselSegment<DIM>::Create(nodes[pSegmentList[j - 1]],nodes[pSegmentList[j]]);
            if(unsigned(radii.size())> pSegmentList[j])
            {
                p_segment->SetRadius(radii[pSegmentList[j]] * mReferenceLength);
            }
            segments.push_back(p_segment);
        }
        // Add the resulting vessel to the network
        p_network->AddVessel(Vessel<DIM>::Create(segments));
    }
    return p_network;
}

template<unsigned DIM>
void VesselNetworkReader<DIM>::SetFileName(const std::string& rFileName)
{
    mFileName = rFileName;
}

//Explicit instantiation
template class VesselNetworkReader<2> ;
template class VesselNetworkReader<3> ;
