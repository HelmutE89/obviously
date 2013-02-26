/*
 * Device3DFilter.cpp
 *
 *  Created on: 21.02.2013
 *      Author: phil
 */

#include <cstdlib>
#include <vector>
#include <cstring>
#include "Device3DFilter.h"
#include "obcore/math/mathbase.h"

namespace obvious
{

Device3DData::Device3DData(void)
{
	_pointCloud=NULL;
	_normals=NULL;
	_mask=NULL;
	_depthMap=NULL;
	_validPoints=0;
}

Device3DData::~Device3DData(void)
{
   if(_pointCloud)
      delete _pointCloud;
   if(_normals)
      delete _normals;
   if(_mask)
      delete _mask;
   if(_depthMap)
      delete _depthMap;
}

void Device3DData::setNormals(double* normals)
{
	//delete _normals;
	_normals=normals;
}

/****************************************************************************************************/

Device3DFilter::Device3DFilter(void)
{
   _input=NULL;
   _output=NULL;
   _threshold=0.0;
}

Device3DFilter::~Device3DFilter(void)
{
   if(_input)
      delete(_input);
   if(_output)
      delete(_output);
}

FILRETVAL Device3DFilter::applyFilter(void)
{
   if((!_input))//||(!_output))
   {
      std::cout<<"\nError! Invalid input pointer!\n";
      return(FILTER_ERROR);
   }
   if((_threshold<0.001)&&(_threshold>-0.001))
   {
      std::cout<<"\nError! Threshold not initialized!\n";
      return(FILTER_ERROR);
   }
   if(_output)
      delete _output;

   double* inPcl=_input->getPointCloud();
   double* inNorm=_input->getNormals();

   std::vector<double> outPcl;
   std::vector<double> outNorm;
   double* outDepthmap=new double[640*480];
   bool* outMask=new bool[640*480];
   double abs=0.0;

   for(unsigned int i=0;i<640*480;i++)
   {
      abs =  euklideanDistance<double>(inPcl, NULL, 3);//Point3D::sizeP);
      if(abs>_threshold)
      {
         outMask[i]=false;                  //set mask to false
         outDepthmap[i]=0.0;                //set DepthMap to NaN
         inPcl+=3;
         inNorm+=3;
         continue;
      }
      else
      {
         for(unsigned int j=0;j<3;j++)
         {
            outPcl.push_back(*inPcl);
            inPcl++;
            outNorm.push_back(*inNorm);
            inNorm++;
         }
         outMask[i]=_input->getMask()[i];
         outDepthmap[i]=_input->getDepthMap()[i];
      }
   }
   if(outPcl.size()!=outNorm.size())
   {
   	std::cout<<"\nError! Nbr. of points differs nbr. of normals\n";
   	return(FILTER_ERROR);
   }

   //Generate filtered data object
   _output=new Device3DData();
   double* outputPcl=new double[outPcl.size()];
   double* outputNorm=new double[outNorm.size()];
   memcpy(outputPcl,outPcl.data(),outPcl.size()*sizeof(double));
   memcpy(outputNorm,outNorm.data(),outNorm.size()*sizeof(double));
   _output->setPointCloud(outputPcl);
   _output->setNormals(outputNorm);
   _output->setDepthMap(outDepthmap);
   _output->setMask(outMask);
   _output->setValidPoints(outPcl.size()/3);

   return(FILTER_OK);
}

Device3DData* Device3DFilter::getOutput(void)
{
	if(!_output)
	{
		std::cout<<"\nError! Output is empty!\n";
		return(NULL);
	}
	return(_output);
}


}