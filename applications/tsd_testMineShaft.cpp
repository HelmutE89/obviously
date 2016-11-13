#include <iostream>
#include "obgraphic/Obvious3D.h"
#include "obcore/base/tools.h"
#include "obcore/base/System.h"
#include "obvision/reconstruct/space/TsdSpace.h"
#include "obvision/reconstruct/space/SensorProjective3D.h"
#include "obvision/reconstruct/space/RayCast3D.h"
#include "obvision/reconstruct/space/RayCastAxisAligned3D.h"
#include "obvision/reconstruct/space/SensorPolar2DWith3DPose.h"
#include "obcore/base/Logger.h"

#include "obcore/math/mathbase.h"

using namespace std;
using namespace obvious;

Obvious3D* _viewer;
SensorPolar2DWith3DPose* _sensor;
TsdSpace* _space;
VtkCloud* _vcloud;

double NOISE_RANGE = 0.02;
double WALL_DISTANCE = 0.3;
int WOBBLE_RANGE_DEG = 5;
double Z_STEP_SIZE = 0.01 / 16;
double ANGULAR_RES_DEG = 0.125;


void generateSyntheticData(SensorPolar2DWith3DPose& sensor)
{

  int beams = sensor.getRealMeasurementSize();
  //double angularRes = sensor.getAngularResolution();
  //double minPhi = sensor.getPhiMin();
  double maxRange = sensor.getMaximumRange();
  //double minRange             = sensor.getMaximumRange();
  //double lowReflectivityRange = sensor.getLowReflectivityRange();

  // Sample data, to be replaced with real measurements
  double* data = new double[beams];
  unsigned char* rgb = new unsigned char[beams * 3];

  Matrix T = sensor.getTransformation();

  Matrix* rays = sensor.getNormalizedRayMap(1);


  for(int i = 0; i < beams; i++)
  {
    double n[2] = {0,0};

    //wall on x
    n[0] = std::abs(WALL_DISTANCE / (*rays)(0, i));

    //wall on y
    n[1] = std::abs(WALL_DISTANCE / (*rays)(1, i));

    //take nearest wall for distance
    double min = n[0] < n[1] ? n[0] : n[1];
    cout << "min1 " << min << endl;
    min += -NOISE_RANGE + ((double)rand() / RAND_MAX) * NOISE_RANGE * 2;
    cout << "min2 " << min << endl;

    //check maxRange
    data[i] = min <= maxRange ? min : NAN;

    //color it red
    rgb[i * 3] = 255;
  }
  sensor.setRealMeasurementData(data);
  sensor.setRealMeasurementRGB(rgb);
  sensor.setStandardMask();
}

void extractEulerAngleXYZ(Matrix t, double& rotXangle, double& rotYangle, double& rotZangle)
{
  rotXangle = atan2(-t(1, 2), t(2, 2));
  double cosYangle = sqrt(pow(t(0, 0), 2) + pow(t(0, 1), 2));
  rotYangle = atan2(t(0, 2), cosYangle);
  double sinXangle = sin(rotXangle);
  double cosXangle = cos(rotXangle);
  rotZangle = atan2(cosXangle * t(1, 0) + sinXangle * t(2, 0), cosXangle * t(1, 1) + sinXangle * t(2, 1));
}

void _cbRegNewImage(void)
{
  Matrix currentT = _sensor->getTransformation();
  double rotXangle, rotYangle, rotZangle;
  extractEulerAngleXYZ(currentT,rotXangle, rotYangle, rotZangle);

  rotYangle = rotYangle / M_PI * 180;
  rotYangle *= -1;
  rotYangle += -WOBBLE_RANGE_DEG + rand()%(WOBBLE_RANGE_DEG * 2 + 1);
  rotYangle = rotYangle * M_PI / 180;

  rotXangle = rotXangle / M_PI * 180;
  rotXangle *= -1;
  rotXangle += -WOBBLE_RANGE_DEG + rand()%(WOBBLE_RANGE_DEG * 2 + 1);
  rotXangle = rotXangle * M_PI / 180;


  double tf[16] = {
      cos(rotYangle),   sin(rotYangle) * sin(rotXangle),  sin(rotYangle) * cos(rotXangle),  0,
      0,                cos(rotXangle),                   -sin(rotXangle),                  0,
      -sin(rotYangle),  cos(rotYangle) * sin(rotXangle),  cos(rotYangle)*cos(rotXangle),    0,
      0,                0,                                0,                                1
  };

  Matrix T(4, 4);
  T.setData(tf);
  _sensor->transform(&T);

  currentT= _sensor->getTransformation();
  currentT(2,3) += Z_STEP_SIZE;
  _sensor->setTransformation(currentT);

  generateSyntheticData(*_sensor);
  _space->pushForward(_sensor);
  _viewer->showSensorPose(currentT);

  unsigned int cnt;

  unsigned int cells = _space->getXDimension() * _space->getYDimension() * _space->getZDimension();
  double* coords = new double[cells * 3];
  double* normals = new double[cells * 3];
  unsigned char* rgb = new unsigned char[cells * 3];
  RayCastAxisAligned3D raycaster;
  raycaster.calcCoords(_space, coords, NULL, NULL, &cnt);

  _vcloud->setCoords(coords, cnt / 3, 3, NULL);
  //_vcloud->setColors(rgb, cnt / 3 ,3);
  _viewer->update();

  delete[] coords;
  delete[] normals;
}

int main(void)
{
  LOGMSG_CONF("tsd_test.log", Logger::file_off | Logger::screen_on, DBG_DEBUG, DBG_DEBUG);

  obfloat voxelSize = 0.01;
  _space = new TsdSpace(voxelSize, LAYOUT_8x8x8, LAYOUT_256x256x256);
  _space->setMaxTruncation(3.0 * voxelSize);

  // translation of sensor
  obfloat tr[3];
  _space->getCentroid(tr);
  tr[2] = 0.1;

  // rotation about y-axis of sensor
  double theta = 0 * M_PI / 180;

  double tf[16] = {cos(theta), 0, sin(theta), tr[0], 0, 1, 0, tr[1], -sin(theta), 0, cos(theta), tr[2], 0, 0, 0, 1};
  Matrix T(4, 4);
  T.setData(tf);

  // Sensor initialization

  int beams = 360 / ANGULAR_RES_DEG;
  double angularResRad = deg2rad(ANGULAR_RES_DEG);
  double minPhi = deg2rad(-180.0);
  double maxRange = 3.0;
  double minRange = 0.3;
  double lowReflectivityRange = 0.5;

  _sensor = new SensorPolar2DWith3DPose(beams, angularResRad, minPhi, maxRange, minRange, lowReflectivityRange);

  _sensor->transform(&T);

  /*unsigned char* buffer = new unsigned char[space.getXDimension()*space.getYDimension()*3];
   for(unsigned int i=0; i<space.getZDimension(); i++)
   {
   char path[64];
   sprintf(path, "/tmp/slice%04d.ppm", i);
   _space->buildSliceImage(i, buffer);
   serializePPM(path, buffer, _space->getXDimension(), _space->getYDimension(), 0);
   }
   delete[] buffer;*/

  _vcloud = new VtkCloud;
  _viewer = new Obvious3D("TSD Cloud");

  _viewer->addAxisAlignedCube(0, _space->getMaxX(), 0, _space->getMaxY(), 0, _space->getMaxZ());
  _viewer->showAxes(true);
  _viewer->addCloud(_vcloud);
  _viewer->registerKeyboardCallback("space", _cbRegNewImage, "Register new image");
  generateSyntheticData(*_sensor);
  _viewer->startRendering();

  //delete[] coords;
  //delete[] normals;
  //if(rgb)
  //  delete[] rgb;
  delete _viewer;
  delete _sensor;
  delete _space;
}

